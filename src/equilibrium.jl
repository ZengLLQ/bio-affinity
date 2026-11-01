"""
    equilibrium.jl

Equilibrium saturation binding, Scatchard transformation, Cheng-Prusoff IC50-to-Ki
transformations, and Morrison tight-binding quadratic formulations.
"""

module Equilibrium

using LinearAlgebra

calc_mean(v) = sum(v) / length(v)

export SaturationBinding, MorrisonParams
export saturation_fraction, scatchard_transform, linear_scatchard_fit
export cheng_prusoff_competitive, cheng_prusoff_uncompetitive, cheng_prusoff_noncompetitive, cheng_prusoff_mixed
export morrison_fractional_velocity, morrison_ic50_apparent, exact_free_ligand

"""
    SaturationBinding

Model parameters for equilibrium saturation:
- `bmax`: Maximum binding capacity Bmax (M or fmol/mg)
- `kd`: Equilibrium dissociation constant Kd (M)
- `hill`: Hill cooperativity coefficient n_H (default 1.0)
"""
struct SaturationBinding
    bmax::Float64
    kd::Float64
    hill::Float64

    function SaturationBinding(bmax::Real, kd::Real; hill::Real=1.0)
        bmax > 0 || throw(ArgumentError("bmax must be positive"))
        kd > 0 || throw(ArgumentError("kd must be positive"))
        hill > 0 || throw(ArgumentError("hill coefficient must be positive"))
        new(Float64(bmax), Float64(kd), Float64(hill))
    end
end

"""
    saturation_fraction(p::SaturationBinding, ligand::Real) -> Float64

Compute specific bound ligand [RL] given free ligand concentration [L]:
    [RL] = Bmax * [L]^nH / (Kd^nH + [L]^nH)
"""
function saturation_fraction(p::SaturationBinding, ligand::Real)
    ligand >= 0 || throw(ArgumentError("ligand concentration must be non-negative"))
    if ligand == 0.0
        return 0.0
    end
    l_h = Float64(ligand)^p.hill
    kd_h = p.kd^p.hill
    return p.bmax * l_h / (kd_h + l_h)
end

"""
    scatchard_transform(bound::AbstractVector{<:Real}, free::AbstractVector{<:Real}) -> Tuple{Vector{Float64}, Vector{Float64}}

Compute Scatchard coordinates (x = Bound, y = Bound / Free).
"""
function scatchard_transform(bound::AbstractVector{<:Real}, free::AbstractVector{<:Real})
    length(bound) == length(free) || throw(DimensionMismatch("bound and free must have equal lengths"))
    x = Float64.(bound)
    y = zeros(Float64, length(bound))
    for i in eachindex(bound)
        free[i] > 0 || throw(ArgumentError("free ligand must be strictly positive for ratio"))
        y[i] = bound[i] / free[i]
    end
    return (x, y)
end

"""
    linear_scatchard_fit(bound::AbstractVector{<:Real}, free::AbstractVector{<:Real}) -> Tuple{Float64, Float64, Float64}

Fit Scatchard line via linear regression y = m*x + c.
Returns `(Kd_est, Bmax_est, r_squared)` where slope m = -1/Kd and x-intercept c/(-m) = Bmax.
"""
function linear_scatchard_fit(bound::AbstractVector{<:Real}, free::AbstractVector{<:Real})
    x, y = scatchard_transform(bound, free)
    n = length(x)
    n >= 2 || throw(ArgumentError("At least 2 points required for linear fit"))

    x_mean = calc_mean(x)
    y_mean = calc_mean(y)

    ss_xx = sum((xi - x_mean)^2 for xi in x)
    ss_xy = sum((x[i] - x_mean) * (y[i] - y_mean) for i in 1:n)
    ss_yy = sum((yi - y_mean)^2 for yi in y)

    ss_xx > 0 || throw(ArgumentError("Degenerate data: zero variance in bound ligand"))

    slope = ss_xy / ss_xx
    intercept = y_mean - slope * x_mean

    slope < 0 || throw(ArgumentError("Scatchard slope must be negative for valid single-site binding"))

    kd_est = -1.0 / slope
    bmax_est = intercept / (-slope)

    r_sq = (ss_xy^2) / (ss_xx * ss_yy)

    return (kd_est, bmax_est, r_sq)
end

"""
    cheng_prusoff_competitive(ic50::Real, ligand::Real, kd::Real) -> Float64

Calculate Ki from IC50 for a competitive inhibitor:
    Ki = IC50 / (1 + [L] / Kd)
"""
function cheng_prusoff_competitive(ic50::Real, ligand::Real, kd::Real)
    ic50 > 0 || throw(ArgumentError("ic50 must be positive"))
    ligand >= 0 || throw(ArgumentError("ligand concentration must be non-negative"))
    kd > 0 || throw(ArgumentError("kd must be positive"))
    return Float64(ic50) / (1.0 + Float64(ligand) / Float64(kd))
end

"""
    cheng_prusoff_uncompetitive(ic50::Real, ligand::Real, kd::Real) -> Float64

Calculate Ki from IC50 for an uncompetitive inhibitor:
    Ki = IC50 / (1 + Kd / [L])
"""
function cheng_prusoff_uncompetitive(ic50::Real, ligand::Real, kd::Real)
    ic50 > 0 || throw(ArgumentError("ic50 must be positive"))
    ligand > 0 || throw(ArgumentError("ligand concentration must be positive for uncompetitive"))
    kd > 0 || throw(ArgumentError("kd must be positive"))
    return Float64(ic50) / (1.0 + Float64(kd) / Float64(ligand))
end

"""
    cheng_prusoff_noncompetitive(ic50::Real) -> Float64

Calculate Ki from IC50 for a pure non-competitive inhibitor (Ki = IC50).
"""
function cheng_prusoff_noncompetitive(ic50::Real)
    ic50 > 0 || throw(ArgumentError("ic50 must be positive"))
    return Float64(ic50)
end

"""
    cheng_prusoff_mixed(ic50::Real, substrate::Real, km::Real, kic::Real, kiu::Real) -> Float64

General relationship for mixed competitive-uncompetitive enzyme inhibition.
Returns expected IC50 given substrate [S], Km, Kic (competitive), and Kiu (uncompetitive):
    IC50 = ([S] + Km) / (Km / Kic + [S] / Kiu)
"""
function cheng_prusoff_mixed(substrate::Real, km::Real, kic::Real, kiu::Real)
    substrate >= 0 || throw(ArgumentError("substrate must be non-negative"))
    km > 0 || throw(ArgumentError("km must be positive"))
    kic > 0 || throw(ArgumentError("kic must be positive"))
    kiu > 0 || throw(ArgumentError("kiu must be positive"))
    num = Float64(substrate) + Float64(km)
    denom = (Float64(km) / Float64(kic)) + (Float64(substrate) / Float64(kiu))
    return num / denom
end

"""
    MorrisonParams

Parameters for tight-binding Morrison equation:
- `etot`: Total enzyme or receptor concentration [E_tot] (M)
- `ki`: True inhibitor inhibition constant Ki (M)
- `substrate`: Substrate or tracer concentration [S] (M)
- `km`: Michaelis-Menten constant or tracer Kd (M)
"""
struct MorrisonParams
    etot::Float64
    ki::Float64
    substrate::Float64
    km::Float64

    function MorrisonParams(etot::Real, ki::Real; substrate::Real=0.0, km::Real=1.0)
        etot > 0 || throw(ArgumentError("etot must be positive"))
        ki > 0 || throw(ArgumentError("ki must be positive"))
        substrate >= 0 || throw(ArgumentError("substrate must be non-negative"))
        km > 0 || throw(ArgumentError("km must be positive"))
        new(Float64(etot), Float64(ki), Float64(substrate), Float64(km))
    end
end

"""
    morrison_ic50_apparent(p::MorrisonParams) -> Float64

Calculate apparent tight-binding IC50:
    IC50 = Ki_app + 0.5 * [E_tot]
where Ki_app = Ki * (1 + [S] / Km).
"""
function morrison_ic50_apparent(p::MorrisonParams)
    kiapp = p.ki * (1.0 + p.substrate / p.km)
    return kiapp + 0.5 * p.etot
end

"""
    morrison_fractional_velocity(p::MorrisonParams, itot::Real) -> Float64

Compute fractional enzyme activity or receptor binding v / v0 using the quadratic Morrison equation:
    v / v0 = 1 - [ (Etot + Itot + Ki_app) - sqrt((Etot + Itot + Ki_app)^2 - 4 * Etot * Itot) ] / (2 * Etot)
"""
function morrison_fractional_velocity(p::MorrisonParams, itot::Real)
    itot >= 0 || throw(ArgumentError("itot must be non-negative"))
    if itot == 0.0
        return 1.0
    end

    kiapp = p.ki * (1.0 + p.substrate / p.km)
    sum_term = p.etot + Float64(itot) + kiapp
    disc = sum_term^2 - 4.0 * p.etot * Float64(itot)
    root_disc = sqrt(max(disc, 0.0))

    bound_fraction = (sum_term - root_disc) / (2.0 * p.etot)
    return clamp(1.0 - bound_fraction, 0.0, 1.0)
end

"""
    exact_free_ligand(ltot::Real, rtot::Real, kd::Real) -> Float64

Calculate exact free ligand concentration [L] accounting for ligand depletion:
Solves quadratic: [L]^2 + (Kd + [R_tot] - [L_tot]) * [L] - Kd * [L_tot] = 0.
"""
function exact_free_ligand(ltot::Real, rtot::Real, kd::Real)
    ltot >= 0 || throw(ArgumentError("ltot must be non-negative"))
    rtot >= 0 || throw(ArgumentError("rtot must be non-negative"))
    kd > 0 || throw(ArgumentError("kd must be positive"))

    if ltot == 0.0
        return 0.0
    end

    b = kd + rtot - ltot
    c = -kd * ltot
    disc = b^2 - 4.0 * c
    # Return positive physical root
    return (-b + sqrt(disc)) / 2.0
end

end # module Equilibrium
