"""
    allostery.jl

Ehlert & Hall Allosteric Ternary Complex Model (ATCM), binding cooperativity (alpha),
signaling efficacy cooperativity (beta), thermodynamic cycle consistency, and the
Monod-Wyman-Changeux (MWC) concerted transition model.
"""

module Allostery

using LinearAlgebra

export AllostericTernaryParams, MWCParams
export allosteric_orthosteric_occupancy, allosteric_apparent_kd, allosteric_functional_response
export thermodynamic_cycle_free_energy, mwc_fractional_saturation, mwc_apparent_hill

"""
    AllostericTernaryParams

Parameters for Ehlert & Hall Allosteric Ternary Complex Model:
- `ka`: Orthosteric ligand equilibrium dissociation constant KA (M)
- `kb`: Allosteric modulator equilibrium dissociation constant KB (M)
- `alpha`: Binding cooperativity factor (α > 1 PAM, α < 1 NAM, α = 1 NAL)
- `beta`: Signaling efficacy cooperativity factor (β > 1 positive, β < 1 negative)
"""
struct AllostericTernaryParams
    ka::Float64
    kb::Float64
    alpha::Float64
    beta::Float64

    function AllostericTernaryParams(ka::Real, kb::Real; alpha::Real=1.0, beta::Real=1.0)
        ka > 0 || throw(ArgumentError("ka must be positive"))
        kb > 0 || throw(ArgumentError("kb must be positive"))
        alpha >= 0 || throw(ArgumentError("alpha must be non-negative"))
        beta >= 0 || throw(ArgumentError("beta must be non-negative"))
        new(Float64(ka), Float64(kb), Float64(alpha), Float64(beta))
    end
end

"""
    allosteric_orthosteric_occupancy(p::AllostericTernaryParams, a::Real, b::Real) -> Float64

Fraction of receptors with orthosteric site occupied (both AR and ARB species):
    YA = ([A]/KA * (1 + alpha * [B]/KB)) / (1 + [A]/KA + [B]/KB + alpha * [A]*[B]/(KA*KB))
"""
function allosteric_orthosteric_occupancy(p::AllostericTernaryParams, a::Real, b::Real)
    a >= 0 || throw(ArgumentError("orthosteric concentration [A] must be non-negative"))
    b >= 0 || throw(ArgumentError("allosteric concentration [B] must be non-negative"))

    if a == 0.0
        return 0.0
    end

    norm_a = Float64(a) / p.ka
    norm_b = Float64(b) / p.kb

    num = norm_a * (1.0 + p.alpha * norm_b)
    denom = 1.0 + norm_a + norm_b + p.alpha * norm_a * norm_b

    return num / denom
end

"""
    allosteric_apparent_kd(p::AllostericTernaryParams, b::Real) -> Float64

Calculate apparent orthosteric dissociation constant in presence of modulator [B]:
    KA_app = KA * (1 + [B]/KB) / (1 + alpha * [B]/KB)
"""
function allosteric_apparent_kd(p::AllostericTernaryParams, b::Real)
    b >= 0 || throw(ArgumentError("allosteric concentration [B] must be non-negative"))
    norm_b = Float64(b) / p.kb
    return p.ka * (1.0 + norm_b) / (1.0 + p.alpha * norm_b)
end

"""
    allosteric_functional_response(p::AllostericTernaryParams, a::Real, b::Real) -> Float64

Normalized functional response where species AR has efficacy 1.0 and ternary complex ARB
has efficacy beta:
    Response = ([A]/KA + alpha * beta * [A]*[B]/(KA*KB)) / (1 + [A]/KA + [B]/KB + alpha * [A]*[B]/(KA*KB))
"""
function allosteric_functional_response(p::AllostericTernaryParams, a::Real, b::Real)
    a >= 0 || throw(ArgumentError("orthosteric concentration [A] must be non-negative"))
    b >= 0 || throw(ArgumentError("allosteric concentration [B] must be non-negative"))

    if a == 0.0
        return 0.0
    end

    norm_a = Float64(a) / p.ka
    norm_b = Float64(b) / p.kb

    num = norm_a + p.alpha * p.beta * norm_a * norm_b
    denom = 1.0 + norm_a + norm_b + p.alpha * norm_a * norm_b

    return num / denom
end

"""
    thermodynamic_cycle_free_energy(p::AllostericTernaryParams; temp_k::Real=310.15) -> NamedTuple

Verify path-independence of allosteric ternary cycle free energy:
    ΔG° = -R * T * ln(K_eq)
Returns `(delta_G_path1, delta_G_path2, cooperativity_delta_G)`.
"""
function thermodynamic_cycle_free_energy(p::AllostericTernaryParams; temp_k::Real=310.15)
    temp_k > 0 || throw(ArgumentError("temperature must be positive in Kelvin"))
    R = 8.31446261815324 # J / (mol * K)

    # Path 1: R -> AR -> ARB
    # ΔG_A = -R*T*ln(1/KA) = R*T*ln(KA)
    # ΔG_B_given_A = R*T*ln(KB / alpha)
    dg1 = R * temp_k * (log(p.ka) + log(p.kb / p.alpha))

    # Path 2: R -> RB -> ARB
    # ΔG_B = R*T*ln(KB)
    # ΔG_A_given_B = R*T*ln(KA / alpha)
    dg2 = R * temp_k * (log(p.kb) + log(p.ka / p.alpha))

    # ΔΔG_alpha: coupling free energy = -R*T*ln(alpha)
    ddg_alpha = -R * temp_k * log(p.alpha)

    return (delta_G_path1=dg1, delta_G_path2=dg2, cooperativity_delta_G=ddg_alpha)
end

"""
    MWCParams

Parameters for Monod-Wyman-Changeux (MWC) concerted transition allosteric model:
- `n_subunits`: Number of interacting allosteric subunits (e.g., 4 for tetramer)
- `l0`: Allosteric equilibrium constant in unliganded state L0 = [T0] / [R0]
- `kr`: Microscopic dissociation constant of Relaxed (R) state (M)
- `kt`: Microscopic dissociation constant of Tense (T) state (M)
"""
struct MWCParams
    n_subunits::Int
    l0::Float64
    kr::Float64
    kt::Float64

    function MWCParams(n_subunits::Integer, l0::Real, kr::Real, kt::Real)
        n_subunits >= 1 || throw(ArgumentError("n_subunits must be at least 1"))
        l0 > 0 || throw(ArgumentError("l0 must be positive"))
        kr > 0 || throw(ArgumentError("kr must be positive"))
        kt > 0 || throw(ArgumentError("kt must be positive"))
        new(Int(n_subunits), Float64(l0), Float64(kr), Float64(kt))
    end
end

"""
    mwc_fractional_saturation(p::MWCParams, ligand::Real) -> Float64

Fractional saturation Y_bar of the MWC allosteric oligomer:
    Y_bar = [ alpha_L * (1 + alpha_L)^(n-1) + L0 * c * alpha_L * (1 + c * alpha_L)^(n-1) ] /
            [ (1 + alpha_L)^n + L0 * (1 + c * alpha_L)^n ]
where alpha_L = [L] / KR and c = KR / KT.
"""
function mwc_fractional_saturation(p::MWCParams, ligand::Real)
    ligand >= 0 || throw(ArgumentError("ligand concentration must be non-negative"))
    if ligand == 0.0
        return 0.0
    end

    alpha_l = Float64(ligand) / p.kr
    c = p.kr / p.kt
    n = p.n_subunits

    term_r = alpha_l * (1.0 + alpha_l)^(n - 1)
    term_t = p.l0 * c * alpha_l * (1.0 + c * alpha_l)^(n - 1)

    denom_r = (1.0 + alpha_l)^n
    denom_t = p.l0 * (1.0 + c * alpha_l)^n

    return (term_r + term_t) / (denom_r + denom_t)
end

"""
    mwc_apparent_hill(p::MWCParams; ligand_grid_points::Int=500) -> Float64

Estimate maximum apparent Hill coefficient nH for MWC saturation curve:
    nH = max d(ln(Y / (1 - Y))) / d(ln [L])
"""
function mwc_apparent_hill(p::MWCParams; ligand_grid_points::Int=500)
    # Search around midpoint Kd range
    mid_l = sqrt(p.kr * p.kt)
    log_min = log10(mid_l) - 3.0
    log_max = log10(mid_l) + 3.0
    log_ligands = range(log_min, stop=log_max, length=ligand_grid_points)

    max_nh = 1.0

    for i in 2:(length(log_ligands)-1)
        l_prev = 10.0^(log_ligands[i-1])
        l_next = 10.0^(log_ligands[i+1])

        y_prev = mwc_fractional_saturation(p, l_prev)
        y_next = mwc_fractional_saturation(p, l_next)

        if 0.05 < y_prev < 0.95 && 0.05 < y_next < 0.95
            logit_prev = log(y_prev / (1.0 - y_prev))
            logit_next = log(y_next / (1.0 - y_next))
            d_log_l = log(l_next) - log(l_prev)
            local_nh = (logit_next - logit_prev) / d_log_l
            if local_nh > max_nh
                max_nh = local_nh
            end
        end
    end

    return max_nh
end

end # module Allostery
