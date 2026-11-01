"""
    kinetics.jl

Core binding kinetics, association/dissociation time-courses, drug-target residence
time estimation, and the Motulsky & Mahan model of competitive association kinetics.
"""

module Kinetics

using LinearAlgebra
using Printf

export BindingKinetics1Order, MotulskyMahanParams, KineticsResult
export association_1st_order, dissociation_1st_order, residence_time, half_life
export motulsky_mahan_analytical, motulsky_mahan_ode, simulate_competitive_kinetics

"""
    BindingKinetics1Order

Parameters for pseudo-first order single ligand binding:
- `kon`: association rate constant (M^-1 s^-1)
- `koff`: dissociation rate constant (s^-1)
- `rtot`: total receptor concentration (M)
- `ligand`: ligand concentration (M)
"""
struct BindingKinetics1Order
    kon::Float64
    koff::Float64
    rtot::Float64
    ligand::Float64

    function BindingKinetics1Order(kon::Real, koff::Real, rtot::Real, ligand::Real)
        kon > 0 || throw(ArgumentError("kon must be positive"))
        koff > 0 || throw(ArgumentError("koff must be positive"))
        rtot > 0 || throw(ArgumentError("rtot must be positive"))
        ligand >= 0 || throw(ArgumentError("ligand must be non-negative"))
        new(Float64(kon), Float64(koff), Float64(rtot), Float64(ligand))
    end
end

"""
    residence_time(koff::Real) -> Float64

Calculate drug-target residence time τ = 1 / koff in seconds.
"""
function residence_time(koff::Real)
    koff > 0 || throw(ArgumentError("koff must be positive"))
    return 1.0 / Float64(koff)
end

"""
    half_life(koff::Real) -> Float64

Calculate receptor-ligand complex dissociation half-life t_1/2 = ln(2) / koff in seconds.
"""
function half_life(koff::Real)
    koff > 0 || throw(ArgumentError("koff must be positive"))
    return log(2.0) / Float64(koff)
end

"""
    association_1st_order(p::BindingKinetics1Order, t::Real) -> Float64

Compute bound complex [RL](t) under pseudo-first order conditions:
    [RL](t) = (kon * [L] * [R_tot] / k_obs) * (1 - exp(-k_obs * t))
where k_obs = kon * [L] + koff.
"""
function association_1st_order(p::BindingKinetics1Order, t::Real)
    t >= 0 || throw(ArgumentError("t must be non-negative"))
    kobs = p.kon * p.ligand + p.koff
    req = (p.kon * p.ligand * p.rtot) / kobs
    return req * (1.0 - exp(-kobs * Float64(t)))
end

"""
    dissociation_1st_order(rl0::Real, koff::Real, t::Real) -> Float64

Compute remaining bound complex [RL](t) = [RL]_0 * exp(-koff * t).
"""
function dissociation_1st_order(rl0::Real, koff::Real, t::Real)
    t >= 0 || throw(ArgumentError("t must be non-negative"))
    koff > 0 || throw(ArgumentError("koff must be positive"))
    return Float64(rl0) * exp(-Float64(koff) * Float64(t))
end

"""
    MotulskyMahanParams

Parameters for Motulsky & Mahan competitive association kinetics:
- `k1`: association rate of labeled tracer (M^-1 s^-1)
- `k2`: dissociation rate of labeled tracer (s^-1)
- `k3`: association rate of unlabeled competitor (M^-1 s^-1)
- `k4`: dissociation rate of unlabeled competitor (s^-1)
- `rtot`: total receptor concentration Bmax (M)
- `tracer`: tracer concentration [L] (M)
- `competitor`: competitor concentration [I] (M)
"""
struct MotulskyMahanParams
    k1::Float64
    k2::Float64
    k3::Float64
    k4::Float64
    rtot::Float64
    tracer::Float64
    competitor::Float64

    function MotulskyMahanParams(k1::Real, k2::Real, k3::Real, k4::Real,
                                rtot::Real, tracer::Real, competitor::Real)
        k1 > 0 || throw(ArgumentError("k1 must be positive"))
        k2 > 0 || throw(ArgumentError("k2 must be positive"))
        k3 > 0 || throw(ArgumentError("k3 must be positive"))
        k4 > 0 || throw(ArgumentError("k4 must be positive"))
        rtot > 0 || throw(ArgumentError("rtot must be positive"))
        tracer >= 0 || throw(ArgumentError("tracer must be non-negative"))
        competitor >= 0 || throw(ArgumentError("competitor must be non-negative"))
        new(Float64(k1), Float64(k2), Float64(k3), Float64(k4),
            Float64(rtot), Float64(tracer), Float64(competitor))
    end
end

"""
    motulsky_mahan_analytical(p::MotulskyMahanParams, t::Real) -> Tuple{Float64, Float64}

Analytical closed-form solution of Motulsky & Mahan competitive kinetics.
Returns `(RL_bound, RI_bound)` at time `t`.

Formulation:
    KA = k1 * [L] + k2
    KB = k3 * [I] + k4
    KF = 0.5 * ((KA + KB) + sqrt((KA - KB)^2 + 4 * k1 * [L] * k3 * [I]))
    KS = 0.5 * ((KA + KB) - sqrt((KA - KB)^2 + 4 * k1 * [L] * k3 * [I]))

    RL(t) = [R_tot * k1 * [L] / (KF - KS)] * [ ((k4 - KS)/KS)*(1 - exp(-KS*t)) - ((k4 - KF)/KF)*(1 - exp(-KF*t)) ]
    RI(t) = [R_tot * k3 * [I] / (KF - KS)] * [ ((k2 - KS)/KS)*(1 - exp(-KS*t)) - ((k2 - KF)/KF)*(1 - exp(-KF*t)) ]
"""
function motulsky_mahan_analytical(p::MotulskyMahanParams, t::Real)
    t >= 0 || throw(ArgumentError("t must be non-negative"))
    if t == 0.0
        return (0.0, 0.0)
    end

    KA = p.k1 * p.tracer + p.k2
    KB = p.k3 * p.competitor + p.k4

    disc = (KA - KB)^2 + 4.0 * p.k1 * p.tracer * p.k3 * p.competitor
    root_disc = sqrt(max(disc, 0.0))

    KF = 0.5 * ((KA + KB) + root_disc)
    KS = 0.5 * ((KA + KB) - root_disc)

    denom = KF - KS
    if denom < 1e-15
        # Degenerate case when rates coincide
        denom = 1e-15
    end

    # Tracer binding [RL]
    term_s_rl = ((p.k4 - KS) / KS) * (1.0 - exp(-KS * Float64(t)))
    term_f_rl = ((p.k4 - KF) / KF) * (1.0 - exp(-KF * Float64(t)))
    rl = (p.rtot * p.k1 * p.tracer / denom) * (term_s_rl - term_f_rl)

    # Competitor binding [RI]
    term_s_ri = ((p.k2 - KS) / KS) * (1.0 - exp(-KS * Float64(t)))
    term_f_ri = ((p.k2 - KF) / KF) * (1.0 - exp(-KF * Float64(t)))
    ri = (p.rtot * p.k3 * p.competitor / denom) * (term_s_ri - term_f_ri)

    return (max(0.0, rl), max(0.0, ri))
end

"""
    motulsky_mahan_ode(p::MotulskyMahanParams, t_span::Real; dt::Real=0.1, depletion::Bool=false)

Numerical 4th-order Runge-Kutta (RK4) integration for competitive binding.
If `depletion=true`, accounts for receptor depletion and ligand consumption.
Returns `(times, rl_series, ri_series)`.
"""
function motulsky_mahan_ode(p::MotulskyMahanParams, t_span::Real; dt::Real=0.1, depletion::Bool=false)
    t_span > 0 || throw(ArgumentError("t_span must be positive"))
    dt > 0 || throw(ArgumentError("dt must be positive"))

    n_steps = ceil(Int, t_span / dt)
    actual_dt = t_span / n_steps

    times = collect(range(0.0, stop=t_span, length=n_steps+1))
    rl_arr = zeros(Float64, n_steps + 1)
    ri_arr = zeros(Float64, n_steps + 1)

    rl = 0.0
    ri = 0.0

    function rhs(y_rl, y_ri)
        if depletion
            r_free = max(0.0, p.rtot - y_rl - y_ri)
            l_free = max(0.0, p.tracer - y_rl)
            i_free = max(0.0, p.competitor - y_ri)
            drl = p.k1 * l_free * r_free - p.k2 * y_rl
            dri = p.k3 * i_free * r_free - p.k4 * y_ri
        else
            r_free = max(0.0, p.rtot - y_rl - y_ri)
            drl = p.k1 * p.tracer * r_free - p.k2 * y_rl
            dri = p.k3 * p.competitor * r_free - p.k4 * y_ri
        end
        return (drl, dri)
    end

    for i in 1:n_steps
        # RK4 step
        k1_rl, k1_ri = rhs(rl, ri)

        k2_rl, k2_ri = rhs(rl + 0.5 * actual_dt * k1_rl, ri + 0.5 * actual_dt * k1_ri)

        k3_rl, k3_ri = rhs(rl + 0.5 * actual_dt * k2_rl, ri + 0.5 * actual_dt * k2_ri)

        k4_rl, k4_ri = rhs(rl + actual_dt * k3_rl, ri + actual_dt * k3_ri)

        rl += (actual_dt / 6.0) * (k1_rl + 2.0 * k2_rl + 2.0 * k3_rl + k4_rl)
        ri += (actual_dt / 6.0) * (k1_ri + 2.0 * k2_ri + 2.0 * k3_ri + k4_ri)

        rl_arr[i+1] = max(0.0, rl)
        ri_arr[i+1] = max(0.0, ri)
    end

    return (times, rl_arr, ri_arr)
end

"""
    KineticsResult

Container for simulated or fitted binding kinetics diagnostics.
"""
struct KineticsResult
    times::Vector{Float64}
    bound_tracer::Vector{Float64}
    bound_competitor::Vector{Float64}
    residence_time_tracer::Float64
    residence_time_competitor::Float64
    peak_tracer_overshoot::Float64
    equilibrium_fraction_tracer::Float64
end

"""
    simulate_competitive_kinetics(p::MotulskyMahanParams, t_span::Real; n_points::Int=200) -> KineticsResult

Generate complete time-series and diagnostic metrics for competitive association.
"""
function simulate_competitive_kinetics(p::MotulskyMahanParams, t_span::Real; n_points::Int=200)
    times = collect(range(0.0, stop=t_span, length=n_points))
    rl = zeros(Float64, n_points)
    ri = zeros(Float64, n_points)

    for (idx, t) in enumerate(times)
        val_rl, val_ri = motulsky_mahan_analytical(p, t)
        rl[idx] = val_rl
        ri[idx] = val_ri
    end

    tau_tracer = residence_time(p.k2)
    tau_comp = residence_time(p.k4)

    # Overshoot peak detection
    peak_rl = maximum(rl)
    eq_rl = rl[end]

    return KineticsResult(times, rl, ri, tau_tracer, tau_comp, peak_rl, eq_rl / p.rtot)
end

end # module Kinetics
