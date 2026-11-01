"""
    simulator.jl

Concentration-response sigmoidal curve simulation, 4-Parameter Logistic (4PL) modeling,
and Schild regression analysis for competitive receptor antagonists.
"""

module Simulator

using LinearAlgebra

calc_mean(v) = sum(v) / length(v)

export SigmoidCurve, SchildResult
export generate_sigmoidal_curve, evaluate_4pl, fit_sigmoidal_4pl, perform_schild_analysis

"""
    SigmoidCurve

Parameters for 4-Parameter Logistic (4PL) sigmoidal curve:
- `bottom`: Baseline response
- `top`: Plateau maximal response
- `log_ec50`: Log10 of EC50 (log10(M))
- `hill`: Hill slope
"""
struct SigmoidCurve
    bottom::Float64
    top::Float64
    log_ec50::Float64
    hill::Float64

    function SigmoidCurve(bottom::Real, top::Real, log_ec50::Real; hill::Real=1.0)
        hill != 0 || throw(ArgumentError("hill slope cannot be zero"))
        new(Float64(bottom), Float64(top), Float64(log_ec50), Float64(hill))
    end
end

"""
    evaluate_4pl(p::SigmoidCurve, conc::Real) -> Float64

Evaluate 4PL response at molar concentration `conc`:
    Response = Bottom + (Top - Bottom) / (1 + 10^((log_ec50 - log10(conc)) * hill))
"""
function evaluate_4pl(p::SigmoidCurve, conc::Real)
    conc >= 0 || throw(ArgumentError("concentration must be non-negative"))
    if conc == 0.0
        return p.bottom
    end

    log_c = log10(Float64(conc))
    exp_term = 10.0^((p.log_ec50 - log_c) * p.hill)
    return p.bottom + (p.top - p.bottom) / (1.0 + exp_term)
end

"""
    generate_sigmoidal_curve(p::SigmoidCurve; log_min::Real=-12.0, log_max::Real=-3.0, n_points::Int=100)

Generate (concentrations, responses) array across logarithmic concentration spectrum.
"""
function generate_sigmoidal_curve(p::SigmoidCurve; log_min::Real=-12.0, log_max::Real=-3.0, n_points::Int=100)
    log_concs = range(Float64(log_min), stop=Float64(log_max), length=n_points)
    concs = [10.0^lc for lc in log_concs]
    responses = [evaluate_4pl(p, c) for c in concs]
    return (concs, responses)
end

"""
    fit_sigmoidal_4pl(concs::AbstractVector{<:Real}, responses::AbstractVector{<:Real}) -> SigmoidCurve

Fit 4PL sigmoidal parameters via non-linear optimization grid search & refinement.
"""
function fit_sigmoidal_4pl(concs::AbstractVector{<:Real}, responses::AbstractVector{<:Real})
    length(concs) == length(responses) || throw(DimensionMismatch("concs and responses must match in length"))
    length(concs) >= 4 || throw(ArgumentError("At least 4 points required to fit 4PL curve"))

    c_vals = Float64.(concs)
    r_vals = Float64.(responses)

    # Initial estimations
    bottom_est = minimum(r_vals)
    top_est = maximum(r_vals)
    mid_resp = 0.5 * (bottom_est + top_est)

    # Find concentration closest to mid-response
    closest_idx = argmin(abs.(r_vals .- mid_resp))
    log_ec50_est = log10(max(c_vals[closest_idx], 1e-15))
    hill_est = 1.0

    # Grid refinement
    best_loss = Inf
    best_curve = SigmoidCurve(bottom_est, top_est, log_ec50_est; hill=hill_est)

    log_range = range(log_ec50_est - 1.5, stop=log_ec50_est + 1.5, length=40)
    hill_range = range(0.5, stop=2.5, length=20)

    for l_ec in log_range
        for h in hill_range
            sc = SigmoidCurve(bottom_est, top_est, l_ec; hill=h)
            loss = 0.0
            for i in eachindex(c_vals)
                pred = evaluate_4pl(sc, c_vals[i])
                loss += (pred - r_vals[i])^2
            end
            if loss < best_loss
                best_loss = loss
                best_curve = sc
            end
        end
    end

    return best_curve
end

"""
    SchildResult

Outcome of Schild regression analysis:
- `antagonist_concs`: Vector of antagonist concentrations [B] (M)
- `dose_ratios`: Vector of Dose Ratios (DR = EC50_antag / EC50_control)
- `log_b`: log10([B])
- `log_dr_minus_1`: log10(DR - 1)
- `slope`: Schild slope (expected ~ 1.0 for competitive antagonism)
- `pa2`: pA2 = -log10(KB)
- `kb`: Antagonist equilibrium dissociation constant KB (M)
- `r_squared`: Regression coefficient of determination
"""
struct SchildResult
    antagonist_concs::Vector{Float64}
    dose_ratios::Vector{Float64}
    log_b::Vector{Float64}
    log_dr_minus_1::Vector{Float64}
    slope::Float64
    pa2::Float64
    kb::Float64
    r_squared::Float64
end

"""
    perform_schild_analysis(ec50_control::Real, antagonist_concs::AbstractVector{<:Real},
                           ec50_shifted::AbstractVector{<:Real}) -> SchildResult

Perform classical Arunlakshana & Schild regression analysis:
    log10(DR - 1) = slope * log10([B]) + intercept
where DR = ec50_shifted / ec50_control.
Returns `SchildResult` with slope, pA2, and KB.
"""
function perform_schild_analysis(ec50_control::Real, antagonist_concs::AbstractVector{<:Real},
                               ec50_shifted::AbstractVector{<:Real})
    ec50_control > 0 || throw(ArgumentError("ec50_control must be positive"))
    length(antagonist_concs) == length(ec50_shifted) || throw(DimensionMismatch("concs and shifted EC50s must match"))
    length(antagonist_concs) >= 2 || throw(ArgumentError("At least 2 antagonist concentrations required for Schild regression"))

    drs = Float64[]
    log_b = Float64[]
    log_dr_1 = Float64[]

    for i in eachindex(antagonist_concs)
        b = Float64(antagonist_concs[i])
        ec50_s = Float64(ec50_shifted[i])
        b > 0 || throw(ArgumentError("Antagonist concentrations must be positive"))
        ec50_s > ec50_control || throw(ArgumentError("Shifted EC50 must be greater than control for right-shift"))

        dr = ec50_s / Float64(ec50_control)
        push!(drs, dr)
        push!(log_b, log10(b))
        push!(log_dr_1, log10(dr - 1.0))
    end

    # Linear regression: log(DR-1) = m * log(B) + c
    x = log_b
    y = log_dr_1
    n = length(x)

    x_mean = calc_mean(x)
    y_mean = calc_mean(y)

    ss_xx = sum((xi - x_mean)^2 for xi in x)
    ss_xy = sum((x[i] - x_mean) * (y[i] - y_mean) for i in 1:n)
    ss_yy = sum((yi - y_mean)^2 for yi in y)

    slope = ss_xy / ss_xx
    intercept = y_mean - slope * x_mean

    r_sq = (ss_xy^2) / (ss_xx * ss_yy)

    # When log(DR - 1) = 0: x_int = -intercept / slope
    # pA2 = -x_int = intercept / slope
    pa2 = -(-intercept / slope)
    kb = 10.0^(-pa2)

    return SchildResult(Float64.(antagonist_concs), drs, log_b, log_dr_1, slope, pa2, kb, r_sq)
end

end # module Simulator
