"""
    biased_agonism.jl

Black-Leff operational model of pharmacological agonism, transducer ratio (tau / KA),
spare receptor reserve quantification, and Kenakin pathway bias factor (ΔΔlog(tau/KA)).
"""

module BiasedAgonism

using LinearAlgebra

export OperationalModelParams, PathwayParams, AgonistProfile, BiasAnalysisResult
export operational_effect, operational_ec50, operational_emax_observed, spare_receptor_fraction
export transducer_ratio, pathway_delta_log_tau_ka, calculate_bias_factor

"""
    OperationalModelParams

Parameters for Black-Leff Operational Model of Agonism:
- `emax_system`: Theoretical ceiling maximal effect of the tissue/system
- `ka`: Agonist equilibrium dissociation constant (M)
- `tau`: Transducer efficacy ratio tau = [R_tot] / KE (dimensionless)
- `slope`: Transducer transducer slope (default 1.0)
"""
struct OperationalModelParams
    emax_system::Float64
    ka::Float64
    tau::Float64
    slope::Float64

    function OperationalModelParams(emax_system::Real, ka::Real, tau::Real; slope::Real=1.0)
        emax_system > 0 || throw(ArgumentError("emax_system must be positive"))
        ka > 0 || throw(ArgumentError("ka must be positive"))
        tau > 0 || throw(ArgumentError("tau must be positive"))
        slope > 0 || throw(ArgumentError("slope must be positive"))
        new(Float64(emax_system), Float64(ka), Float64(tau), Float64(slope))
    end
end

"""
    operational_effect(p::OperationalModelParams, agonist::Real) -> Float64

Compute functional pharmacological effect using the Black-Leff operational model:
    Effect([A]) = (Emax_sys * tau^n * [A]^n) / (KA^n + (1 + tau^n) * [A]^n)
"""
function operational_effect(p::OperationalModelParams, agonist::Real)
    agonist >= 0 || throw(ArgumentError("agonist concentration must be non-negative"))
    if agonist == 0.0
        return 0.0
    end

    a_n = Float64(agonist)^p.slope
    ka_n = p.ka^p.slope
    tau_n = p.tau^p.slope

    num = p.emax_system * tau_n * a_n
    denom = ka_n + (1.0 + tau_n) * a_n

    return num / denom
end

"""
    operational_ec50(p::OperationalModelParams) -> Float64

Calculate apparent potency EC50 under the operational model:
    EC50 = KA / (1 + tau^n)^(1/n)
"""
function operational_ec50(p::OperationalModelParams)
    tau_n = p.tau^p.slope
    return p.ka / ((1.0 + tau_n)^(1.0 / p.slope))
end

"""
    operational_emax_observed(p::OperationalModelParams) -> Float64

Calculate observable asymptotic maximum response E_obs,max:
    E_obs,max = Emax_sys * tau^n / (1 + tau^n)
"""
function operational_emax_observed(p::OperationalModelParams)
    tau_n = p.tau^p.slope
    return p.emax_system * tau_n / (1.0 + tau_n)
end

"""
    spare_receptor_fraction(p::OperationalModelParams) -> Float64

Fraction of total receptors occupied at the EC50 concentration (Y50):
    Y50 = EC50 / (KA + EC50) = 1 / (1 + (1 + tau^n)^(1/n))
When tau >> 1, Y50 << 0.5, revealing high spare receptor reserve.
"""
function spare_receptor_fraction(p::OperationalModelParams)
    ec50 = operational_ec50(p)
    return ec50 / (p.ka + ec50)
end

"""
    transducer_ratio(tau::Real, ka::Real) -> Float64

Compute log10(tau / KA) for a given pathway.
"""
function transducer_ratio(tau::Real, ka::Real)
    tau > 0 || throw(ArgumentError("tau must be positive"))
    ka > 0 || throw(ArgumentError("ka must be positive"))
    return log10(Float64(tau) / Float64(ka))
end

"""
    PathwayParams

Operational parameters for an agonist on a specific signaling pathway:
- `name`: Pathway name (e.g., "G-protein / cAMP", "beta-Arrestin recruitment")
- `ka`: Dissociation constant KA (M)
- `tau`: Operational efficacy tau
"""
struct PathwayParams
    name::String
    ka::Float64
    tau::Float64

    function PathwayParams(name::AbstractString, ka::Real, tau::Real)
        ka > 0 || throw(ArgumentError("ka must be positive"))
        tau > 0 || throw(ArgumentError("tau must be positive"))
        new(String(name), Float64(ka), Float64(tau))
    end
end

"""
    AgonistProfile

Multi-pathway profile for a pharmacological compound.
"""
struct AgonistProfile
    name::String
    pathways::Vector{PathwayParams}
end

"""
    pathway_delta_log_tau_ka(test_p::PathwayParams, ref_p::PathwayParams) -> Float64

Compute relative transducer performance compared to reference:
    Δlog10(tau / KA) = log10(tau_test / KA_test) - log10(tau_ref / KA_ref)
"""
function pathway_delta_log_tau_ka(test_p::PathwayParams, ref_p::PathwayParams)
    test_ratio = transducer_ratio(test_p.tau, test_p.ka)
    ref_ratio = transducer_ratio(ref_p.tau, ref_p.ka)
    return test_ratio - ref_ratio
end

"""
    BiasAnalysisResult

Result of Kenakin bias analysis across two signaling pathways.
"""
struct BiasAnalysisResult
    agonist_name::String
    reference_name::String
    pathway1_name::String
    pathway2_name::String
    delta_log_pathway1::Float64
    delta_log_pathway2::Float64
    delta_delta_log::Float64
    bias_ratio::Float64
    favored_pathway::String
end

"""
    calculate_bias_factor(test::AgonistProfile, ref::AgonistProfile,
                          pathway1_idx::Int=1, pathway2_idx::Int=2) -> BiasAnalysisResult

Calculate Kenakin pathway bias factor ΔΔlog10(tau / KA) and bias ratio:
    ΔΔlog(tau/KA) = Δlog(tau/KA)_pathway1 - Δlog(tau/KA)_pathway2
    Bias Ratio = 10^(ΔΔlog(tau/KA))
"""
function calculate_bias_factor(test::AgonistProfile, ref::AgonistProfile,
                               pathway1_idx::Int=1, pathway2_idx::Int=2)
    p1_test = test.pathways[pathway1_idx]
    p1_ref = ref.pathways[pathway1_idx]
    p2_test = test.pathways[pathway2_idx]
    p2_ref = ref.pathways[pathway2_idx]

    d_log1 = pathway_delta_log_tau_ka(p1_test, p1_ref)
    d_log2 = pathway_delta_log_tau_ka(p2_test, p2_ref)

    dd_log = d_log1 - d_log2
    ratio = 10.0^dd_log

    favored = if ratio > 1.05
        p1_test.name
    elseif ratio < 0.95
        p2_test.name
    else
        "Balanced (Unbiased)"
    end

    return BiasAnalysisResult(test.name, ref.name, p1_test.name, p2_test.name,
                              d_log1, d_log2, dd_log, ratio, favored)
end

end # module BiasedAgonism
