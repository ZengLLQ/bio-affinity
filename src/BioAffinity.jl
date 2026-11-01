"""
    BioAffinity.jl

# BioAffinity: Quantitative Receptor Pharmacology, Drug-Target Binding Kinetics & GPCR Biased Agonism Engine

High-performance computational molecular pharmacology framework implemented in pure Julia.
Provides exact mathematical solutions and high-order numerical algorithms for:
- Association and dissociation kinetics, drug-target residence time (τ = 1 / koff).
- The Motulsky & Mahan model of competitive association kinetics with kinetic overshoot modeling.
- Equilibrium saturation binding, Scatchard transformation, and Cheng-Prusoff transformations.
- Morrison tight-binding quadratic formulations for high-affinity picomolar targets.
- Ehlert & Hall Allosteric Ternary Complex Model (ATCM) with binding (α) and efficacy (β) cooperativity.
- Monod-Wyman-Changeux (MWC) concerted transition allosteric oligomer models.
- Black-Leff operational model of agonism and receptor reserve / spare receptor quantification.
- Kenakin functional GPCR pathway biased agonism factor (ΔΔlog(τ / KA)).
- Arunlakshana & Schild regression analysis for competitive receptor antagonists.
"""
module BioAffinity

include("kinetics.jl")
include("equilibrium.jl")
include("allostery.jl")
include("biased_agonism.jl")
include("simulator.jl")
include("cli.jl")

using .Kinetics
using .Equilibrium
using .Allostery
using .BiasedAgonism
using .Simulator
using .CLI

# Re-export key interfaces
export Kinetics, Equilibrium, Allostery, BiasedAgonism, Simulator, CLI

# Kinetics exports
export BindingKinetics1Order, MotulskyMahanParams, KineticsResult
export association_1st_order, dissociation_1st_order, residence_time, half_life
export motulsky_mahan_analytical, motulsky_mahan_ode, simulate_competitive_kinetics

# Equilibrium exports
export SaturationBinding, MorrisonParams
export saturation_fraction, scatchard_transform, linear_scatchard_fit
export cheng_prusoff_competitive, cheng_prusoff_uncompetitive, cheng_prusoff_noncompetitive, cheng_prusoff_mixed
export morrison_fractional_velocity, morrison_ic50_apparent, exact_free_ligand

# Allostery exports
export AllostericTernaryParams, MWCParams
export allosteric_orthosteric_occupancy, allosteric_apparent_kd, allosteric_functional_response
export thermodynamic_cycle_free_energy, mwc_fractional_saturation, mwc_apparent_hill

# BiasedAgonism exports
export OperationalModelParams, PathwayParams, AgonistProfile, BiasAnalysisResult
export operational_effect, operational_ec50, operational_emax_observed, spare_receptor_fraction
export transducer_ratio, pathway_delta_log_tau_ka, calculate_bias_factor

# Simulator exports
export SigmoidCurve, SchildResult
export generate_sigmoidal_curve, evaluate_4pl, fit_sigmoidal_4pl, perform_schild_analysis

# CLI exports
export run_cli

"""
    main(args=ARGS)

Entry point for command-line execution.
"""
function main(args=ARGS)
    CLI.run_cli(args)
end

end # module BioAffinity
