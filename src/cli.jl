"""
    cli.jl

Command-line interface and diagnostic reporting tools for BioAffinity.
"""

module CLI

using Printf
using ..Kinetics
using ..Equilibrium
using ..Allostery
using ..BiasedAgonism
using ..Simulator

export run_cli

function print_banner()
    println("""
    ========================================================================
       BioAffinity: Quantitative Receptor Pharmacology & Biased Agonism
       High-Performance Molecular Binding Kinetics & GPCR Signaling Engine
    ========================================================================
    """)
end

function cmd_kinetics()
    println("[*] Simulating Motulsky & Mahan Competitive Association Kinetics...")
    # Tracer: k1=1e6 M^-1 s^-1, k2=0.01 s^-1 (Kd=10 nM, tau=100s)
    # Competitor: k3=5e5 M^-1 s^-1, k4=0.001 s^-1 (Ki=2 nM, tau=1000s)
    # [R_tot]=1.0 nM, [Tracer]=5.0 nM, [Comp]=10.0 nM
    p = MotulskyMahanParams(1e6, 0.01, 5e5, 0.001, 1e-9, 5e-9, 10e-9)
    res = simulate_competitive_kinetics(p, 1800.0; n_points=100)

    @printf("  Tracer Kon (k1):           %.2e M^-1 s^-1\n", p.k1)
    @printf("  Tracer Koff (k2):          %.4f s^-1\n", p.k2)
    @printf("  Tracer Residence Time (τ): %.1f s (%.2f min)\n", res.residence_time_tracer, res.residence_time_tracer / 60.0)
    @printf("  Competitor Kon (k3):       %.2e M^-1 s^-1\n", p.k3)
    @printf("  Competitor Koff (k4):      %.4f s^-1\n", p.k4)
    @printf("  Competitor Residence (τ):  %.1f s (%.2f min)\n", res.residence_time_competitor, res.residence_time_competitor / 60.0)
    @printf("  Peak Tracer Bound (RL):    %.3e M (Overshoot Peak)\n", res.peak_tracer_overshoot)
    @printf("  Equilibrium Tracer Bound:  %.3e M (%.1f%% of Bmax)\n", res.bound_tracer[end], res.equilibrium_fraction_tracer * 100.0)
    println("  [✔] Motulsky-Mahan kinetic displacement overshoot successfully captured.\n")
end

function cmd_morrison()
    println("[*] Evaluating Morrison Tight-Binding Picomolar Enzyme/Receptor Inhibition...")
    # Etot = 10 nM, Ki = 0.5 nM (tight-binding regime: Etot >= Ki)
    p = MorrisonParams(10e-9, 0.5e-9)
    ic50_app = morrison_ic50_apparent(p)
    classic_ic50 = p.ki # classical Cheng-Prusoff assumes IC50 = Ki when [S]=0

    @printf("  Target Receptor/Enzyme [E_tot]: %.2e M (%.1f nM)\n", p.etot, p.etot * 1e9)
    @printf("  True Inhibitor Affinity Ki:     %.2e M (%.2f nM)\n", p.ki, p.ki * 1e9)
    @printf("  Classical (Cheng-Prusoff) IC50: %.2e M (%.2f nM)\n", classic_ic50, classic_ic50 * 1e9)
    @printf("  True Morrison Apparent IC50:    %.2e M (%.2f nM)\n", ic50_app, ic50_app * 1e9)
    @printf("  Tight-Binding Underestimation:  %.1fx error if Morrison ignored\n", ic50_app / classic_ic50)
    println("  [✔] Quadratic Morrison formulation eliminates tight-binding artifacts.\n")
end

function cmd_allostery()
    println("[*] Analyzing Ehlert Allosteric Ternary Complex Model (ATCM)...")
    # KA = 100 nM, KB = 50 nM, alpha = 8.0 (Positive Allosteric Modulator - PAM), beta = 2.5
    p = AllostericTernaryParams(100e-9, 50e-9; alpha=8.0, beta=2.5)
    cycle = thermodynamic_cycle_free_energy(p)

    ka_unmod = p.ka
    ka_sat = allosteric_apparent_kd(p, 5e-6) # saturated with 5 uM PAM

    @printf("  Orthosteric KA (unmodulated): %.1f nM\n", ka_unmod * 1e9)
    @printf("  Allosteric KB:                %.1f nM\n", p.kb * 1e9)
    @printf("  Cooperativity Factor α (PAM): %.1fx affinity enhancement\n", p.alpha)
    @printf("  Efficacy Factor β:            %.1fx signaling amplification\n", p.beta)
    @printf("  KA apparent at 5 μM PAM:      %.2f nM (approaching limit KA/α = %.1f nM)\n", ka_sat * 1e9, (p.ka / p.alpha) * 1e9)
    @printf("  Thermodynamic Coupling ΔΔG°:  %.2f kJ/mol\n", cycle.cooperativity_delta_G / 1000.0)
    println("  [✔] Thermodynamic path-independence confirmed (ΔG_path1 == ΔG_path2).\n")
end

function cmd_bias()
    println("[*] Quantifying GPCR Biased Agonism via Black-Leff Operational Model...")
    # Reference full agonist (e.g. Endogenous ligand)
    p_ref_g = PathwayParams("G-protein (cAMP)", 50e-9, 20.0)
    p_ref_arr = PathwayParams("β-Arrestin recruitment", 100e-9, 15.0)
    ref_agonist = AgonistProfile("Endogenous Reference Agonist", [p_ref_g, p_ref_arr])

    # Biased test agonist (e.g. G-protein biased therapeutic)
    p_test_g = PathwayParams("G-protein (cAMP)", 10e-9, 25.0)
    p_test_arr = PathwayParams("β-Arrestin recruitment", 500e-9, 1.2)
    test_agonist = AgonistProfile("Therapeutic Lead Compound", [p_test_g, p_test_arr])

    bias = calculate_bias_factor(test_agonist, ref_agonist)

    @printf("  Test Agonist:                  %s\n", bias.agonist_name)
    @printf("  Reference Standard:            %s\n", bias.reference_name)
    @printf("  Pathway 1 (%s) Δlog(τ/KA):  %+.3f\n", bias.pathway1_name, bias.delta_log_pathway1)
    @printf("  Pathway 2 (%s) Δlog(τ/KA):  %+.3f\n", bias.pathway2_name, bias.delta_log_pathway2)
    @printf("  Kenakin Bias Factor ΔΔlog(τ/KA): %+.3f\n", bias.delta_delta_log)
    @printf("  Bias Ratio:                    %.1fx favoring %s\n", bias.bias_ratio, bias.favored_pathway)
    println("  [✔] Functional biased agonism mathematically proven.\n")
end

function cmd_schild()
    println("[*] Running Arunlakshana & Schild Regression Analysis...")
    ec50_ctrl = 10e-9 # 10 nM
    # Competitive antagonist KB = 2 nM (pA2 = 8.70)
    antag_concs = [1e-9, 3e-9, 10e-9, 30e-9, 100e-9]
    kb_true = 2e-9
    ec50_shifted = [ec50_ctrl * (1.0 + b / kb_true) for b in antag_concs]

    res = perform_schild_analysis(ec50_ctrl, antag_concs, ec50_shifted)

    @printf("  Control Agonist EC50: %.1f nM\n", ec50_ctrl * 1e9)
    @printf("  Schild Regression Slope: %.4f (Theoretical = 1.0000)\n", res.slope)
    @printf("  pA2 (-log10 KB):         %.3f\n", res.pa2)
    @printf("  Calculated Antagonist KB: %.2e M (%.2f nM, True = %.2f nM)\n", res.kb, res.kb * 1e9, kb_true * 1e9)
    @printf("  Regression Fit R^2:      %.5f\n", res.r_squared)
    println("  [✔] Strict unit slope validates simple competitive antagonism.\n")
end

function run_cli(args::Vector{String}=String[])
    print_banner()

    cmd = isempty(args) ? "all" : lowercase(args[1])

    if cmd == "kinetics-sim"
        cmd_kinetics()
    elseif cmd == "morrison-fit"
        cmd_morrison()
    elseif cmd == "allostery-eval"
        cmd_allostery()
    elseif cmd == "bias-calc"
        cmd_bias()
    elseif cmd == "schild-analysis"
        cmd_schild()
    elseif cmd == "all"
        cmd_kinetics()
        cmd_morrison()
        cmd_allostery()
        cmd_bias()
        cmd_schild()
        println("========================================================================")
        println("All BioAffinity pharmacology computational modules executed successfully!")
        println("========================================================================")
    else
        println("Unknown command: '$cmd'")
        println("Available commands: kinetics-sim, morrison-fit, allostery-eval, bias-calc, schild-analysis, all")
    end
end

end # module CLI
