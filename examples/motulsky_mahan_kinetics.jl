"""
    motulsky_mahan_kinetics.jl

Example script demonstrating the Motulsky & Mahan model of competitive association kinetics.
Enables extraction of unlabeled competitor on-rate (k_on), off-rate (k_off), and drug-target
residence time (tau) by monitoring tracer kinetic displacement.
"""

using BioAffinity
using Printf

function run_example()
    println("====================================================================")
    println(" Motulsky & Mahan Competitive Binding Kinetics: Residence Time Extraction")
    println("====================================================================")

    # System parameters:
    # Fluorescent / Radioligand Tracer:
    #   k1 = 1.0e6 M^-1 s^-1, k2 = 0.02 s^-1 (Kd = 20 nM, tau = 50 s)
    # Unlabeled Drug Candidate (Slow Dissociation / Long Residence Time):
    #   k3 = 4.0e5 M^-1 s^-1, k4 = 0.0005 s^-1 (Ki = 1.25 nM, tau = 2000 s = 33.3 min)
    # Receptor concentration Bmax = 2.0 nM
    # Tracer concentration [L] = 10.0 nM
    # Competitor concentrations: 0 nM (control), 2 nM, 10 nM, 50 nM

    bmax = 2.0e-9
    tracer_conc = 10.0e-9
    k1 = 1.0e6
    k2 = 0.02
    k3 = 4.0e5
    k4 = 0.0005

    println("Experimental Assay Setup:")
    @printf("  Receptor Total (Bmax):           %.2f nM\n", bmax * 1e9)
    @printf("  Tracer Concentration [L]:        %.2f nM\n", tracer_conc * 1e9)
    @printf("  Tracer Kd:                       %.2f nM  (tau = %.1f s)\n", (k2 / k1) * 1e9, 1.0 / k2)
    @printf("  Candidate Drug Ki:               %.2f nM  (tau = %.1f s = %.1f min)\n",
            (k4 / k3) * 1e9, 1.0 / k4, (1.0 / k4) / 60.0)
    println()

    competitor_doses = [0.0, 2.0e-9, 10.0e-9, 50.0e-9]
    t_span = 3600.0 # 1 hour assay
    time_points = [60.0, 180.0, 300.0, 600.0, 1200.0, 1800.0, 3600.0]

    println("Tracer Binding Time Course [RL](t) (nM) Across Competitor Doses:")
    print("Time (min) | ")
    for c in competitor_doses
        @printf("Comp=%-5.1fnM | ", c * 1e9)
    end
    println()
    println("-"^65)

    for t in time_points
        @printf("%9.1f | ", t / 60.0)
        for c in competitor_doses
            if c == 0.0
                # Single-ligand association
                kin1 = BindingKinetics1Order(k1, k2, bmax, tracer_conc)
                val_rl = association_1st_order(kin1, t)
            else
                p_mm = MotulskyMahanParams(k1, k2, k3, k4, bmax, tracer_conc, c)
                val_rl, _ = motulsky_mahan_analytical(p_mm, t)
            end
            @printf("  %7.4f     | ", val_rl * 1e9)
        end
        println()
    end

    println("-"^65)
    println("\nKey Pharmacological Insight:")
    println("Notice the kinetic 'overshoot': in the presence of 10 nM and 50 nM competitor,")
    println("the tracer associates rapidly due to high Kon, reaching a peak bound state,")
    println("before the slow-dissociating competitor gradually replaces it at equilibrium.")
    println("This characteristic displacement curve confirms high drug-target residence time.")
end

run_example()
