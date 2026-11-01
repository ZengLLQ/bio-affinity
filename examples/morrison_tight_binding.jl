"""
    morrison_tight_binding.jl

Example script demonstrating the Morrison tight-binding equation for high-affinity
picomolar enzyme inhibitors or receptor blockers.
Contrasts the exact quadratic Morrison formulation with classical Cheng-Prusoff IC50.
"""

using BioAffinity
using Printf

function run_example()
    println("====================================================================")
    println(" Morrison Tight-Binding Picomolar Target Inhibition vs Cheng-Prusoff")
    println("====================================================================")

    # Enzyme / Receptor assay setup:
    # Target concentration [E_tot] = 10.0 nM
    # Substrate concentration [S] = 50.0 uM, Km = 25.0 uM ([S]/Km = 2.0)
    # Three inhibitors spanning weak to ultra-potent tight-binding:
    #   Inhibitor 1: Ki = 500.0 nM (weak binding: Ki >> [E_tot])
    #   Inhibitor 2: Ki = 10.0 nM (intermediate: Ki ≈ [E_tot])
    #   Inhibitor 3: Ki = 0.1 nM = 100 pM (tight binding: Ki << [E_tot])

    etot = 10.0e-9
    sub = 50.0e-6
    km = 25.0e-6

    inhibitors = [
        ("Weak Inhibitor (Classical)", 500.0e-9),
        ("Moderate Inhibitor (Borderline)", 10.0e-9),
        ("Ultra-Potent Inhibitor (Tight-Binding)", 0.1e-9)
    ]

    println("Assay Conditions:")
    @printf("  Enzyme/Receptor Concentration [E_tot]: %.2f nM\n", etot * 1e9)
    @printf("  Substrate [S] / Km Ratio:              %.1f\n", sub / km)
    println()

    @printf("%-40s | True Ki (nM) | Cheng-Prusoff IC50 | Morrison IC50 | Error Factor\n", "Inhibitor Profile")
    println("-"^95)

    for (name, ki) in inhibitors
        p = MorrisonParams(etot, ki; substrate=sub, km=km)
        # Classical Cheng-Prusoff: IC50 = Ki * (1 + [S]/Km)
        ic50_classic = ki * (1.0 + sub / km)
        # Morrison tight-binding apparent IC50: IC50 = Ki*(1 + [S]/Km) + 0.5 * [E_tot]
        ic50_morrison = morrison_ic50_apparent(p)
        error_ratio = ic50_morrison / ic50_classic

        @printf("%-40s |   %9.3f  |      %8.3f nM   |   %8.3f nM  |   %.2fx\n",
                name, ki * 1e9, ic50_classic * 1e9, ic50_morrison * 1e9, error_ratio)
    end
    println("-"^95)
    println("\nKey Pharmacological Insight:")
    println("For the ultra-potent inhibitor (Ki = 0.1 nM), the classical Cheng-Prusoff equation")
    println("predicts an IC50 of 0.30 nM. However, because [E_tot] = 10 nM, half the enzyme cannot")
    println("be inhibited until at least 5 nM of inhibitor is added! The true IC50 is 5.30 nM.")
    println("Neglecting Morrison tight-binding leads to an apparent 17.67-fold underestimation of potency.")
end

run_example()
