"""
    gpcr_biased_agonism.jl

Example script demonstrating GPCR biased agonism quantification using the Black-Leff
operational model of agonism and the Kenakin transducer ratio (ΔΔlog(tau/KA)).
"""

using BioAffinity
using Printf

function run_example()
    println("====================================================================")
    println(" GPCR Biased Agonism Analysis: G-Protein vs β-Arrestin Signaling")
    println("====================================================================")

    # Reference full agonist (e.g., Endogenous Opioid / Catecholamine):
    # Balanced transducer coupling across both signaling cascades
    ref_g = PathwayParams("G-protein Signaling (cAMP)", 50.0e-9, 20.0)
    ref_arr = PathwayParams("β-Arrestin-2 Recruitment", 80.0e-9, 15.0)
    ref_profile = AgonistProfile("Endogenous Reference Agonist", [ref_g, ref_arr])

    # Candidate 1: Balanced Partial Agonist
    cand1_g = PathwayParams("G-protein Signaling (cAMP)", 25.0e-9, 5.0)
    cand1_arr = PathwayParams("β-Arrestin-2 Recruitment", 40.0e-9, 3.75)
    cand1_profile = AgonistProfile("Candidate-A (Balanced Partial)", [cand1_g, cand1_arr])

    # Candidate 2: G-Protein Biased Agonist (Arrestin Sparing)
    # Designed to retain therapeutic analgesia while minimizing arrestin-mediated respiratory depression
    cand2_g = PathwayParams("G-protein Signaling (cAMP)", 10.0e-9, 30.0)
    cand2_arr = PathwayParams("β-Arrestin-2 Recruitment", 600.0e-9, 0.8)
    cand2_profile = AgonistProfile("Candidate-B (G-Protein Biased Lead)", [cand2_g, cand2_arr])

    # Candidate 3: Arrestin-Biased Ligand (G-protein Sparing)
    cand3_g = PathwayParams("G-protein Signaling (cAMP)", 800.0e-9, 0.5)
    cand3_arr = PathwayParams("β-Arrestin-2 Recruitment", 15.0e-9, 25.0)
    cand3_profile = AgonistProfile("Candidate-C (Arrestin Biased Tool)", [cand3_g, cand3_arr])

    candidates = [cand1_profile, cand2_profile, cand3_profile]

    println("Agonist Profiles & Black-Leff Operational Parameters:")
    @printf("%-35s | G-prot KA (nM) | G-prot tau | Arr KA (nM) | Arr tau\n", "Compound Name")
    println("-"^75)
    @printf("%-35s |   %8.2f     |  %8.2f  |  %8.2f   |  %8.2f\n",
            ref_profile.name, ref_g.ka * 1e9, ref_g.tau, ref_arr.ka * 1e9, ref_arr.tau)

    for cand in candidates
        p_g = cand.pathways[1]
        p_arr = cand.pathways[2]
        @printf("%-35s |   %8.2f     |  %8.2f  |  %8.2f   |  %8.2f\n",
                cand.name, p_g.ka * 1e9, p_g.tau, p_arr.ka * 1e9, p_arr.tau)
    end
    println("-"^75)
    println()

    println("Kenakin Pathway Bias Factor Evaluation (Relative to Reference):")
    @printf("%-35s | Δlog(G-prot) | Δlog(Arr) | ΔΔlog(Bias) | Bias Ratio | Favored Pathway\n", "Compound Name")
    println("-"^95)

    for cand in candidates
        res = calculate_bias_factor(cand, ref_profile)
        @printf("%-35s |    %+6.3f    |   %+6.3f  |   %+6.3f   |   %7.1fx  | %s\n",
                cand.name, res.delta_log_pathway1, res.delta_log_pathway2,
                res.delta_delta_log, res.bias_ratio, res.favored_pathway)
    end
    println("-"^95)
end

run_example()
