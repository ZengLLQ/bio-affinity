using Test
using BioAffinity
using LinearAlgebra

@testset "BioAffinity Test Suite" begin

    @testset "Module 1: Binding Kinetics & Residence Time" begin
        # Test 1: Single ligand association & equilibrium asymptote
        kon = 1e6      # 1e6 M^-1 s^-1
        koff = 0.01    # 0.01 s^-1 (Kd = 10 nM, tau = 100 s)
        rtot = 1e-9    # 1.0 nM Bmax
        ligand = 10e-9 # 10.0 nM ligand (= Kd)

        kin = BindingKinetics1Order(kon, koff, rtot, ligand)
        @test kin.kon == 1e6
        @test kin.koff == 0.01
        @test residence_time(kin.koff) ≈ 100.0
        @test half_life(kin.koff) ≈ 100.0 * log(2.0)

        # At t = 0, RL = 0
        @test association_1st_order(kin, 0.0) == 0.0

        # At t -> infinity, RL -> Bmax * [L] / (Kd + [L]) = 1.0 nM * 10 / (10 + 10) = 0.5 nM
        rl_inf = association_1st_order(kin, 1e6)
        @test isapprox(rl_inf, 0.5e-9; rtol=1e-5)

        # Dissociation from 1.0 nM
        @test dissociation_1st_order(1.0e-9, koff, 0.0) == 1.0e-9
        @test isapprox(dissociation_1st_order(1.0e-9, koff, 100.0), 1.0e-9 * exp(-1.0); rtol=1e-6)

        # Test 2: Motulsky & Mahan Competitive Kinetics
        # Tracer: k1=1e6 M^-1 s^-1, k2=0.01 s^-1
        # Competitor: k3=5e5 M^-1 s^-1, k4=0.001 s^-1
        # Bmax = 1.0 nM, Tracer = 5.0 nM, Competitor = 10.0 nM
        mm = MotulskyMahanParams(1e6, 0.01, 5e5, 0.001, 1e-9, 5e-9, 10e-9)

        # Check t = 0
        rl0, ri0 = motulsky_mahan_analytical(mm, 0.0)
        @test rl0 == 0.0
        @test ri0 == 0.0

        # Compare analytical solution with RK4 numerical ODE integration
        t_test = 600.0
        rl_ana, ri_ana = motulsky_mahan_analytical(mm, t_test)

        ode_times, ode_rl, ode_ri = motulsky_mahan_ode(mm, t_test; dt=0.05, depletion=false)
        @test isapprox(rl_ana, ode_rl[end]; rtol=1e-4)
        @test isapprox(ri_ana, ode_ri[end]; rtol=1e-4)

        # Check receptor conservation: RL + RI <= Rtot
        for t in [10.0, 50.0, 100.0, 300.0, 1000.0, 5000.0]
            r_l, r_i = motulsky_mahan_analytical(mm, t)
            @test (r_l + r_i) <= mm.rtot + 1e-15
            @test r_l >= 0.0
            @test r_i >= 0.0
        end

        # Test kinetic overshoot: tracer binding peaks before equilibrium
        sim_res = simulate_competitive_kinetics(mm, 2000.0; n_points=200)
        @test sim_res.peak_tracer_overshoot > sim_res.bound_tracer[end]
        @test sim_res.residence_time_competitor ≈ 1000.0
        @test sim_res.residence_time_tracer ≈ 100.0
    end

    @testset "Module 2: Equilibrium Saturation & Morrison Formulation" begin
        # Saturation binding
        bmax = 100.0 # fmol/mg
        kd = 5.0e-9  # 5 nM
        sat = SaturationBinding(bmax, kd; hill=1.0)

        @test saturation_fraction(sat, 0.0) == 0.0
        # At [L] = Kd, bound = 0.5 * Bmax
        @test saturation_fraction(sat, 5.0e-9) ≈ 50.0
        # At [L] = 9 * Kd, bound = 0.9 * Bmax
        @test isapprox(saturation_fraction(sat, 45.0e-9), 90.0; rtol=1e-6)

        # Scatchard transformation
        ligands = [1e-9, 2e-9, 5e-9, 10e-9, 20e-9, 50e-9]
        bound = [saturation_fraction(sat, l) for l in ligands]
        x_scat, y_scat = scatchard_transform(bound, ligands)
        @test length(x_scat) == length(ligands)
        @test all(y_scat .> 0)

        # Linear Scatchard fit
        kd_fit, bmax_fit, r2 = linear_scatchard_fit(bound, ligands)
        @test isapprox(kd_fit, kd; rtol=1e-5)
        @test isapprox(bmax_fit, bmax; rtol=1e-5)
        @test r2 > 0.9999

        # Cheng-Prusoff equations
        ic50 = 50e-9
        l_conc = 10e-9
        kd_val = 5e-9
        ki_comp = cheng_prusoff_competitive(ic50, l_conc, kd_val)
        # Ki = 50 / (1 + 10/5) = 50 / 3 = 16.6667 nM
        @test isapprox(ki_comp, 50e-9 / 3.0; rtol=1e-6)

        ki_uncomp = cheng_prusoff_uncompetitive(ic50, l_conc, kd_val)
        # Ki = 50 / (1 + 5/10) = 50 / 1.5 = 33.3333 nM
        @test isapprox(ki_uncomp, 50e-9 / 1.5; rtol=1e-6)

        @test cheng_prusoff_noncompetitive(ic50) == ic50

        # Morrison tight-binding equation
        # Etot = 10 nM, Ki = 1 nM (Etot >= Ki -> tight binding)
        morr = MorrisonParams(10e-9, 1e-9)
        # At Itot = 0, fractional velocity = 1.0
        @test morrison_fractional_velocity(morr, 0.0) == 1.0
        # At large Itot, fractional velocity -> 0
        @test morrison_fractional_velocity(morr, 1e-6) < 0.01

        # Apparent IC50 = Ki_app + 0.5 * Etot = 1 nM + 5 nM = 6 nM
        @test isapprox(morrison_ic50_apparent(morr), 6e-9; rtol=1e-6)

        # At Itot = IC50_apparent = 6 nM, v/v0 should be approximately 0.5
        v_at_ic50 = morrison_fractional_velocity(morr, 6e-9)
        @test isapprox(v_at_ic50, 0.5; atol=0.05)

        # Free ligand calculation with depletion
        l_tot = 10e-9
        r_tot = 5e-9
        kd_exact = 2e-9
        l_free = exact_free_ligand(l_tot, r_tot, kd_exact)
        rl_bound = r_tot * l_free / (kd_exact + l_free)
        # Invariant: l_free + rl_bound == l_tot
        @test isapprox(l_free + rl_bound, l_tot; rtol=1e-6)
    end

    @testset "Module 3: Allosteric Ternary Complex Model (ATCM) & MWC" begin
        # ATCM PAM with alpha = 10.0, beta = 2.0
        ka = 100e-9
        kb = 50e-9
        p_pam = AllostericTernaryParams(ka, kb; alpha=10.0, beta=2.0)

        # Free energy cycle path-independence
        cyc = thermodynamic_cycle_free_energy(p_pam)
        @test isapprox(cyc.delta_G_path1, cyc.delta_G_path2; atol=1e-10)
        # Coupling free energy is negative for PAM (favorable cooperativity)
        @test cyc.cooperativity_delta_G < 0.0

        # In absence of allosteric modulator (B = 0), apparent Kd is KA
        @test allosteric_apparent_kd(p_pam, 0.0) == ka

        # In presence of saturating PAM (B >> KB), apparent Kd approaches KA / alpha = 10 nM
        ka_sat = allosteric_apparent_kd(p_pam, 1e-3)
        @test isapprox(ka_sat, ka / 10.0; rtol=1e-3)

        # ATCM NAM with alpha = 0.1
        p_nam = AllostericTernaryParams(ka, kb; alpha=0.1, beta=0.5)
        cyc_nam = thermodynamic_cycle_free_energy(p_nam)
        @test isapprox(cyc_nam.delta_G_path1, cyc_nam.delta_G_path2; atol=1e-10)
        @test cyc_nam.cooperativity_delta_G > 0.0
        # In presence of NAM, apparent Kd increases (lower affinity)
        @test allosteric_apparent_kd(p_nam, 50e-9) > ka

        # MWC model (n = 4 tetramer, L0 = 1000, KR = 1e-6, KT = 1e-4)
        mwc = MWCParams(4, 1000.0, 1e-6, 1e-4)
        @test mwc_fractional_saturation(mwc, 0.0) == 0.0
        @test isapprox(mwc_fractional_saturation(mwc, 1.0), 1.0; atol=1e-4)

        hill_coeff = mwc_apparent_hill(mwc)
        @test hill_coeff > 1.5 # Significant positive cooperativity
    end

    @testset "Module 4: GPCR Biased Agonism & Operational Model" begin
        # Black-Leff operational model
        emax_sys = 100.0
        ka = 100e-9 # 100 nM
        tau_spare = 19.0 # large reserve

        op = OperationalModelParams(emax_sys, ka, tau_spare)
        @test operational_effect(op, 0.0) == 0.0

        # Observed Emax = Emax_sys * tau / (1 + tau) = 100 * 19 / 20 = 95.0
        @test operational_emax_observed(op) ≈ 95.0

        # EC50 = KA / (1 + tau) = 100 nM / 20 = 5 nM
        @test operational_ec50(op) ≈ 5.0e-9

        # Spare receptor occupancy at EC50: Y50 = EC50 / (KA + EC50) = 5 / 105 ≈ 4.76%
        @test isapprox(spare_receptor_fraction(op), 5e-9 / 105e-9; rtol=1e-5)

        # Kenakin Biased Agonism calculation
        p_ref_g = PathwayParams("cAMP", 50e-9, 10.0)
        p_ref_arr = PathwayParams("Arrestin", 50e-9, 10.0)
        ref_ag = AgonistProfile("Ref", [p_ref_g, p_ref_arr])

        # Test compound 1: perfectly balanced
        test_balanced = AgonistProfile("Balanced", [PathwayParams("cAMP", 10e-9, 2.0),
                                                   PathwayParams("Arrestin", 10e-9, 2.0)])
        bias_bal = calculate_bias_factor(test_balanced, ref_ag)
        @test isapprox(bias_bal.delta_delta_log, 0.0; atol=1e-6)
        @test isapprox(bias_bal.bias_ratio, 1.0; rtol=1e-5)
        @test bias_bal.favored_pathway == "Balanced (Unbiased)"

        # Test compound 2: 100-fold G-protein biased
        # cAMP: KA = 10 nM, tau = 20.0 -> tau/KA = 2e9
        # Arrestin: KA = 100 nM, tau = 0.2 -> tau/KA = 2e6
        # Bias ratio = (2e9 / 2e8) / (2e6 / 2e8) = 10 / 0.01 = 1000
        test_biased = AgonistProfile("BiasedDrug", [PathwayParams("cAMP", 10e-9, 20.0),
                                                    PathwayParams("Arrestin", 100e-9, 0.2)])
        bias_res = calculate_bias_factor(test_biased, ref_ag)
        @test bias_res.bias_ratio > 100.0
        @test bias_res.favored_pathway == "cAMP"
    end

    @testset "Module 5: Simulator & Schild Regression Analysis" begin
        # 4PL Curve
        sig = SigmoidCurve(0.0, 100.0, -8.0; hill=1.0)
        # At conc = 10^-8 M (10 nM), response should be 50.0
        @test isapprox(evaluate_4pl(sig, 10e-9), 50.0; atol=1e-6)

        concs, resps = generate_sigmoidal_curve(sig; log_min=-11.0, log_max=-5.0, n_points=50)
        @test length(concs) == 50
        @test resps[1] < 1.0
        @test resps[end] > 99.0

        # Fit 4PL curve
        fitted_sig = fit_sigmoidal_4pl(concs, resps)
        @test isapprox(fitted_sig.log_ec50, -8.0; atol=0.05)
        @test isapprox(fitted_sig.top, 100.0; atol=1.0)
        @test isapprox(fitted_sig.bottom, 0.0; atol=1.0)

        # Arunlakshana & Schild Regression
        ec50_control = 5.0e-9
        kb_true = 1.0e-9 # 1.0 nM antagonist dissociation constant
        antagonist_concs = [0.5e-9, 1.0e-9, 2.0e-9, 5.0e-9, 10.0e-9]
        # Shifted EC50 = EC50_ctrl * (1 + [B] / KB)
        ec50_shifted = [ec50_control * (1.0 + b / kb_true) for b in antagonist_concs]

        schild = perform_schild_analysis(ec50_control, antagonist_concs, ec50_shifted)
        # Unit slope criterion
        @test isapprox(schild.slope, 1.0000; atol=1e-4)
        # KB matches true value 1.0 nM (pA2 = 9.0)
        @test isapprox(schild.pa2, 9.0; atol=1e-4)
        @test isapprox(schild.kb, kb_true; rtol=1e-4)
        @test schild.r_squared > 0.9999
    end

    @testset "Module 6: CLI & Full Execution" begin
        # Verify CLI execution runs cleanly with all subcommands
        @test_nowarn BioAffinity.main(["kinetics-sim"])
        @test_nowarn BioAffinity.main(["morrison-fit"])
        @test_nowarn BioAffinity.main(["allostery-eval"])
        @test_nowarn BioAffinity.main(["bias-calc"])
        @test_nowarn BioAffinity.main(["schild-analysis"])
        @test_nowarn BioAffinity.main(["all"])
    end
end
