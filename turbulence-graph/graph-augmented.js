// Augmented nodes and edges from research agents
// Appended to GRAPH_DATA after graph-data.js loads

(function() {
  const extra_nodes = [

    // ── From Assumptions Agent ───────────────────────────────────────────
    {
      id: "k41_dimensional_analysis",
      label: "K41 Dimensional Analysis Engine",
      type: "framework",
      angle: "assumptions",
      importance: 9,
      detail: "The logical core of K41: given universality + locality + isotropy + stationarity + ergodicity, the ONLY combination of ε [m²/s³] and k [1/m] with dimensions of energy spectral density [m³/s²] is ε^(2/3) k^(−5/3). The −5/3 law is a theorem ABOUT THE ASSUMPTIONS, not an independent discovery. Any measured deviation from −5/3 (empirical range −1.67 to −1.72) is either finite-Re artifact, measurement error, or evidence that at least one assumption fails.",
    },
    {
      id: "universality_violations",
      label: "Universality Violations (C_K varies 30%)",
      type: "empirical",
      angle: "assumptions",
      importance: 8,
      detail: "The Kolmogorov 'constant' C_K in E(k) = C_K ε^(2/3) k^(−5/3) varies from ~1.4 to 2.2 across different flow configurations (pipe, jet, wake, atmosphere, DNS). Sreenivasan (1995) review: C_K = 1.62 ± 0.17. This ~30% variation is not explained by finite-Re effects — it persists at Re_λ > 1000. Polymer additives alter small-scale statistics at ppm concentration. Near-wall turbulence maintains anisotropy down to viscous scales. C_K variation is the empirical fingerprint of universality failure.",
      refs: ["Sreenivasan 1995"]
    },
    {
      id: "closure_underdetermination",
      label: "Closure is Provably Unsolvable (No Scale Sep.)",
      type: "proven",
      angle: "assumptions",
      importance: 9,
      detail: "Unlike the BBGKY hierarchy in kinetic theory — closed by molecular chaos (Stosszahlansatz) because molecular collisions are local and instantaneous — turbulent interactions are non-local, multi-scale, and persistent. There is no turbulence analog of 'molecular chaos' to justify truncation. Kraichnan's Direct Interaction Approximation (DIA, 1959) is the most sophisticated closure attempt; it fails to reproduce −5/3 in Eulerian form because it cannot handle non-local sweeping. This failure is provably tied to the sweeping problem — a fundamentally non-local effect that invalidates any local closure.",
      refs: ["Kraichnan 1959", "McComb 1990"]
    },
    {
      id: "backscatter_magnitude",
      label: "Backscatter: 30–40% of Spacetime",
      type: "empirical",
      angle: "assumptions",
      importance: 8,
      detail: "DNS studies (Piomelli et al. 1991, Carati et al. 2001) show instantaneous backscatter events (negative subgrid energy flux) occur 30–40% of the time in space-time, and the instantaneous backscatter flux can be comparable in magnitude to the forward flux. The net forward cascade is the difference of two large numbers. Germano's dynamic Smagorinsky model implicitly handles this by allowing a sign change in the constant — a practical acknowledgment that the forward-only assumption fails locally.",
      refs: ["Piomelli et al. 1991", "Carati et al. 2001"]
    },
    {
      id: "re_scaling_inertial_range",
      label: "Clean Inertial Range Needs Re > 3500",
      type: "proven",
      angle: "assumptions",
      importance: 8,
      detail: "K41 itself gives L/η ~ Re^(3/4). The inertial range is contaminated: forcing effects extend down to ~0.1L from above; viscous effects extend up to ~60η from below. A clean decade of inertial range requires L/η > 600, i.e., Re > 600^(4/3) ≈ 3500. The world-record DNS (Kaneda et al. 2003, 4096³ grid, Re_λ ≈ 1200) shows approximately 1.5 decades of −5/3 scaling — barely enough for statistical tests of the theory it is testing. The 'inertial range' is partly idealization.",
      refs: ["Kaneda et al. 2003"]
    },
    {
      id: "kolmogorov_constant_problem",
      label: "Kolmogorov Constant C_K ≈ 1.5 (Empirical)",
      type: "assumption",
      angle: "assumptions",
      importance: 7,
      detail: "The Kolmogorov constant C_K in E(k) = C_K ε^(2/3) k^(−5/3) is not derived from theory — it is measured empirically. K41 universality predicts it should be universal. Its measured value C_K ≈ 1.5–1.7 (Sreenivasan 1995) but with ~30% variation across flows. This is analogous to the fine structure constant before QED — a measured parameter that a deeper theory should derive. No mechanism exists to compute C_K from N-S.",
    },
    {
      id: "sweeping_decorrelation",
      label: "Random Sweeping Effect (Kraichnan 1964)",
      type: "empirical",
      angle: "cross_domain",
      importance: 8,
      year: 1964,
      detail: "Kraichnan (1964): large-scale velocity sweeps small eddies without distorting them, introducing spurious decorrelation in Eulerian statistics. The Eulerian frequency spectrum has ω^(−2) exponent (not ω^(−5/3)) due to sweeping. Lagrangian spectrum recovers K41 scaling. Physical implication: small scales 'remember' the large-scale velocity that sweeps them — a non-Markovian correlation across scale space. This is the physical mechanism behind the Markovian cascade assumption's failure and why DIA fails in Eulerian form.",
      refs: ["Kraichnan 1964"]
    },

    // ── From Geometry/Topology Agent ────────────────────────────────────
    {
      id: "helicity_cascade_dual",
      label: "Dual Cascade: Energy + Helicity",
      type: "empirical",
      angle: "geometry_topology",
      importance: 7,
      detail: "In 3D turbulence, both energy and helicity cascade forward simultaneously. Biferale et al. showed that same-helicity-sign triads transfer energy backward (a hidden inverse cascade channel) while opposite-sign triads transfer forward. Standard K41 ignores helicity dynamics entirely. A complete theory must account for this dual-cascade structure: the interplay between energy and helicity cascades is not captured by ε alone, requiring a second invariant.",
      refs: ["Biferale et al. 2012"]
    },
    {
      id: "sle_percolation_universality",
      label: "SLE(6) = Critical Percolation Universality",
      type: "proven",
      angle: "geometry_topology",
      importance: 9,
      detail: "SLE(κ) with κ = 6 is the universality class of critical site percolation on triangular lattice (Smirnov 2001, Fields Medal). κ = 6 has the 'locality' property — the curve doesn't 'feel' the domain boundary. That 2D turbulence vorticity cluster boundaries match SLE(6)/percolation means the stream function in the inverse cascade behaves like a critical random field of percolation type. This is unexplained from N-S first principles — why should turbulence select the percolation universality class?",
      refs: ["Smirnov 2001", "Bernard et al. 2006"]
    },
    {
      id: "knotted_vortex_experiment",
      label: "Knotted Vortex Experiments (Kleckner-Irvine 2013)",
      type: "empirical",
      angle: "geometry_topology",
      importance: 8,
      year: 2013,
      detail: "Kleckner and Irvine created trefoil knot vortices and linked vortex rings in water. Key finding: knotted vortices inevitably untie through reconnections within a few crossing times, exciting Kelvin waves on resulting unknotted rings. Total helicity is approximately conserved through this process — writhe (topology) converts to twist (geometry). Knottedness is not a stable constraint: knots are transient in viscous flows. Higher knot invariants (Jones polynomial) do not provide robust conserved quantities.",
      refs: ["Kleckner & Irvine 2013"]
    },
    {
      id: "beltrami_fields",
      label: "Beltrami Fields: Any Knot Type Realizable",
      type: "proven",
      angle: "geometry_topology",
      importance: 7,
      year: 2015,
      detail: "Enciso and Peralta-Salas (2015, Ann. Math.): for ANY prescribed finite collection of knots and links, there exists a smooth Beltrami field on ℝ³ (∇×u = λu) whose vortex lines realize that topology. Beltrami fields are the only non-trivial steady Euler solutions in bounded domains. This means topologically complex vortex configurations are geometrically possible in Euler flow — but stability and relevance to turbulent N-S is open.",
      refs: ["Enciso & Peralta-Salas 2015"]
    },
    {
      id: "constantin_fefferman_criterion",
      label: "Constantin-Fefferman Geometric Regularity Criterion",
      type: "proven",
      angle: "geometry_topology",
      importance: 9,
      year: 1993,
      detail: "Constantin and Fefferman (1993): if the vorticity direction ξ = ω/|ω| is Lipschitz continuous in regions of intense vorticity (|ω| > Λ), then the N-S solution remains regular. Geometric criterion: vortex lines cannot form 'kinks' too rapidly. Singularity formation requires vorticity directions to become geometrically incoherent on the set where |ω| is large. Novel connection: potential singularities require BOTH Hölder exponent h < 1/3 (Onsager) AND angular incoherence (Constantin-Fefferman) concentrated on a set of Hausdorff dimension → 0 — possibly a self-contradictory geometric condition.",
      refs: ["Constantin & Fefferman 1993"]
    },
    {
      id: "geometric_cascade_mechanism",
      label: "Geometric Mechanism of Cascade (Open)",
      type: "open_problem",
      angle: "geometry_topology",
      importance: 9,
      detail: "A precise geometric description of WHY energy cascades forward in 3D turbulence remains lacking. The vortex stretching narrative (ω is amplified by (ω·∇)u) is physical but not rigorous. DNS shows preferential alignment of vorticity with the intermediate eigenvector of the strain tensor (Betchov relations), but why this specific geometric alignment drives forward cascade (not backward) in 3D is not rigorously established. The connection to Arnold's negative-curvature SDiff geometry is suggestive but not quantitative.",
    },
    {
      id: "gaussian_multiplicative_chaos",
      label: "Gaussian Multiplicative Chaos (GMC)",
      type: "proven",
      angle: "cross_domain",
      importance: 9,
      year: 1985,
      detail: "Kahane (1985), Rhodes-Vargas (2010s): for a log-correlated Gaussian field φ(x) with ⟨φ(x)φ(y)⟩ ~ −log|x−y|, the measure M_γ = e^{γφ(x) − γ²⟨φ²⟩/2}dx is rigorously well-defined for γ < √2 (subcritical). Energy dissipation in turbulence is conjectured to be in the GMC universality class. GMC gives ζ_p = p − γ²p²/2 (log-normal multifractal). Critical GMC (γ = √2) connects to KPZ/Liouville theory. This is the most mathematically developed rigorous framework touching turbulence intermittency.",
      refs: ["Kahane 1985", "Rhodes & Vargas 2014"]
    },
    {
      id: "ldt_multifractal_rigorous",
      label: "Multifractal = Large Deviations (Rigorous)",
      type: "proven",
      angle: "cross_domain",
      importance: 9,
      detail: "The equivalence between multifractal scaling and large deviations theory is mathematically rigorous for multiplicative cascade models. The multifractal spectrum D(h) IS the Cramér rate function of the cascade multiplier distribution, connected via Legendre transform. For GMC, this is an identity. For log-Poisson (She-Lévêque), the Lévy multiplier distribution gives D(h) via non-Gaussian Cramér theory. This is the deepest rigorous cross-domain connection in turbulence — not an analogy but an exact equivalence in the cascade framework.",
      refs: ["Kahane & Peyrière 1976", "Barral & Mandelbrot 2002"]
    },
    {
      id: "instanton_turbulence",
      label: "Instantons for Extreme Events",
      type: "framework",
      angle: "rg_qft",
      importance: 8,
      year: 1996,
      detail: "In the MSR path integral, instantons are saddle points dominating the path integral for rare/extreme events. Falkovich, Kolokolov, Lebedev, Migdal (1996): computed instantons for Burgers turbulence. For 3D N-S, the instanton equation is a coupled time-reversed Euler equation for the optimal fluctuation leading to extreme vorticity. Grafke, Grauer et al. (2015+) computed these numerically. This is the only systematic nonperturbative method for turbulence extremes — it directly computes the large-deviation rate function for extreme events from the N-S action.",
      refs: ["Falkovich et al. 1996", "Grafke et al. 2015"]
    },
    {
      id: "gallavotti_cohen",
      label: "Gallavotti-Cohen Fluctuation Theorem",
      type: "framework",
      angle: "cross_domain",
      importance: 6,
      year: 1995,
      detail: "For a chaotic system with time-reversible dynamics and SRB measure, the probability ratio of entropy production σ satisfies P(σ=A)/P(σ=−A) = e^{At} for large t. For turbulence: entropy production rate is ε; GC predicts a specific symmetry in the PDF of time-averaged dissipation fluctuations. Gallavotti (2000s) explicitly proposed this for turbulence. Testing it requires extremely long DNS runs to access backward-flux tail statistics. Provides weak inequality constraints on dissipation PDF, not the ζ_p values.",
      refs: ["Gallavotti & Cohen 1995"]
    },
    {
      id: "directed_percolation_transition",
      label: "Directed Percolation at Laminar-Turbulent Transition",
      type: "empirical",
      angle: "cross_domain",
      importance: 7,
      year: 2016,
      detail: "Lemoult et al. (2016, Nature Physics): the laminar-turbulent transition in quasi-1D channel flow belongs to the directed percolation (DP) universality class — experimentally verified with DP exponents. Turbulent puffs spread and decay with DP-like critical exponents. This is a genuine application of critical phenomena universality, but it concerns the ONSET of turbulence, not the statistical structure of fully developed turbulence. The inertial range scaling is unrelated to this DP transition.",
      refs: ["Lemoult et al. 2016"]
    },
    {
      id: "levy_stable_cascade",
      label: "Lévy-Stable Cascade Family",
      type: "novel",
      angle: "cross_domain",
      importance: 8,
      detail: "The multiplicative cascade can be parametrized by the Lévy-Khintchine representation: log ε_r = ∫ dΛ where Λ is a Lévy process in log-scale. Gaussian component → log-normal (K62). Poisson component → log-Poisson (She-Lévêque). Intermediate stable distributions → family interpolating between them. The Lévy exponent β directly maps to the multifractal spectrum via Legendre transform. Memory (non-Markovian) extensions → fractional Lévy processes. This framework unifies K62, She-Lévêque, and GMC as special cases, and suggests a Lévy exponent measurement from DNS as a new way to characterize turbulence.",
    },
    {
      id: "muzy_bacry_mrw",
      label: "Muzy-Bacry Multifractal Random Walk",
      type: "framework",
      angle: "cross_domain",
      importance: 7,
      year: 2003,
      detail: "Bacry and Muzy (2003): multifractal random walk with log-correlated Gaussian volatility ⟨ω(t)ω(s)⟩ = λ² log(T/|t−s|). Gives exact ζ_p = p/2 − λ²p²/8. Extended with fractional kernel ⟨ω(t)ω(s)⟩ ~ λ²|t−s|^{2H−1} encodes cascade memory via Hurst exponent H. This is the concrete mathematical framework for non-Markovian cascade memory in 1D that has not been fully extended to 3D turbulence structure functions — the fractional extension would directly test whether cascade memory explains anomalous scaling.",
      refs: ["Bacry & Muzy 2003"]
    },
    {
      id: "hopf_functional_equation",
      label: "Hopf Functional Equation (Exact, Unsolvable)",
      type: "framework",
      angle: "rg_qft",
      importance: 7,
      year: 1952,
      detail: "Hopf (1952): the exact equation for the characteristic functional Z[J] = ⟨exp(i∫J·u dx)⟩ of the turbulent velocity field. This is the analog of the Schrödinger equation for turbulence — linear in Z but with functional derivatives. All turbulence closures can be viewed as approximations to Hopf's equation. The MSR generating functional IS Hopf's functional written as a path integral. In principle, solving Hopf gives everything; in practice, it is as unsolvable as N-S.",
      refs: ["Hopf 1952"]
    },
    {
      id: "ot_benamou_brenier",
      label: "Benamou-Brenier Formulation in Scale Space",
      type: "novel",
      angle: "cross_domain",
      importance: 8,
      year: 2000,
      detail: "NOVEL ABSTRACTION: The Benamou-Brenier (2000) dynamical OT: W_2²(μ₀,μ₁) = inf_{ρ,v} ∫₀¹∫|v|²ρ dx dt subject to ∂_t ρ + ∇·(ρv) = 0. In turbulence scale space: μ_s = energy density at log-scale s = log(L/r), v_s = energy transfer rate. Optimality conditions → Hamilton-Jacobi equation in scale space. This makes explicit the sense in which the cascade is 'optimal' — deviations from optimality quantify intermittency 'waste'. The connection between non-geodesic transport and entropy production could tie OT to the Gallavotti-Cohen theorem. Completely unexplored in turbulence literature.",
      refs: ["Benamou & Brenier 2000"]
    },
    {
      id: "kraichnan_2d_duality",
      label: "Kraichnan 2D Inverse Cascade / Enstrophy",
      type: "proven",
      angle: "cross_domain",
      importance: 7,
      year: 1967,
      detail: "Kraichnan (1967): 2D turbulence has two conserved quantities — energy E = ∫|u|²/2 dx and enstrophy Z = ∫|ω|²/2 dx. This produces a dual cascade: energy flows upscale (inverse cascade, E(k) ~ k^{−5/3}) and enstrophy downscale (E(k) ~ k^{−3}). This is better understood than 3D turbulence — the additional conservation law provides the extra constraint that closes the cascade problem at leading order. The lesson: a second conserved invariant dramatically constrains the cascade. In 3D, helicity is a candidate second invariant but its role in constraining ζ_p is unclear.",
      refs: ["Kraichnan 1967"]
    },
    {
      id: "order_parameter_turbulence",
      label: "No Order Parameter for Turbulence",
      type: "open_problem",
      angle: "cross_domain",
      importance: 7,
      detail: "What is the order parameter for the 'turbulent state' analogous to magnetization in ferromagnetism? The Reynolds stress ⟨u_i u_j⟩ or energy spectrum E(k) are candidates but neither transforms like a Landau order parameter under symmetry breaking. The laminar-turbulent transition in pipe flow is subcritical (first-order like, no continuous bifurcation) making standard Landau theory inapplicable. Directed percolation works at the transition. For fully developed turbulence, there is no agreed order parameter — the field lacks a Landau-Ginzburg description.",
    },
    {
      id: "absolute_equilibrium_kraichnan",
      label: "Absolute Equilibrium of Truncated Euler",
      type: "framework",
      angle: "cross_domain",
      importance: 6,
      year: 1973,
      detail: "Kraichnan (1973): the Galerkin-truncated Euler equations have an exact equilibrium Gibbs measure P ~ exp(−αE − βH) where E = energy, H = helicity (3D) or enstrophy (2D). Predicts equipartition E(k) ~ k² (Rayleigh-Jeans). The transition from this equilibrium to the non-equilibrium cascade as viscosity → 0 (with k_max → ∞) is a genuine non-equilibrium phase transition whose universality class is unclassified. This transition might be where the RG fixed point governing turbulence 'emerges'.",
      refs: ["Kraichnan 1973"]
    },

  ];

  const extra_edges = [
    // New connections
    { source: "k41_hypotheses", target: "k41_dimensional_analysis", type: "derives_from", label: "assumptions enable dim. analysis" },
    { source: "universality_violations", target: "universality_assumption", type: "contradicts", label: "C_K varies 30% empirically" },
    { source: "closure_underdetermination", target: "closure_problem", type: "extends", label: "proves not just hard but structurally impossible" },
    { source: "sweeping_decorrelation", target: "locality_assumption", type: "contradicts", label: "large scales sweep small non-locally" },
    { source: "sweeping_decorrelation", target: "markovian_cascade", type: "critiques", label: "physical mechanism of memory" },
    { source: "backscatter_magnitude", target: "forward_cascade_assumption", type: "contradicts", label: "30-40% backscatter in spacetime" },
    { source: "re_scaling_inertial_range", target: "inertial_subrange", type: "critiques", label: "clean range needs Re > 3500" },
    { source: "helicity_cascade_dual", target: "helicity", type: "extends", label: "dual cascade structure" },
    { source: "helicity_cascade_dual", target: "k41_hypotheses", type: "critiques", label: "K41 ignores helicity dynamics" },
    { source: "sle_percolation_universality", target: "sle_2d_turbulence", type: "derives_from", label: "identifies universality class" },
    { source: "sle_percolation_universality", target: "novel_sle_3d", type: "enables", label: "if 2D=percolation, what is 3D?" },
    { source: "knotted_vortex_experiment", target: "helicity", type: "extends", label: "writhe→twist conversion observed" },
    { source: "knotted_vortex_experiment", target: "novel_topological_qft", type: "critiques", label: "knots transient; higher invariants don't persist" },
    { source: "beltrami_fields", target: "novel_topological_qft", type: "enables", label: "any knot realizable in Euler" },
    { source: "constantin_fefferman_criterion", target: "millennium_prize", type: "enables", label: "geometric path to regularity" },
    { source: "constantin_fefferman_criterion", target: "onsager_conjecture", type: "analogous_to", label: "both give geometric thresholds for singularity" },
    { source: "geometric_cascade_mechanism", target: "k41_hypotheses", type: "open_question", label: "K41 assumes but doesn't explain forward cascade" },
    { source: "gaussian_multiplicative_chaos", target: "multifractal", type: "derives_from", label: "GMC is rigorous foundation" },
    { source: "gaussian_multiplicative_chaos", target: "ldt_multifractal_rigorous", type: "derives_from", label: "GMC makes LDT=multifractal rigorous" },
    { source: "ldt_multifractal_rigorous", target: "anomalous_scaling", type: "enables", label: "rigorous framework for ζ_p computation" },
    { source: "ldt_multifractal_rigorous", target: "large_deviations", type: "derives_from", label: "the identity is exact in cascade models" },
    { source: "levy_stable_cascade", target: "log_normal", type: "extends", label: "log-normal is Gaussian special case" },
    { source: "levy_stable_cascade", target: "she_leveque", type: "extends", label: "log-Poisson is Lévy special case" },
    { source: "levy_stable_cascade", target: "gaussian_multiplicative_chaos", type: "extends", label: "GMC is Gaussian Lévy case" },
    { source: "muzy_bacry_mrw", target: "novel_nonmarkovian", type: "extends", label: "concrete 1D model of cascade memory" },
    { source: "muzy_bacry_mrw", target: "gaussian_multiplicative_chaos", type: "extends", label: "MRW is temporal GMC" },
    { source: "instanton_turbulence", target: "msr_field_theory", type: "derives_from", label: "instantons are MSR saddle points" },
    { source: "instanton_turbulence", target: "large_deviations", type: "enables", label: "instantons compute rate functions" },
    { source: "instanton_turbulence", target: "anomalous_scaling", type: "open_question", label: "can instanton ensemble give ζ_p?" },
    { source: "hopf_functional_equation", target: "navier_stokes", type: "derives_from", label: "exact reformulation" },
    { source: "hopf_functional_equation", target: "msr_field_theory", type: "analogous_to", label: "MSR = Hopf as path integral" },
    { source: "hopf_functional_equation", target: "closure_problem", type: "derives_from", label: "moment hierarchy is Hopf expanded" },
    { source: "directed_percolation_transition", target: "phase_transition_analogy", type: "extends", label: "DP is confirmed for transition; not inertial range" },
    { source: "gallavotti_cohen", target: "non_equilibrium_thermo", type: "extends", label: "GC is steady-state fluctuation theorem" },
    { source: "gallavotti_cohen", target: "novel_arrow_of_time", type: "enables", label: "GC constrains PDF of entropy production" },
    { source: "kraichnan_2d_duality", target: "helicity_cascade_dual", type: "analogous_to", label: "2D has enstrophy; 3D has helicity as 2nd invariant" },
    { source: "order_parameter_turbulence", target: "phase_transition_analogy", type: "critiques", label: "no order parameter → Landau theory inapplicable" },
    { source: "absolute_equilibrium_kraichnan", target: "rg_fixed_point", type: "open_question", label: "does cascade fixed point emerge from equilibrium transition?" },
    { source: "ot_benamou_brenier", target: "novel_optimal_transport", type: "extends", label: "dynamical OT formulation in scale space" },
    { source: "ot_benamou_brenier", target: "gallavotti_cohen", type: "open_question", label: "OT waste = entropy production?" },
    { source: "kolmogorov_constant_problem", target: "universality_violations", type: "derives_from", label: "C_K variation is the measurement" },
    { source: "novel_nonmarkovian", target: "sweeping_decorrelation", type: "derives_from", label: "sweeping is the physical memory mechanism" },
    // RG/QFT agent additions
    { source: "msr_action_detail", target: "msr_field_theory", type: "extends", label: "exact action specification" },
    { source: "nprg_fixed_point", target: "rg_fixed_point", type: "extends", label: "nonperturbative confirmation" },
    { source: "nprg_fixed_point", target: "anomalous_scaling", type: "open_question", label: "non-scale-invariant FP = intermittency mechanism" },
    { source: "perturbative_obstruction", target: "anomalous_scaling", type: "blocks", label: "no small parameter to compute ζ_p" },
    { source: "zero_modes_mechanism", target: "novel_nonmarkovian", type: "enables", label: "zero modes = large-scale memory at small scales" },
    { source: "logarithmic_cft_conjecture", target: "rg_fixed_point", type: "extends", label: "FP may be logarithmic/non-diagonalizable CFT" },
    { source: "logarithmic_cft_conjecture", target: "anomalous_scaling", type: "open_question", label: "Jordan block mixing → multiscaling?" },
    { source: "four_fifths_law_exact", target: "karman_howarth", type: "derives_from", label: "exact from N-S + stationarity" },
    { source: "kraichnan_ns_gap", target: "zero_modes", type: "open_question", label: "do NS zero modes exist?" },
    { source: "ope_fusion_rules_detail", target: "ope_turbulence", type: "extends", label: "L'vov-Procaccia + Falkovich-Zamolodchikov" },
    { source: "nprg_fixed_point", target: "novel_info_geometry", type: "enables", label: "NPRG flow = geodesic on theory space" },
  ];

  // ── From RG/QFT Agent ────────────────────────────────────────────────────
  extra_nodes.push(
    {
      id: "nprg_fixed_point",
      label: "NPRG Fixed Point (Canet-Delamotte-Wschebor)",
      type: "proven",
      angle: "rg_qft",
      importance: 9,
      year: 2016,
      detail: "Canet, Delamotte, Wschebor (2011–2022): using the Wetterinck effective average action (nonperturbative RG), a fixed point of the 3D NS turbulence RG flow is confirmed to exist. Crucially: this fixed point does NOT have classical scale invariance — the fixed-point effective action Γ* has non-trivial momentum dependence. This non-scale-invariance IS the QFT mechanism for intermittency. 2022 JFM work derived exact multi-time multi-point Eulerian correlation functions at large k. ζ_p values remain outside current numerical precision.",
      refs: ["Canet et al. 2016", "Canet et al. 2022"]
    },
    {
      id: "perturbative_obstruction",
      label: "Structural Obstruction: No Small Parameter at g* ~ O(1)",
      type: "proven",
      angle: "rg_qft",
      importance: 10,
      detail: "Physical 3D turbulence requires ε = 2 in the force-exponent expansion, giving effective coupling g* ~ O(1) at the Kolmogorov fixed point. No small parameter exists to organize perturbation theory. The anomalous exponents ζ_p − p/3 are genuinely non-perturbative in g. The composite operators governing ζ_p have NEGATIVE naive scaling dimensions — they are 'infrared-relevant' and contaminate all loop orders simultaneously. This is not a technical limitation: the anomalous exponents are inherently non-perturbative quantities. This is the deepest mathematical statement about WHY ζ_p has not been computed.",
    },
    {
      id: "logarithmic_cft_conjecture",
      label: "Turbulence as Logarithmic CFT (Novel)",
      type: "novel",
      angle: "rg_qft",
      importance: 8,
      detail: "NOVEL ABSTRACTION: The NPRG result that the NS fixed point lacks classical scale invariance is structurally similar to logarithmic CFTs, where the dilatation operator has Jordan-block structure (non-diagonalizable). In log-CFTs, operator mixing produces log(r) corrections to power-law scaling — similar to intermittency corrections. If turbulence IS a logarithmic CFT, then modern bootstrap methods for non-unitary/logarithmic CFTs (Hogervorst-Rychkov-van Rees program) become applicable. This would be an entirely new line of attack: no paper has explicitly framed turbulence anomalous scaling as log-CFT.",
    },
    {
      id: "msr_action_detail",
      label: "MSR Action: S[v,ṽ] Explicit Form",
      type: "framework",
      angle: "rg_qft",
      importance: 7,
      detail: "S[v,ṽ] = ∫ dt d³x { ṽᵢ[∂ₜvᵢ + (v·∇)vᵢ + ∇ᵢp − ν∇²vᵢ] − (1/2)ṽᵢDᵢⱼṽⱼ } where ṽᵢ is the auxiliary response field and Dᵢⱼ(k) ∝ k^{4−d−2ε} is the force correlator. Cubic vertex ṽᵢvⱼ∂ⱼvᵢ is the interaction. Galilean invariance → Ward identity constraining only zero-mode of vertex. BRST-like symmetry encodes causality. The action is real (unlike Euclidean QFT) and the coupling g ∝ (force amplitude)^{1/2}. Physical 3D turbulence: ε = 2, giving g* ~ O(1).",
    },
    {
      id: "kraichnan_ns_gap",
      label: "Kraichnan-to-NS Transfer Gap",
      type: "open_problem",
      angle: "rg_qft",
      importance: 10,
      detail: "Kraichnan model succeeds because: (1) velocity is prescribed (no self-advection), (2) delta-in-time correlation eliminates memory and closes the hierarchy, (3) no pressure. In NS: velocity self-advects (nonlinear feedback), pressure enforces incompressibility nonlocally, velocity is non-Gaussian and has time correlations. Attempts to treat NS velocity as approximately Gaussian fail because pressure and non-Gaussianity are O(1), not small corrections. There is NO known way to close the NS n-point hierarchy without uncontrolled approximations. The zero-mode mechanism is non-transferable in any rigorous sense.",
    },
    {
      id: "ope_fusion_rules_detail",
      label: "Fusion Rules / OPE for NS Turbulence",
      type: "framework",
      angle: "rg_qft",
      importance: 7,
      year: 1996,
      detail: "L'vov and Procaccia (1996): as two spatial arguments coalesce, n-point NS correlation functions factorize: F_n(r₁₂ → 0, R) ~ (r₁₂/R)^{ζ_p − ζ_q} F_m(R). Falkovich-Zamolodchikov (2015): genuine OPE language applied to energy cascade correlations. Fusion rules confirmed experimentally at first order. Fundamental limitation: OPE gives consistency relations between exponents but cannot compute them — the operator dimensions ζ_p are inputs, not outputs. Analogous situation to 2D CFT where OPE is powerful given the Virasoro algebra, but for turbulence the analog of the Virasoro algebra is unknown.",
      refs: ["L'vov & Procaccia 1996", "Falkovich & Zamolodchikov 2015"]
    }
  );

  // Merge into GRAPH_DATA
  const existingIds = new Set(GRAPH_DATA.nodes.map(n => n.id));
  extra_nodes.forEach(n => {
    if (!existingIds.has(n.id)) {
      GRAPH_DATA.nodes.push(n);
      existingIds.add(n.id);
    }
  });

  const existingEdges = new Set(GRAPH_DATA.edges.map(e => `${e.source}->${e.target}`));
  extra_edges.forEach(e => {
    const key = `${e.source}->${e.target}`;
    if (!existingEdges.has(key)) {
      GRAPH_DATA.edges.push(e);
      existingEdges.add(key);
    }
  });

})();
