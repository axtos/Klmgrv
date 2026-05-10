// Turbulence Knowledge Graph — Kolmogorov & Beyond
// Node types: proven | empirical | conjecture | framework | open_problem | analogy | novel | historical | assumption
// Edge types: derives_from | contradicts | extends | analogous_to | enables | blocks | assumption_of | critiques | open_question

const GRAPH_DATA = {
  nodes: [

    // ── HISTORICAL ──────────────────────────────────────────────────────────
    {
      id: "navier_stokes",
      label: "Navier-Stokes Equations",
      type: "framework",
      angle: "core",
      importance: 10,
      year: 1845,
      detail: "∂u/∂t + (u·∇)u = −∇p/ρ + ν∇²u, ∇·u = 0. The fundamental PDE for incompressible viscous flow. Derived by Navier (1822) and Stokes (1845). Existence and smoothness of smooth 3D solutions for all time is unsolved — one of the seven Millennium Prize Problems. Every other node in this graph depends on this foundation.",
      refs: ["Navier 1822", "Stokes 1845"]
    },
    {
      id: "reynolds_number",
      label: "Reynolds Number Re",
      type: "framework",
      angle: "core",
      importance: 9,
      year: 1883,
      detail: "Re = UL/ν — ratio of inertial to viscous forces. At low Re flow is laminar; turbulence emerges beyond a critical Re (varies by geometry, typically ~1000–10000). The cost of DNS scales as Re^(9/4) per spatial dimension, making high-Re simulation computationally prohibitive. Engineering flows reach Re ~ 10^7–10^9.",
      refs: ["Reynolds 1883"]
    },
    {
      id: "richardson_cascade",
      label: "Richardson Energy Cascade",
      type: "historical",
      angle: "core",
      importance: 8,
      year: 1922,
      detail: "'Big whorls have little whorls...' Richardson (1922) proposed the qualitative picture of energy flowing from large eddies to successively smaller ones until viscosity dissipates it as heat. This is the conceptual backbone of K41 but was entirely qualitative — no quantitative predictions.",
      refs: ["Richardson 1922"]
    },

    // ── K41 CORE ──────────────────────────────────────────────────────────
    {
      id: "k41_hypotheses",
      label: "K41 Similarity Hypotheses",
      type: "framework",
      angle: "core",
      importance: 10,
      year: 1941,
      detail: "Kolmogorov 1941: (1) At sufficiently small scales, turbulence is locally isotropic. (2) In the equilibrium range, statistics depend only on ε (dissipation rate) and ν (viscosity). (3) In the inertial subrange (η ≪ r ≪ L), statistics depend only on ε and r. The power of this: dimensional analysis alone gives quantitative predictions. The weakness: assumes ε is a single global constant — contradicted by intermittency.",
      refs: ["Kolmogorov 1941a", "Kolmogorov 1941b"]
    },
    {
      id: "kolmogorov_scales",
      label: "Kolmogorov Microscales",
      type: "proven",
      angle: "core",
      importance: 9,
      year: 1941,
      detail: "From ε and ν by dimensional analysis: η = (ν³/ε)^(1/4) (length), τ_η = (ν/ε)^(1/2) (time), u_η = (νε)^(1/4) (velocity). These are exact — no free parameters. They define the scale where viscosity wins. For aircraft boundary layer: η ~ 0.1mm. For DNS you must resolve all scales from domain size L down to η, giving (L/η)^3 ~ Re^(9/4) grid points.",
      refs: ["Kolmogorov 1941a"]
    },
    {
      id: "minus_five_thirds",
      label: "−5/3 Energy Spectrum",
      type: "empirical",
      angle: "core",
      importance: 10,
      year: 1941,
      detail: "E(k) = C_K · ε^(2/3) · k^(−5/3) in the inertial subrange. Derived purely by dimensional analysis from ε. Confirmed experimentally across atmospheres, oceans, wind tunnels, pipe flows. Kolmogorov constant C_K ≈ 1.5 (empirical, not derived). The −5/3 exponent is equivalent to ζ_2 = 2/3 for second-order structure functions. Robust but does not uniquely determine higher-order statistics.",
      refs: ["Kolmogorov 1941b", "Grant et al. 1962"]
    },
    {
      id: "four_fifths_law",
      label: "Kolmogorov 4/5 Law",
      type: "proven",
      angle: "core",
      importance: 10,
      year: 1941,
      detail: "S_3(r) = ⟨(δu_L)³⟩ = −(4/5)εr in the inertial range. This is the ONLY non-trivial exact result derived directly from Navier-Stokes for turbulence. It requires no closure assumptions and no similarity hypotheses — it follows from the Kármán-Howarth equation in the inertial limit. All other turbulence scaling results are approximate or empirical. The exactness of this result constrains ζ_3 = 1 exactly.",
      refs: ["Kolmogorov 1941c", "Monin & Yaglom 1975"]
    },
    {
      id: "structure_functions",
      label: "Structure Functions S_p(r)",
      type: "framework",
      angle: "core",
      importance: 9,
      detail: "S_p(r) = ⟨|u(x+r) − u(x)|^p⟩ ~ r^ζ_p. K41 predicts ζ_p = p/3. Experiments show: ζ_2 ≈ 0.696 (close to 2/3), ζ_3 = 1 (exact), ζ_4 ≈ 1.28 (K41 predicts 4/3 ≈ 1.33), ζ_6 ≈ 1.78 (K41 predicts 2.0). The deviation from p/3 grows with p — this is anomalous scaling. Computing ζ_p from N-S for p ≠ 2, 3 is an open problem.",
      refs: ["Frisch 1995"]
    },
    {
      id: "inertial_subrange",
      label: "Inertial Subrange",
      type: "empirical",
      angle: "core",
      importance: 8,
      detail: "The range of scales η ≪ r ≪ L where neither viscosity nor forcing directly act. Energy cascades through this range without creation or dissipation. Requires Re ≫ 1 to exist. In practice, a clean inertial range requires Re > ~10^4 and even then spans only a few decades. At finite Re there is contamination from both ends — the 'inertial range' is somewhat idealized.",
    },
    {
      id: "dissipation_rate",
      label: "Energy Dissipation Rate ε",
      type: "framework",
      angle: "core",
      importance: 9,
      detail: "ε = ν ⟨(∂u_i/∂x_j)²⟩ — rate at which kinetic energy converts to heat via viscosity. K41 treats this as a single constant characterizing the flow. K62 recognizes ε_r (spatially coarse-grained) fluctuates. The spatial distribution of ε is highly intermittent — concentrated in thin sheets and filaments (vortex sheets, worm-like vortex tubes) filling only a small fraction of the volume.",
    },

    // ── K62 / INTERMITTENCY ──────────────────────────────────────────────
    {
      id: "k62",
      label: "K62 Refined Similarity",
      type: "framework",
      angle: "core",
      importance: 9,
      year: 1962,
      detail: "Kolmogorov 1962 — response to Landau's criticism that ε fluctuates. Replaces global ε with local average ε_r over ball of radius r. Refined hypothesis: statistics at scale r depend on ε_r and r. K62 predicts ζ_p = p/3 − μp(p−3)/18 with μ ≈ 0.25 (intermittency exponent). Improves on K41 but still doesn't derive the corrections from first principles — μ is empirical.",
      refs: ["Kolmogorov 1962", "Oboukhov 1962"]
    },
    {
      id: "intermittency",
      label: "Intermittency",
      type: "empirical",
      angle: "core",
      importance: 10,
      detail: "Turbulent dissipation is not uniformly distributed — it concentrates in regions of intense vorticity (vortex tubes, sheets) separated by relatively quiet regions. Volume fraction of intense regions decreases with scale. This spatial organization is the root cause of anomalous scaling. K41's key failure: it assumed ε is a global constant when it's actually a spatially fluctuating random field with long-range correlations.",
      refs: ["Batchelor & Townsend 1949", "Frisch 1995"]
    },
    {
      id: "anomalous_scaling",
      label: "Anomalous Scaling / ζ_p Problem",
      type: "open_problem",
      angle: "core",
      importance: 10,
      detail: "ζ_p ≠ p/3 for p > 3 (and p < 3). The deviation from linear scaling in p is the central unsolved quantitative problem of turbulence theory. K41 predicts ζ_p = p/3 exactly. Experiments give ζ_6 ≈ 1.78 vs K41's 2.0. Nobody has derived ζ_p from N-S for p ≠ 2, 3. The multifractal model parametrizes these exponents geometrically but does not derive them.",
    },
    {
      id: "landau_criticism",
      label: "Landau's Intermittency Criticism",
      type: "historical",
      angle: "core",
      importance: 7,
      year: 1944,
      detail: "Landau pointed out (privately, ~1944, published in Fluid Mechanics textbook) that K41 fails because ε is not a universal constant — it fluctuates depending on the large-scale flow realization. The Kolmogorov constant C_K must then fluctuate too, breaking universality. This was the seed of what became intermittency theory. Kolmogorov took 18 years to respond with K62.",
      refs: ["Landau & Lifshitz 1959"]
    },

    // ── MULTIFRACTAL ────────────────────────────────────────────────────
    {
      id: "multifractal",
      label: "Multifractal Formalism",
      type: "framework",
      angle: "geometry_topology",
      importance: 9,
      year: 1985,
      detail: "Parisi & Frisch (1985): turbulent velocity field has Hölder exponent h at each point x, with the set of points where velocity is h-Hölder having Hausdorff dimension D(h). The singularity spectrum f(α) = D(h)|_{h=α}. Structure function exponents: ζ_p = min_h [ph + 3 − D(h)] (Legendre transform). Consistent with anomalous scaling empirically but D(h) must be measured, not derived. The framework is purely descriptive.",
      refs: ["Parisi & Frisch 1985", "Frisch 1995"]
    },
    {
      id: "log_normal",
      label: "Log-Normal Distribution of ε_r",
      type: "empirical",
      angle: "core",
      importance: 7,
      year: 1962,
      detail: "K62 proposes ln(ε_r) is normally distributed with variance σ² = A + μ ln(L/r). Equivalently, ε_r is log-normally distributed. Empirically: better than constant ε but deviations visible at high orders. Log-Poisson model (She-Lévêque 1994) fits ζ_p better than log-normal. Neither is derived from N-S — both are phenomenological ansätze for the distribution of ε_r.",
      refs: ["Kolmogorov 1962", "She & Lévêque 1994"]
    },
    {
      id: "she_leveque",
      label: "She-Lévêque Log-Poisson Model",
      type: "empirical",
      angle: "core",
      importance: 7,
      year: 1994,
      detail: "ζ_p = p/9 + 2[1 − (2/3)^(p/3)]. Fits experimental ζ_p remarkably well. Based on the assumption that the most intense dissipative structures are 1D (filaments), giving a log-Poisson distribution for ε_r. Better fit than log-normal but again purely phenomenological. The model has structural assumptions (1D filaments as 'most singular') that haven't been derived from N-S.",
      refs: ["She & Lévêque 1994"]
    },

    // ── MATHEMATICAL/PDE ─────────────────────────────────────────────────
    {
      id: "millennium_prize",
      label: "Millennium Prize: N-S Existence & Smoothness",
      type: "open_problem",
      angle: "core",
      importance: 9,
      year: 2000,
      detail: "Do smooth solutions to 3D N-S exist for all time given smooth initial data? In 2D: yes (Leray 1934 for weak solutions; smooth solutions exist). In 3D: only weak solutions (Leray-Hopf) guaranteed. Smooth solutions could develop a finite-time singularity (blow-up). If singularity occurs, N-S loses meaning. The question is whether real turbulence is 'tamed' by the equations or whether the equations themselves break down.",
      refs: ["Fefferman 2000", "Leray 1934"]
    },
    {
      id: "leray_hopf",
      label: "Leray-Hopf Weak Solutions",
      type: "proven",
      angle: "core",
      importance: 7,
      year: 1934,
      detail: "Leray (1934) and Hopf (1951) proved existence of 'weak' solutions to 3D N-S for all time. These satisfy N-S in an integrated sense but may not be smooth — they can have discontinuities in space-time. Whether these weak solutions are unique is also open. The turbulent velocity field may be such a weak solution.",
      refs: ["Leray 1934", "Hopf 1951"]
    },
    {
      id: "onsager_conjecture",
      label: "Onsager Conjecture (Proved 2019)",
      type: "proven",
      angle: "core",
      importance: 8,
      year: 1949,
      detail: "Onsager (1949) conjectured: if the Hölder regularity h > 1/3, Euler equation conserves energy; if h ≤ 1/3, energy dissipation is possible even without viscosity ('anomalous dissipation'). Proved in two parts: conservation for h > 1/3 (Constantin et al. 1994); existence of solutions with h ≤ 1/3 that dissipate energy (Isett 2018, De Lellis & Székelyhidi). Relevant: K41 implies h = 1/3 at the threshold.",
      refs: ["Onsager 1949", "Isett 2018", "De Lellis & Székelyhidi 2019"]
    },
    {
      id: "karman_howarth",
      label: "Kármán-Howarth Equation",
      type: "proven",
      angle: "core",
      importance: 7,
      year: 1938,
      detail: "Exact evolution equation for two-point velocity correlations, derived from N-S. Contains third-order correlations — the first instance of the closure problem. The 4/5 law is derived from this in the inertial limit. Hierarchy: N-S → equation for 2-point stat contains 3-point stat → equation for 3-point stat contains 4-point stat → ... This infinite hierarchy is the closure problem in its purest form.",
      refs: ["Kármán & Howarth 1938"]
    },
    {
      id: "closure_problem",
      label: "The Closure Problem",
      type: "open_problem",
      angle: "core",
      importance: 10,
      detail: "The statistical moments of turbulence satisfy an infinite hierarchy: equation for ⟨u^n⟩ involves ⟨u^{n+1}⟩. There is no natural truncation. Every turbulence model (RANS k-ε, LES subgrid models) is a closure approximation. The mathematical reason: N-S is nonlinear, so the mean of a product ≠ product of means. This is not a technical obstacle — it is fundamental to nonlinear PDEs. It cannot be eliminated without additional physical input.",
    },

    // ── RG / FIELD THEORY ────────────────────────────────────────────────
    {
      id: "msr_field_theory",
      label: "MSR Field Theory of Turbulence",
      type: "framework",
      angle: "rg_qft",
      importance: 8,
      year: 1973,
      detail: "Martin-Siggia-Rose (1973) reformulated stochastic N-S as a path integral / field theory. The action S[u, p̂] couples velocity field u to response field p̂. In principle, all correlation functions computable as Feynman diagrams. In practice: the theory is strongly coupled and non-perturbative at Re ≫ 1. Perturbative expansion fails — coupling constant is effectively infinite at the fixed point governing the inertial range.",
      refs: ["Martin, Siggia & Rose 1973", "Wyld 1961"]
    },
    {
      id: "rg_fixed_point",
      label: "RG Fixed Point of Turbulence",
      type: "conjecture",
      angle: "rg_qft",
      importance: 9,
      detail: "The inertial range is scale-invariant → it should correspond to a fixed point of the renormalization group. This fixed point is what controls the anomalous scaling exponents ζ_p. It is known to exist (the −5/3 spectrum is robust), but its precise nature — what operators are relevant/irrelevant, what the anomalous dimensions are — cannot be computed by perturbative methods. The fixed point is strongly coupled.",
    },
    {
      id: "epsilon_expansion_failure",
      label: "ε-Expansion Failure in 3D",
      type: "empirical",
      angle: "rg_qft",
      importance: 7,
      detail: "In critical phenomena, the d = 4 − ε expansion organizes perturbation theory in small ε, giving good results for d=3 (ε=1). For turbulence, the analogous expansion around d=2 or d=4 fails because the coupling is not small at d=3. The turbulence RG fixed point is infrared, not ultraviolet, and strongly coupled. Perturbative RG gives the −5/3 law at leading order but cannot compute intermittency corrections.",
      refs: ["Forster, Nelson & Stephen 1977"]
    },
    {
      id: "kraichnan_model",
      label: "Kraichnan Model (Passive Scalar)",
      type: "proven",
      angle: "rg_qft",
      importance: 8,
      year: 1994,
      detail: "Kraichnan (1968, rigorous results from 1994): passive scalar advected by a Gaussian, white-in-time, self-similar velocity field. EXACT anomalous scaling exponents derivable via RG / exact diagonalization of a Schrödinger-like operator. The anomalous exponents arise from 'zero modes' of the operator. Key lesson: intermittency/anomalous scaling can be EXACT in a solvable model. The obstacle is extending this to active fields (where velocity feeds back on itself).",
      refs: ["Kraichnan 1968", "Gawedzki & Kupiainen 1995", "Falkovich, Gawedzki & Procaccia 2001"]
    },
    {
      id: "zero_modes",
      label: "Zero Modes / Statistically Conserved Quantities",
      type: "framework",
      angle: "rg_qft",
      importance: 8,
      detail: "In the Kraichnan model, anomalous scaling exponents are eigenvalues of an operator whose zero modes correspond to correlation functions that are statistically conserved despite the cascade. These zero modes SURVIVE in the limit of large scale separation — they are the fingerprint of the large-scale forcing that persists to small scales, breaking simple dimensional analysis. The key insight: anomalous scaling = survival of information about large-scale forcing at small scales.",
      refs: ["Gawedzki & Kupiainen 1995", "Bernard et al. 1998"]
    },
    {
      id: "ope_turbulence",
      label: "OPE Approach (Falkovich et al.)",
      type: "framework",
      angle: "rg_qft",
      importance: 7,
      detail: "Operator Product Expansion applied to turbulence correlation functions. In CFT, OPE gives exact recursion: O_i(x)O_j(y) = Σ_k C_{ijk}|x−y|^{Δ_k−Δ_i−Δ_j} O_k((x+y)/2). For turbulence, an analogous structure exists but the spectrum of operators and their dimensions Δ_k are unknown. OPE has been used to derive constraint equations on exponents but not to compute them.",
      refs: ["Falkovich, Gawedzki & Procaccia 2001"]
    },

    // ── GEOMETRY / TOPOLOGY ──────────────────────────────────────────────
    {
      id: "helicity",
      label: "Helicity as Topological Invariant",
      type: "proven",
      angle: "geometry_topology",
      importance: 7,
      year: 1969,
      detail: "H = ∫ u·ω dV where ω = ∇×u. Moffatt (1969): helicity measures the average linking number of vortex lines — a topological quantity. In ideal (inviscid) flow, helicity is exactly conserved. In turbulence with viscosity, helicity is not conserved but its rate of change involves triple correlations. Helicity constraint: flows with H≠0 have vortex lines that are topologically linked, which constrains reconnection events.",
      refs: ["Moffatt 1969"]
    },
    {
      id: "vortex_tubes",
      label: "Vortex Tubes / Worms",
      type: "empirical",
      angle: "geometry_topology",
      importance: 8,
      detail: "DNS and experiments reveal intense vorticity concentrated in tube-like structures ('worms') of diameter ~ η and length ~ L^{1/4}η^{3/4} (intermediate scale). These are the 'skeleton' of intermittency. Their geometry — Hausdorff dimension between 1 and 2, fractal cross-section — is not derived from theory. The She-Lévêque model assumes 1D filaments as most singular structures — partially but not fully supported by DNS.",
      refs: ["Siggia 1981", "Jiménez et al. 1993"]
    },
    {
      id: "vortex_reconnection",
      label: "Vortex Reconnection",
      type: "empirical",
      angle: "geometry_topology",
      importance: 7,
      detail: "Vortex tubes change topology through reconnection events — two tubes approach, break, and rejoin with different connectivity. This requires viscosity and changes helicity locally. Reconnection is suspected to drive the energy cascade (Kerr, Brenner) and may be responsible for finite-time singularity formation in Euler equations. The rate of reconnection and its statistical properties in turbulence are not analytically characterized.",
    },
    {
      id: "arnold_geodesic",
      label: "Arnold: Euler Flow as Geodesic",
      type: "framework",
      angle: "geometry_topology",
      importance: 7,
      year: 1966,
      detail: "Arnold (1966): ideal fluid flow is a geodesic on the infinite-dimensional Lie group SDiff(M) of volume-preserving diffeomorphisms, with the L² metric. Turbulence onset corresponds to geodesics becoming unstable (negative sectional curvature). Viscosity corresponds to adding a non-conservative force. This geometric picture is beautiful but has not yielded quantitative turbulence predictions — the geometry of SDiff is too complex.",
      refs: ["Arnold 1966", "Arnold & Khesin 1998"]
    },
    {
      id: "sle_2d_turbulence",
      label: "SLE(6) in 2D Turbulence",
      type: "empirical",
      angle: "geometry_topology",
      importance: 8,
      year: 2006,
      detail: "Bernard et al. (2006): boundaries of vorticity clusters in 2D inverse-cascade turbulence follow SLE(6) — the same process as percolation cluster boundaries. This is a remarkable connection: 2D turbulence at the inverse cascade regime is described by a conformal field theory with c = 0 (percolation CFT). SLE(6) implies the boundaries are non-differentiable curves with Hausdorff dimension 7/4. This has NOT been established for 3D turbulence.",
      refs: ["Bernard et al. 2006"]
    },
    {
      id: "conformal_2d",
      label: "Conformal Symmetry in 2D Turbulence",
      type: "empirical",
      angle: "rg_qft",
      importance: 7,
      detail: "2D turbulence (inverse energy cascade) exhibits approximate conformal invariance of vorticity level sets. This suggests a CFT description. The central charge c = 0 (Bernard et al.) corresponding to percolation. 3D turbulence is less symmetric — conformal invariance is broken by the forward cascade direction and anisotropic geometry. Whether 3D turbulence has hidden conformal symmetry is open.",
    },

    // ── CROSS-DOMAIN ANALOGIES ───────────────────────────────────────────
    {
      id: "phase_transition_analogy",
      label: "Turbulence ↔ Phase Transitions",
      type: "analogy",
      angle: "cross_domain",
      importance: 8,
      detail: "Critical phenomena (phase transitions) are also characterized by scale invariance and power laws. Wilson RG succeeded there by identifying the relevant order parameter and finding a weakly-coupled fixed point (d=4−ε expansion). For turbulence: (a) there is no obvious order parameter, (b) the fixed point is strongly coupled, (c) the symmetry breaking pattern is less clear. The analogy is structural but the tools don't transfer cleanly.",
    },
    {
      id: "large_deviations",
      label: "Large Deviations Theory",
      type: "framework",
      angle: "cross_domain",
      importance: 7,
      detail: "Large deviations theory characterizes the probability of rare, extreme events: P(X > x) ~ e^{−nI(x)} where I is the rate function. Applied to turbulence: probability of extreme dissipation events, velocity increments far in the tails. Multifractal formalism is essentially large deviations theory applied to scale-by-scale statistics. Provides a consistent framework for intermittency but doesn't derive the rate function from N-S.",
    },
    {
      id: "non_equilibrium_thermo",
      label: "Non-Equilibrium Thermodynamics",
      type: "analogy",
      angle: "cross_domain",
      importance: 6,
      detail: "Turbulence is a non-equilibrium steady state: energy injected at large scales, dissipated at small scales. Entropy production rate is exactly ε. Fluctuation theorems (Jarzynski, Crooks) constrain the statistics of work done in non-equilibrium processes. Applied to cascade: energy transfer events between scales are 'work' done by inertial forces. These constraints exist in principle but are weak — they give inequalities, not the exact ζ_p values.",
    },
    {
      id: "random_matrix",
      label: "Random Matrix Theory",
      type: "analogy",
      angle: "cross_domain",
      importance: 5,
      detail: "Velocity gradient tensor A_{ij} = ∂u_i/∂x_j has eigenvalue statistics that have been measured in DNS. Some universality with RMT has been observed for the strain-rate tensor. However, the velocity gradient tensor is not symmetric and the turbulent environment is correlated, so RMT universality is approximate and partial. Limited connection to the ζ_p problem.",
    },

    // ── COMPUTATIONAL / ML ───────────────────────────────────────────────
    {
      id: "dns",
      label: "Direct Numerical Simulation (DNS)",
      type: "framework",
      angle: "computational",
      importance: 9,
      detail: "DNS resolves all scales from L to η without modeling. Cost ~ Re^(9/4). Maximum DNS today: Re_λ ~ 2000 (Taylor-scale Reynolds number), corresponding to a few decades of inertial range. This is sufficient to test K41 but not to observe clean asymptotic scaling of high-order structure functions. DNS provides 'truth data' for testing theories and training ML models.",
    },
    {
      id: "neural_operators",
      label: "Neural Operators (FNO, DeepONet)",
      type: "framework",
      angle: "computational",
      importance: 7,
      year: 2021,
      detail: "Fourier Neural Operator (Li et al. 2021): learns the solution operator G: f ↦ u of a PDE, mapping forcing to solution in function space. Resolution-invariant, mesh-independent. Trained on DNS data, can predict turbulent fields in milliseconds. Limitation: accuracy degrades for Re outside training distribution; extrapolation to high Re is unreliable. Does not provide understanding of ζ_p.",
      refs: ["Li et al. 2021", "Lu et al. 2021"]
    },
    {
      id: "ml_turbulence_models",
      label: "ML-Augmented Turbulence Closure",
      type: "framework",
      angle: "computational",
      importance: 7,
      detail: "Neural networks trained on DNS data to correct RANS/LES closure errors. Approaches: Ling et al. (tensor-basis neural networks respecting Galilean invariance), Schmelzer et al. (sparse regression), Duraisamy et al. (field inversion + ML). Improves RANS accuracy for specific flow types but lacks generalizability and physical interpretability. Does not help with the fundamental ζ_p problem.",
    },

    // ── ASSUMPTIONS (CRITICAL) ────────────────────────────────────────────
    {
      id: "ergodicity_assumption",
      label: "Ergodicity Assumption",
      type: "assumption",
      angle: "assumptions",
      importance: 8,
      detail: "K41 implicitly assumes ergodicity: time averages equal ensemble averages. For stationary, homogeneous turbulence this is standard, but: (1) turbulence near walls is non-homogeneous, (2) turbulence near onset is non-stationary, (3) intermittent extreme events may be non-ergodic (they sample rare parts of phase space). If turbulence is weakly non-ergodic, ensemble statistics may not capture the physics of individual realizations.",
    },
    {
      id: "locality_assumption",
      label: "Locality of Energy Cascade",
      type: "assumption",
      angle: "assumptions",
      importance: 9,
      detail: "K41 assumes energy transfer is local in scale space — only neighboring triads (k, q, p with |k−q| ~ p) dominate. Non-local transfers (sweeping of small scales by large scales, distant triads) exist but are assumed secondary. DNS evidence: local transfers dominate but non-local contributions are ~20−30% of total transfer in some flows. If non-locality is fundamental, K41's scale-by-scale framework misses something important.",
    },
    {
      id: "markovian_cascade",
      label: "Markovian Cascade Assumption",
      type: "assumption",
      angle: "assumptions",
      importance: 9,
      detail: "NOVEL ASSUMPTION: K41 implicitly assumes the cascade is Markovian — statistics at scale r depend only on the statistics at scale r, not on the history of how the cascade arrived there from scale L. This means: the cascade has no 'memory' of large-scale structure beyond what's captured in ε_r. But if vortex topology (helicity, knottedness) imprints large-scale correlations that survive to small scales, the cascade has non-Markovian character. Zero modes in Kraichnan model suggest this is exactly what causes anomalous scaling.",
    },
    {
      id: "universality_assumption",
      label: "Small-Scale Universality",
      type: "assumption",
      angle: "assumptions",
      importance: 8,
      detail: "K41 Hypothesis 1: statistics at small scales are universal — independent of large-scale forcing, geometry, and boundary conditions. Evidence for: −5/3 spectrum appears across wildly different flows. Evidence against: C_K varies by ~10% across flow types; higher-order statistics show flow-dependent corrections; near-wall turbulence violates isotropy down to viscous scales. Universality may hold only for second-order statistics (energy spectrum) but not for higher-order/intermittency statistics.",
    },
    {
      id: "forward_cascade_assumption",
      label: "Forward Cascade Dominance",
      type: "assumption",
      angle: "assumptions",
      importance: 7,
      detail: "3D turbulence has a net forward energy cascade (large → small scales). But 'backscatter' — inverse transfer from small to large scales — exists and is ~10−30% of total energy flux in LES. In 2D, cascade is purely inverse (Kraichnan 1967). In quasi-2D (rotating turbulence, stratified flows, MHD), balance shifts. If backscatter is not just a correction but a fundamental feature of 3D cascade dynamics, then the standard picture is incomplete.",
    },

    // ── NOVEL ABSTRACTIONS ────────────────────────────────────────────────
    {
      id: "novel_nonmarkovian",
      label: "Non-Markovian Cascade Framework",
      type: "novel",
      angle: "novel",
      importance: 9,
      detail: "NOVEL ABSTRACTION: Reformulate the energy cascade as a non-Markovian stochastic process in scale space. Instead of: P(statistics at scale r | statistics at scale 2r), consider: P(statistics at scale r | statistics at all scales r' > r). This is a path in scale space with memory. Mathematical tools: fractional Brownian motion (Hurst exponent H ≠ 1/2 encodes memory), generalized Langevin equation with memory kernel, Volterra processes. Key question: can anomalous scaling exponents ζ_p be the Hurst exponents of this process? The Kraichnan zero-mode mechanism suggests yes — information from large scales survives.",
    },
    {
      id: "novel_sle_3d",
      label: "SLE for 3D Turbulence Structures",
      type: "novel",
      angle: "novel",
      importance: 8,
      detail: "NOVEL ABSTRACTION: Bernard et al. showed 2D turbulence vorticity cluster boundaries follow SLE(6). In 3D, vortex tube surfaces and dissipation set boundaries are 2D surfaces embedded in 3D. Could their intersection with a 2D plane follow an SLE process? If yes, this would: (1) give a conformal field theory description of 3D turbulence structure, (2) determine the fractal dimension of dissipation sets from first principles, (3) connect ζ_p to conformal weights. The κ parameter of SLE would encode the intermittency exponent μ.",
    },
    {
      id: "novel_info_geometry",
      label: "Information Geometry of Scale Space",
      type: "novel",
      angle: "novel",
      importance: 8,
      detail: "NOVEL ABSTRACTION: The space of turbulent flow statistics at scale r forms a statistical manifold with the Fisher-Rao metric g_{ij}(r). The RG flow (cascade) is a flow on this manifold. Fixed points are scale-invariant statistics (K41). Anomalous dimensions = geodesic curvature of RG trajectories on this manifold. Key question: is there a conserved quantity along the RG flow (like the c-theorem in 2D CFT, which says the central charge decreases monotonically under RG)? An analogue for 3D turbulence would constrain the accessible fixed points.",
    },
    {
      id: "novel_optimal_transport",
      label: "Optimal Transport in Scale Space",
      type: "novel",
      angle: "novel",
      importance: 7,
      detail: "NOVEL ABSTRACTION: The energy cascade moves energy 'mass' from the forcing scale L to the dissipation scale η in scale space. Optimal transport theory (Monge-Kantorovich problem) asks: what is the minimal-cost transport plan? If turbulence solves an optimal transport problem in scale space, the −5/3 spectrum would be the optimal transport 'density' and the Wasserstein distance would replace ε as the fundamental invariant. This reframes turbulence as a variational problem. Potential connection to gradient flows and dissipative structures.",
    },
    {
      id: "novel_topological_qft",
      label: "Topological Invariants via TQFT",
      type: "novel",
      angle: "novel",
      importance: 7,
      detail: "NOVEL ABSTRACTION: Vortex lines are 1D objects in 3D space — their topology is classified by knot invariants (Jones polynomial, HOMFLY). Helicity is the simplest such invariant. Could higher knot invariants (writhe, self-linking number, Vassiliev invariants) provide new conserved or slowly-varying quantities in turbulence? Topological Quantum Field Theory (Chern-Simons theory) computes these invariants for gauge fields. Vorticity ω is analogous to a gauge field (ω = ∇×u). A Chern-Simons theory of vorticity might give new constraints on turbulent statistics.",
    },
    {
      id: "novel_arrow_of_time",
      label: "Irreversibility as Fundamental Constraint",
      type: "novel",
      angle: "novel",
      importance: 8,
      detail: "NOVEL ABSTRACTION: Turbulence breaks time-reversal symmetry explicitly (energy flows large→small, not reverse). The degree of time-irreversibility is measurable: ⟨(δu)³⟩ ≠ 0 is the Kolmogorov 4/5 law — exactly the third-order asymmetry. Could the anomalous exponents ζ_p (p > 3) be constrained by the requirement that turbulence maximizes entropy production rate at each scale, subject to energy flux conservation? This would be a maximum entropy production principle applied scale-by-scale. Untested but connects ζ_p to thermodynamics.",
    },
    {
      id: "novel_missing_concept",
      label: "The Missing Concept (The Limit Analogue)",
      type: "open_problem",
      angle: "novel",
      importance: 10,
      detail: "The deepest question: is turbulence missing a mathematical concept (like calculus was missing the concept of a limit), or is it just missing computation/cleverness within existing frameworks? Evidence for missing concept: (1) K41 exhausted dimensional analysis, (2) RG fails because coupling is too strong, (3) zero modes in Kraichnan model suggest anomalous scaling = memory/information survival, but this hasn't been made precise for N-S. The candidate missing concept: a rigorous theory of how large-scale topological structure of a turbulent flow imprints on small-scale statistics — beyond what ε captures.",
    },

  ],

  edges: [
    // Core chain
    { source: "navier_stokes", target: "karman_howarth", type: "derives_from", label: "exact derivation" },
    { source: "navier_stokes", target: "closure_problem", type: "derives_from", label: "nonlinearity creates" },
    { source: "navier_stokes", target: "millennium_prize", type: "open_question", label: "existence unknown" },
    { source: "navier_stokes", target: "leray_hopf", type: "derives_from", label: "weak solutions proven" },
    { source: "navier_stokes", target: "onsager_conjecture", type: "derives_from", label: "Euler limit" },
    { source: "karman_howarth", target: "four_fifths_law", type: "derives_from", label: "inertial limit" },
    { source: "karman_howarth", target: "closure_problem", type: "derives_from", label: "hierarchy illustrated" },
    { source: "reynolds_number", target: "inertial_subrange", type: "enables", label: "Re >> 1 required" },
    { source: "reynolds_number", target: "dns", type: "blocks", label: "cost ~ Re^(9/4)" },
    { source: "richardson_cascade", target: "k41_hypotheses", type: "extends", label: "quantified by" },

    // K41 chain
    { source: "k41_hypotheses", target: "kolmogorov_scales", type: "derives_from", label: "dimensional analysis" },
    { source: "k41_hypotheses", target: "minus_five_thirds", type: "derives_from", label: "inertial range limit" },
    { source: "k41_hypotheses", target: "structure_functions", type: "derives_from", label: "predicts ζ_p = p/3" },
    { source: "k41_hypotheses", target: "inertial_subrange", type: "derives_from", label: "defines" },
    { source: "k41_hypotheses", target: "dissipation_rate", type: "assumption_of", label: "ε = global constant" },
    { source: "k41_hypotheses", target: "locality_assumption", type: "assumption_of", label: "local cascade" },
    { source: "k41_hypotheses", target: "universality_assumption", type: "assumption_of", label: "small scales universal" },
    { source: "k41_hypotheses", target: "markovian_cascade", type: "assumption_of", label: "implicit" },
    { source: "k41_hypotheses", target: "ergodicity_assumption", type: "assumption_of", label: "time=ensemble" },
    { source: "four_fifths_law", target: "k41_hypotheses", type: "critiques", label: "ζ_3 = 1 exact: K41 right here" },
    { source: "four_fifths_law", target: "anomalous_scaling", type: "critiques", label: "anchors ζ_3 = 1" },
    { source: "landau_criticism", target: "k41_hypotheses", type: "critiques", label: "ε not universal" },
    { source: "landau_criticism", target: "intermittency", type: "derives_from", label: "seed of" },
    { source: "landau_criticism", target: "k62", type: "enables", label: "motivated" },

    // K62 / intermittency
    { source: "k62", target: "k41_hypotheses", type: "extends", label: "replaces global ε" },
    { source: "k62", target: "log_normal", type: "derives_from", label: "postulates" },
    { source: "k62", target: "intermittency", type: "extends", label: "quantifies" },
    { source: "intermittency", target: "anomalous_scaling", type: "derives_from", label: "causes" },
    { source: "intermittency", target: "vortex_tubes", type: "derives_from", label: "concentrated in" },
    { source: "anomalous_scaling", target: "structure_functions", type: "derives_from", label: "ζ_p ≠ p/3" },
    { source: "log_normal", target: "she_leveque", type: "critiques", label: "log-Poisson fits better" },
    { source: "she_leveque", target: "vortex_tubes", type: "assumption_of", label: "1D filaments assumed" },
    { source: "multifractal", target: "anomalous_scaling", type: "extends", label: "parametrizes" },
    { source: "multifractal", target: "large_deviations", type: "analogous_to", label: "same structure" },
    { source: "multifractal", target: "she_leveque", type: "extends", label: "contains as special case" },

    // RG / QFT
    { source: "msr_field_theory", target: "navier_stokes", type: "derives_from", label: "reformulates stochastically" },
    { source: "msr_field_theory", target: "rg_fixed_point", type: "enables", label: "in principle" },
    { source: "msr_field_theory", target: "epsilon_expansion_failure", type: "derives_from", label: "strong coupling" },
    { source: "rg_fixed_point", target: "anomalous_scaling", type: "open_question", label: "controls exponents" },
    { source: "rg_fixed_point", target: "minus_five_thirds", type: "derives_from", label: "leading order" },
    { source: "epsilon_expansion_failure", target: "anomalous_scaling", type: "blocks", label: "cannot compute" },
    { source: "kraichnan_model", target: "zero_modes", type: "derives_from", label: "exact result" },
    { source: "zero_modes", target: "anomalous_scaling", type: "enables", label: "mechanism in solvable model" },
    { source: "zero_modes", target: "markovian_cascade", type: "critiques", label: "large-scale info survives" },
    { source: "ope_turbulence", target: "anomalous_scaling", type: "open_question", label: "would compute ζ_p" },
    { source: "conformal_2d", target: "sle_2d_turbulence", type: "derives_from", label: "CFT → SLE" },

    // Geometry / Topology
    { source: "helicity", target: "vortex_reconnection", type: "derives_from", label: "topology changes at" },
    { source: "vortex_tubes", target: "intermittency", type: "derives_from", label: "geometric carrier of" },
    { source: "vortex_reconnection", target: "helicity", type: "critiques", label: "breaks conservation" },
    { source: "vortex_reconnection", target: "markovian_cascade", type: "critiques", label: "topological events" },
    { source: "arnold_geodesic", target: "navier_stokes", type: "extends", label: "geometrizes Euler" },
    { source: "arnold_geodesic", target: "anomalous_scaling", type: "open_question", label: "negative curvature → chaos" },
    { source: "sle_2d_turbulence", target: "conformal_2d", type: "derives_from", label: "evidence for" },
    { source: "multifractal", target: "vortex_tubes", type: "extends", label: "describes geometry of" },
    { source: "onsager_conjecture", target: "intermittency", type: "derives_from", label: "h ≤ 1/3 → anomalous dissipation" },

    // Cross-domain
    { source: "phase_transition_analogy", target: "rg_fixed_point", type: "analogous_to", label: "both have fixed points" },
    { source: "phase_transition_analogy", target: "epsilon_expansion_failure", type: "critiques", label: "Wilson works; this doesn't" },
    { source: "large_deviations", target: "multifractal", type: "analogous_to", label: "same math structure" },
    { source: "non_equilibrium_thermo", target: "dissipation_rate", type: "extends", label: "ε = entropy production" },
    { source: "non_equilibrium_thermo", target: "novel_arrow_of_time", type: "enables", label: "motivates" },
    { source: "random_matrix", target: "vortex_tubes", type: "open_question", label: "eigenvalue → vortex stats?" },

    // Computational
    { source: "dns", target: "k41_hypotheses", type: "enables", label: "tests" },
    { source: "dns", target: "neural_operators", type: "enables", label: "training data" },
    { source: "neural_operators", target: "anomalous_scaling", type: "open_question", label: "can it encode ζ_p?" },
    { source: "ml_turbulence_models", target: "closure_problem", type: "enables", label: "practical workaround" },
    { source: "ml_turbulence_models", target: "closure_problem", type: "critiques", label: "not a solution" },

    // Novel abstractions
    { source: "novel_nonmarkovian", target: "markovian_cascade", type: "critiques", label: "relaxes" },
    { source: "novel_nonmarkovian", target: "zero_modes", type: "analogous_to", label: "memory = zero modes" },
    { source: "novel_nonmarkovian", target: "anomalous_scaling", type: "open_question", label: "Hurst exp = ζ_p?" },
    { source: "novel_sle_3d", target: "sle_2d_turbulence", type: "extends", label: "3D generalization" },
    { source: "novel_sle_3d", target: "vortex_tubes", type: "extends", label: "surface statistics" },
    { source: "novel_sle_3d", target: "anomalous_scaling", type: "open_question", label: "κ encodes μ?" },
    { source: "novel_info_geometry", target: "rg_fixed_point", type: "extends", label: "RG as geodesic flow" },
    { source: "novel_info_geometry", target: "anomalous_scaling", type: "open_question", label: "curvature = anomaly?" },
    { source: "novel_optimal_transport", target: "minus_five_thirds", type: "open_question", label: "optimal plan = -5/3?" },
    { source: "novel_optimal_transport", target: "dissipation_rate", type: "extends", label: "ε = transport cost" },
    { source: "novel_topological_qft", target: "helicity", type: "extends", label: "higher invariants" },
    { source: "novel_topological_qft", target: "vortex_reconnection", type: "extends", label: "topological constraints" },
    { source: "novel_arrow_of_time", target: "four_fifths_law", type: "derives_from", label: "irreversibility in ⟨(δu)³⟩" },
    { source: "novel_arrow_of_time", target: "anomalous_scaling", type: "open_question", label: "max entropy production?" },
    { source: "novel_missing_concept", target: "anomalous_scaling", type: "open_question", label: "requires new math" },
    { source: "novel_missing_concept", target: "novel_nonmarkovian", type: "extends", label: "one candidate" },
    { source: "novel_missing_concept", target: "novel_sle_3d", type: "extends", label: "one candidate" },
    { source: "novel_missing_concept", target: "zero_modes", type: "derives_from", label: "hint from Kraichnan" },

    // Assumptions ← critiqued by evidence
    { source: "forward_cascade_assumption", target: "k41_hypotheses", type: "assumption_of", label: "implicit" },
    { source: "locality_assumption", target: "k41_hypotheses", type: "assumption_of", label: "triad locality" },
  ],

  meta: {
    title: "Turbulence & Kolmogorov — Knowledge Graph",
    subtitle: "Mapping the frontier: proven results, open problems, and novel abstractions",
    generated: "2026-05-10",
    angles: ["core", "rg_qft", "geometry_topology", "cross_domain", "computational", "assumptions", "novel"],
    critical_path: [
      "navier_stokes",
      "closure_problem",
      "k41_hypotheses",
      "anomalous_scaling",
      "zero_modes",
      "novel_missing_concept"
    ],
    key_insights: [
      "The 4/5 law is the ONLY exact result from N-S for turbulence. Everything else is approximate.",
      "Zero modes in the Kraichnan model show anomalous scaling = large-scale memory surviving to small scales.",
      "K41's deepest assumption is Markovian cascade — this is what intermittency breaks.",
      "SLE(6) describes 2D turbulence; a 3D extension would be a major breakthrough.",
      "The missing concept may be: how topological structure of large-scale flow imprints on small-scale statistics."
    ]
  }
};
