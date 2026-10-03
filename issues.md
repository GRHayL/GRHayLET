## Issues Found by Agentic LLM Review

Addressed on 10-02-2026. The 21 original findings below are retained as the
pre-repair record. Current dispositions are:

| Issues | Disposition | Implemented change |
| --- | --- | --- |
| 1–3 | Repaired | Fixed neighbor mean/empty-neighbor policy; complete grid recovery seeds with rebuilt u0; mandatory conversion independent of guarded diagnostic cadence. |
| 4 | Guarded; shared correction still required | Startup rejects tabulated entropy evolution. Re-enabling it requires coordinated GRHayL packing, flux, inverse, and atmosphere changes; the dependency is outside this checkout. |
| 5 | Repaired by restriction | Hybrid entropy requires one cold segment with equal cold/thermal Gamma. |
| 6–14 | Repaired | Direct Simple thermodynamics callback; Hybrid-only Font fallback; checked face/legacy EOS errors; initialized inactive fields; raw metric derivatives; complete access metadata; synchronized perturbation-before-curl; canonical B export independent of legacy import. |
| 15 | Repaired with explicit unavailability | Internal proxy is IllinoisGRMHD::hybrid_entropy. HydroBase entropy is NaN where physical normalization/units are unsupported. |
| 16–21 | Repaired | Registered sync alias; joint outflow/Lorentz enforcement with infeasible-state error; retained reset digit; full threshold-selected denominator; current TOV output groups; corrected Lorentz identities. |

Verification uses disposable checks, with no permanent KB test infrastructure or
modified regression observations/parfiles. GCC syntax checks cover all 52 local
C units against inspected GRHayL headers and substitute CCTK declarations. UBSan
caller probes execute the four recovery bodies under controlled inversion/EOS
callbacks, consumer/diagnostic conversion at zero and sparse cadence, direct
Simple thermodynamics with actual GRHayL characteristic/HLLE routines across 31
pressures, 10,000 correlated-metric boundary cases plus an infeasible case, and
actual zero/cubic derivative stencils with floating-point traps enabled. Static
contract checks cover changed CCL access/order, carrier construction, output
requests, and comments. The canonical KB linter passes. These checks do not
establish a configured Cactus build, schedule execution, MPI/restart evolution,
production EOS inversion, or astrophysical/numerical validation. The recovery
probe deliberately sets the noncheckpointed u0 grid seed to NaN and verifies it
is reconstructed from the supplied grid velocities. Startup probes also
reject tabulated entropy, unequal-Gamma hybrid entropy, and piecewise hybrid
entropy; selected-population diagnostic probes observe 0/27 and 9/27. The
follow-up review additionally corrected legacy Prim2Con velocity inputs and
removed an unnecessary incoming validity requirement for the rebuilt proxy.


This meta-ticket collects 21 independently actionable IllinoisGRMHD findings: nine high, eight medium, and four low. Three independent review seats examined the four EOS/entropy families, recovery and conversion, numerical and magnetic routines, boundary conditions, Cactus scheduling and interfaces, relevant GRHayL dependencies, documentation, and shipped configurations. Duplicate findings have been merged and confirmed false-positive claims excluded. Conditional recovery paths, driver-dependent effects, extreme thermodynamic states, and diagnostic defects are explicitly distinguished from ordinary production failures.

<details>
<summary><strong>Issue 1 (high):</strong> Recovery neighbor averaging changes constant conservative data.</summary>

<br>

- **Description**: Recovery retries mix an unnormalized sum of N neighbors with a single central state, then divide by a count that grows on the first three attempts without adding neighbors. Neither the blend nor its final neighbor-only attempt implements the intended neighbor mean.

- **Production trigger and likelihood**: Configured primary/backup recovery fails at a positive-density point and enters the neighboring-state retries. Successful primary recovery avoids this branch; its frequency in production was not measured.

- **Impact**: A successful retry can commit artificially reduced density, energy, momentum, and active entropy/composition conservatives. B is retained separately, so magnetic-to-fluid ratios can also change.

- **Locations**: `IllinoisGRMHD/src/Hybrid/conservs_to_prims.c:129`; `IllinoisGRMHD/src/HybridEntropy/conservs_to_prims.c:134`; `IllinoisGRMHD/src/Tabulated/conservs_to_prims.c:121`; `IllinoisGRMHD/src/TabulatedEntropy/conservs_to_prims.c:124`

- **Repository evidence and checks**: Exact arithmetic with 26 identical neighbors and an identical center U gives `(29/108)U`, `(27/56)U`, `(79/116)U`, and `(26/29)U` on the four attempts. A constant-preserving blend must return U. Solver acceptance and full-run effects were not executed. Origins: `IGM-01`, `IGM2-01`.

- **Confidence**: High for the normalization defect and stated recovery-path scope.

- **Desired postcondition and fix direction**: Use a fixed neighbor mean and a convex blend with the central state. Preserve the original neighbor count and define a zero-neighbor policy.

- **Verification**: Force entry into recovery retries with constant conservatives and require every blend to preserve them. Check known nonconstant blends and truncated/empty neighborhoods separately; verify accepted recovery commits the intended repaired state.

</details>

<details>
<summary><strong>Issue 2 (high, nondefault control):</strong> Guess-off recovery receives uninitialized active primitive fields.</summary>

<br>

- **Description**: Each recovery caller initializes only B in its automatic primitive struct. With primitive guessing disabled, the dispatcher forwards this incomplete object instead of loading current grid primitives.

- **Production trigger and likelihood**: Allowed `calc_primitive_guess=no`, positive conservative density, and a selected solver that reads the incoming guess. The default guess-enabled mode can avoid this defect; exposure differs by solver.

- **Impact**: Indeterminate u0, velocity, pressure, epsilon, or temperature can cause unpredictable roots, unnecessary backups/resets, or invalid calculations. Existing restart/grid primitives are not supplied as the parameter promises.

- **Locations**: `IllinoisGRMHD/src/Hybrid/conservs_to_prims.c:66` and corresponding constructors in all four families; `GRHayL-main/implementations/GRHayLib/param.ccl:58`; `GRHayL-main/GRHayL/Con2Prim/con2prim_multi_method.c:98`; `GRHayL-main/GRHayL/Con2Prim/Hybrid/Noble/initialize_Noble.c:51`

- **Repository evidence and checks**: Caller/callee tracing establishes the missing assignments, conditional guess construction, and actual Noble reads of u0, velocities, pressure, and epsilon. Temperature-preserving table recovery also requires a meaningful seed. No guess-off evolution or full-thorn sanitizer run was executed. Origins: `IGM-02`, `IGM2-02`.

- **Confidence**: High for the incomplete-input mechanism under the stated control and solver conditions.

- **Desired postcondition and fix direction**: Load a complete physically valid current-state guess and declare all grid inputs, or reject guess-off mode until implemented. Zero-filling alone does not implement the current-state guess contract.

- **Verification**: Exercise guess-off recovery with supported Noble and relevant table methods, including restart use. Require defined seeds from the grid and agreement with the intended state; verify guess-enabled behavior separately.

</details>

<details>
<summary><strong>Issue 3 (high):</strong> Mandatory HydroBase conversion obeys diagnostic cadence and admits remainder by zero.</summary>

<br>

- **Description**: The modern converter evaluates `cctk_iteration % Convert_to_HydroBase_every` without guarding zero. Leakage RHS and retained compatibility initialization call it outside the ordinary positive-cadence scheduling guard.

- **Production trigger and likelihood**: An unconditional caller is active, the obsolete conversion thorn is inactive, and the modern cadence is its allowed/default zero. With cadence greater than one, required updates can instead be skipped on intermediate iterations. Ordinary guarded analysis calls do not establish this trigger.

- **Impact**: Zero cadence reaches undefined integer division, commonly terminating execution. Sparse cadence can leave consumer fields stale; numerical leakage effects were not independently demonstrated.

- **Locations**: `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c:7`; `IllinoisGRMHD/param.ccl:9`; `IllinoisGRMHD/schedule.ccl:197`; `IllinoisGRMHD/schedule.ccl:592`

- **Repository evidence and checks**: Parameter, branch, and schedule tracing establishes the unguarded paths. A disposable UBSan probe of the exact modulo expression exited with division-by-zero diagnostics. No complete leakage or compatibility startup was executed. The converter's access metadata is counted separately in Issue 12. Origins: `IGM-03`, `IGM2-03`.

- **Confidence**: High for remainder-by-zero reachability and cadence-controlled skipping; downstream numerical effects remain conditional.

- **Desired postcondition and fix direction**: Separate cadence-independent consumer updates from guarded diagnostic export. Guard zero in the modern exporter regardless of calling context.

- **Verification**: Exercise cadence zero and values greater than one through both unconditional callers. Required consumer fields must refresh at each intended RHS stage; diagnostic export must retain its documented cadence without division by zero.

</details>

<details>
<summary><strong>Issue 4 (high, shared mathematics):</strong> The tabulated entropy current omits baryon density.</summary>

<br>

- **Description**: Table entropy s is specific entropy per baryon, but conservative packing and flux use `sqrt(gamma)*W*s` instead of `rho_star*s`. Matching recovery inverses retain the same missing density factor.

- **Production trigger and likelihood**: Tabulated entropy evolution in smooth compression/expansion. Thermal consequences require an entropy-consuming recovery path; Palenzuela consumes entropy in this way with temperature evolution enabled. Successful energy recovery can mask the mismatch by recomputing entropy.

- **Impact**: The implemented current preserves s/rho along a parcel instead of physical specific entropy s. Entropy recovery can therefore choose incorrect temperature or pressure during otherwise adiabatic flow.

- **Locations**: `IllinoisGRMHD/src/TabulatedEntropy/prims_to_conservs.c:37`; `IllinoisGRMHD/src/TabulatedEntropy/evaluate_sources_rhs.c:23`; `GRHayL-main/GRHayL/Con2Prim/compute_conservs.c:97`; `GRHayL-main/GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.h:32`; `GRHayL-main/GRHayL/Flux_Source/tabulated_entropy/ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy.c:141`; tabulated Palenzuela/Newman entropy inverses

- **Repository evidence and checks**: The table loader supplies per-baryon entropy without density multiplication. Packing, flux, zero-source evolution, and inverses were traced. Combining physical adiabatic entropy transport with baryon continuity requires `rho*s*u^mu`; the implemented continuum equation instead changes s during compression. No production-table compression evolution was run. Origin: `IGM-04`.

- **Confidence**: High for the shared continuum/variable-contract mismatch; simulation impact is unmeasured.

- **Desired postcondition and fix direction**: Coordinate density-weighted entropy packing, flux, recovery, and atmosphere treatment across the thorn and dependency. A local packing-only change must not break the matching upstream inverse.

- **Verification**: Evolve smooth adiabatic table-EOS compression with entropy-consuming recovery and require constant material specific entropy. Check round trips, atmosphere handling, and energy-to-entropy backup transitions under the revised contract.

</details>

<details>
<summary><strong>Issue 5 (high, shared mathematics):</strong> General hybrid entropy evolves the wrong adiabatic invariant.</summary>

<br>

- **Description**: The proxy `S=P/rho^(Gamma_cold-1)` and its zero-source conservative evolution preserve `P/rho^Gamma_cold`. General hybrid thermodynamics does not preserve that quantity when thermal and cold Gammas differ and thermal pressure is positive.

- **Production trigger and likelihood**: Hybrid entropy evolution/recovery with unequal cold and thermal Gammas. A smooth single-segment flow suffices; shocks and piecewise transitions are unnecessary. Constant equal-Gamma setups form a benign subset.

- **Impact**: Entropy-based recovery can impose an incorrect thermal pressure on an adiabatic state. Successful energy recovery can conceal the inconsistent backup predictor by overwriting the proxy.

- **Locations**: `IllinoisGRMHD/src/HybridEntropy/prims_to_conservs.c:37`; `IllinoisGRMHD/src/HybridEntropy/evaluate_sources_rhs.c:24`; `GRHayL-main/GRHayL/EOS/Hybrid/NRPyEOS_compute_entropy_function.c:27`; `GRHayL-main/GRHayL/Con2Prim/Hybrid/Palenzuela1D/hybrid_Palenzuela1D_entropy.c:38`; `GRHayL-main/implementations/GRHayLib/src/initialize_and_shutdown.c:87`

- **Repository evidence and checks**: The hybrid adiabat has `P=K*rho^Gamma_cold+C*rho^Gamma_th`, consistent with the [official hybrid EOS formulas](https://einsteintoolkit.org/thornguide/EinsteinEOS/EOS_Omni/documentation.html). With K=C=1, Gamma_cold=2, Gamma_th=1.5, and rho changing from 1 to 2, the physical pressure is approximately 6.828427; the implemented invariant gives 8. Existing checking restricts one named solver rather than the general mode. No unequal-Gamma evolution was executed. Origin: `IGM-05`.

- **Confidence**: High for the smooth-flow mathematical counterexample and stated EOS range.

- **Desired postcondition and fix direction**: Restrict entropy mode to its valid constant/equal-Gamma boundary or consistently evolve and invert an appropriate thermal invariant, including an explicit cold-segment policy. This differs from the external entropy-publication defect in Issue 15.

- **Verification**: Check smooth positive-thermal-pressure compression with unequal Gammas against the analytic adiabat. Verify the valid equal-Gamma case separately and exercise segment transitions if the general mode remains supported.

</details>

<details>
<summary><strong>Issue 6 (high, extreme thermodynamic ratio):</strong> Simple EOS loses positive sound speed through cancellation.</summary>

<br>

- **Description**: Simple dispatch selects the hybrid thermodynamic helper with an auxiliary K=1 cold curve. Separate cold/thermal expressions cancel mathematically but lose tiny positive pressure when subtracting order-one terms in floating point.

- **Production trigger and likelihood**: Accepted Simple EOS states with pressure tiny relative to the placeholder cold pressure, permitted zero pressure floor, and sufficiently small velocity/B for the lost sound speed to leave zero wave bounds. The tested state is rho=1, P=1e-20, Gamma=2, zero velocity/B, and a flat metric; production incidence was not measured.

- **Impact**: Nonzero physical sound speed can become zero. Stationary unmagnetized HLLE then rejects zero total wave speed and its void wrapper aborts; nonzero flow/B can mask that terminal consequence while retaining inaccurate speeds.

- **Locations**: `GRHayL-main/GRHayL/GRHayL_Core/initialize_eos.c:40`; `GRHayL-main/GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_enthalpy_and_cs2.c:33`; `GRHayL-main/GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn0.c`; `GRHayL-main/GRHayL/Flux_Source/hybrid/ghl_calculate_HLLE_fluxes_dirn0_hybrid.c:19`; `IllinoisGRMHD/src/Hybrid/calculate_fluxes_rhs.c:189`

- **Repository evidence and checks**: Bounded actual-library probes accepted initialization and primitive enforcement, preserved positive primitive epsilon, then returned h=1, cs2=0, zero characteristic bounds, and checked HLLE status 52. The direct ideal-fluid cs2 is approximately 2e-20. The helper recomputes an internal epsilon rather than overwriting `prims.eps`. No conservative round trip, thorn stencil, or evolution was executed. Origin: `IGM-06`.

- **Confidence**: High for the finite-precision library mechanism and bounded state; full-thorn occurrence remains unmeasured.

- **Desired postcondition and fix direction**: Use direct Simple formulas for epsilon, h, and cs2 through a Simple-specific branch. Changing the placeholder K alone only moves the cancellation scale.

- **Verification**: Sweep finite positive pressure across the cancellation regime and compare with direct ideal-fluid results. Require valid nonzero stationary wave bounds and successful flux evaluation; then check conservative round trips and a supported thorn stencil.

</details>

<details>
<summary><strong>Issue 7 (high, exhausted recovery path):</strong> Simple recovery bypasses the Font1D compatibility restriction.</summary>

<br>

- **Description**: Simple shares the Hybrid recovery body, whose final fallback directly invokes Font1D regardless of EOS type or configured backup list. GRHayLib expressly rejects configured Simple/Font combinations.

- **Production trigger and likelihood**: Simple recovery exhausts configured methods and all neighbor attempts. Ordinary successful recovery avoids this path; no forced-exhaustion runtime or frequency measurement was performed.

- **Impact**: The emergency method can select an unintended cold pressure/energy state from Simple's placeholder K=1 curve and report successful Font recovery. This bypasses the validated compatibility policy; it does not prove every returned state violates the Gamma-law algebraic relation.

- **Locations**: `IllinoisGRMHD/schedule.ccl:235`; `IllinoisGRMHD/src/Hybrid/conservs_to_prims.c:180`; `IllinoisGRMHD/src/HybridEntropy/conservs_to_prims.c:189`; `GRHayL-main/implementations/GRHayLib/src/initialize_and_shutdown.c:45`; `GRHayL-main/GRHayL/Con2Prim/Hybrid/Font1D/hybrid_Font1D.c:197`

- **Repository evidence and checks**: EOS dispatch, configured-method rejection, unconditional caller bypass, and Font's cold-curve assignments were traced. Current Simple placeholder arrays are initialized, so an uninitialized-array allegation is excluded. Evidence is static. Origins: `IGM-07`, `IGM2-06`.

- **Confidence**: High for compatibility-policy bypass and unintended repair selection; runtime acceptance remains unexecuted.

- **Desired postcondition and fix direction**: Gate Font on its supported EOS or route every fallback through a validated policy. Give Simple an explicit compatible recovery/atmosphere failure policy.

- **Verification**: Force Simple recovery exhaustion and require only approved fallbacks or the stated terminal policy. Verify that valid Hybrid Font behavior remains available and diagnostics identify the actual repair method.

</details>

<details>
<summary><strong>Issue 8 (high):</strong> Tabulated face pressure-inversion errors are discarded.</summary>

<br>

- **Description**: Both tabulated families discard left/right pressure-inverse statuses. Global table bounds do not guarantee that reconstructed pressure is attainable at the particular reconstructed rho/Ye.

- **Production trigger and likelihood**: Face inversion fails, including a globally bounded pressure above the local attainable range. Incidence for realistic tables/stencils was not measured. Temperature-search failure preserves the seed T; later interpolation failure can occur after T changes.

- **Impact**: A later successful T-based callback can replace reconstructed pressure/epsilon and conceal the earlier inversion failure. The entropy family can retain entropy inconsistent with the replacement state. Final epsilon need not remain uninitialized because that callback can fill it.

- **Locations**: `IllinoisGRMHD/src/Tabulated/calculate_fluxes_rhs.c:185`; `IllinoisGRMHD/src/TabulatedEntropy/calculate_fluxes_rhs.c:187`; `GRHayL-main/GRHayL/EOS/Tabulated/interpolators/NRPyEOS_eps_and_T_from_rho_Ye_P.c`; `GRHayL-main/GRHayL/EOS/Tabulated/interpolators/NRPyEOS_eps_S_and_T_from_rho_Ye_P.c`; `GRHayL-main/GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_enthalpy_and_cs2.c:9`

- **Repository evidence and checks**: Actual-library synthetic-table probes with `P=rho*T` produced a pressure-inversion error for a locally unreachable pressure, then a successful temperature-based evaluation with a different pressure. Source tracing establishes ignored statuses and the later thermal replacement. No production EOS or complete thorn flux loop was executed. Origins: `IGM-08`, `IGM2-04`.

- **Confidence**: High for the conditional error-masking dataflow and bounded reproduction.

- **Desired postcondition and fix direction**: Handle each inverse status before characteristic/flux evaluation. Use a controlled error or documented fallback that recomputes every dependent thermal field consistently.

- **Verification**: Exercise locally unreachable pressure and other inverse errors in both families. Require explicit handling and a thermodynamically consistent fallback, including entropy; verify successful faces preserve their reconstructed pressure.

</details>

<details>
<summary><strong>Issue 9 (high, compatibility initialization):</strong> Legacy initialization ignores a checked EOS failure.</summary>

<br>

- **Description**: The retained compatibility initializer allocates an uninitialized EOS object and discards the checked initializer status. Current initialization validates a temporary candidate before committing output, so rejected parameters leave the allocated object untouched.

- **Production trigger and likelihood**: `ID_converter_ILGRMHD` is active with values allowed by legacy ranges but rejected by current GRHayL, such as Gamma=1 or invalid density bounds. Valid legacy parameters avoid this defect; ordinary modern startup is a different path.

- **Impact**: Subsequent consumers can access indeterminate EOS fields instead of receiving the intended parameter error. Normal GRHayLib initialization is skipped for this compatibility configuration and therefore does not rescue it.

- **Locations**: `IllinoisGRMHD/src/backward_compatible_initialize.c:16`; `IllinoisGRMHD/src/backward_compatible_initialize.c:42`; `IllinoisGRMHD/param.ccl:120`; `IllinoisGRMHD/schedule.ccl:565`; `GRHayL-main/GRHayL/GRHayL_Core/initialize_eos.c:188`; `GRHayL-main/implementations/GRHayLib/schedule.ccl:4`

- **Repository evidence and checks**: Allocation, local parameter ranges, scheduling exclusions, initializer early returns, and successful-output commit were traced. Checked-API probes confirmed rejected inputs leave supplied output unchanged. No obsolete-thorn startup was run or external legacy-thorn availability certified. Origins: `IGM-09`, `IGM2-05`; the higher supported severity is retained.

- **Confidence**: High for checked-error loss and uninitialized output under the stated legacy configuration.

- **Desired postcondition and fix direction**: Check allocation and initialization status and stop before any rejected output is consumed. Defensive initialization supplements rather than replaces the established error policy.

- **Verification**: Start the compatibility route with invalid-but-locally-permitted parameters and require a clear initialization error before consumer execution. Check valid legacy initialization and the mutually exclusive modern route separately.

</details>

<details>
<summary><strong>Issue 10 (medium):</strong> Generic helpers read uninitialized inactive struct fields.</summary>

<br>

- **Description**: Mode-specific primitive/conservative constructors omit inactive entropy or Ye fields, but generic undensitization and packing helpers unconditionally read both. Primitive enforcement does not universally populate these inactive members.

- **Production trigger and likelihood**: Ordinary non-entropy operation or hybrid operation without Ye. This does not require guess-off mode and is distinct from Issue 2's solver-active inputs.

- **Impact**: Indeterminate floating-point arithmetic violates a defined caller/helper boundary and can surface under sanitizers or traps. Unused outputs are normally discarded; default active mass/energy corruption was not established.

- **Locations**: `IllinoisGRMHD/src/Hybrid/prims_to_conservs.c:28`; corresponding constructors in other families and recovery/boundary paths; `GRHayL-main/GRHayL/Con2Prim/undensitize_conservatives.c:41`; `GRHayL-main/GRHayL/Con2Prim/compute_conservs.c:97`; `GRHayL-main/GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c`

- **Repository evidence and checks**: Assignment/read tracing confirms omitted members and unconditional scalar expressions. Particular optimizations or ignored outputs do not establish an initialization contract. No complete-thorn sanitizer or uninitialized-value trapping run was executed. Origins: `IGM-10`, `IGM2-07`.

- **Confidence**: High for the incomplete-struct input mechanism; specific active-output damage is unproved.

- **Desired postcondition and fix direction**: Fully initialize each quantity struct with deterministic inactive placeholders before setting active physical fields. Preserve a physically meaningful guess where required by Issue 2.

- **Verification**: Exercise conversion, recovery, and boundary constructors in all four modes under a suitable uninitialized-value checker. Require defined helper inputs and unchanged intended active outputs.

</details>

<details>
<summary><strong>Issue 11 (medium):</strong> Derivative tensors are initialized as physical metrics.</summary>

<br>

- **Description**: Metric/lapse derivatives are passed to `ghl_initialize_metric`, which computes reciprocal lapse and a spatial-matrix inverse. A derivative tensor need not be an invertible physical metric.

- **Production trigger and likelihood**: Valid zero or singular derivatives, including ordinary flat constant geometry. Nontrapping operation can mask the misuse because current source kernels consume only raw derivative members.

- **Impact**: Unnecessary inverse calculations produce nonfinite unused auxiliaries and floating-point exception flags. Enabled traps can stop a valid calculation; default NaN RHS is not established.

- **Locations**: `IllinoisGRMHD/src/compute_metric_derivs.c:41`; all four source-RHS callers; `GRHayL-main/GRHayL/GRHayL_Core/initialize_metric.c:60`; `GRHayL-main/GRHayL/Flux_Source/ghl_calculate_source_terms.c:37`

- **Repository evidence and checks**: Actual-library probes of zero derivatives recorded `FE_DIVBYZERO=1`, `FE_INVALID=1`, infinite reciprocal lapse, and NaN inverse entries. Source tracing confirms the kernel reads raw derivative lapse/shift/covariant metric, not those auxiliaries. No trapping Cactus evolution was run. Origins: `IGM-11`, `IGM2-08`.

- **Confidence**: High for invalid auxiliary arithmetic and the conditional trap hazard.

- **Desired postcondition and fix direction**: Populate a derivative carrier directly with raw lapse, shift, and symmetric covariant entries; use a derivative-specific type/helper if appropriate. Avoid physical-metric inversion of derivative data.

- **Verification**: Evaluate flat and singular derivative stencils with floating-point traps enabled and require no invalid inverse work. Compare nonconstant source terms with independently calculated derivatives to preserve the intended RHS.

</details>

<details>
<summary><strong>Issue 12 (medium, driver-dependent interface):</strong> Schedule access declarations omit accumulator inputs and converter fields.</summary>

<br>

- **Description**: Tmunu and A RHS routines use `+=` while their schedule occurrences omit accumulated arrays from READS. Leakage's shared converter also accesses metric/B and Lorentz-factor/Bvec fields beyond its declared velocity-only footprint.

- **Production trigger and likelihood**: Execution under a driver/checking/transfer mode that relies on these declarations. Default PreSync-off behavior can mask the mismatch; no affected driver execution was demonstrated.

- **Impact**: Validity, synchronization, poisoning, or transfer machinery can omit required old values or fail to track actual updates. The additive numerical operations themselves are intentional.

- **Locations**: `IllinoisGRMHD/schedule.ccl:169`; `IllinoisGRMHD/src/compute_Tmunu.c:42`; `IllinoisGRMHD/schedule.ccl:208`; `IllinoisGRMHD/src/evaluate_phitilde_and_A_gauge_rhs.c:128`; `IllinoisGRMHD/schedule.ccl:197`; `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c:69`

- **Repository evidence and checks**: Source-to-declaration comparison establishes concrete missing inputs/outputs for each occurrence. The [official Cactus/ET guide](https://einsteintoolkit.org/usersguide/UsersGuide.html) describes the declarations' validity and synchronization roles. No strict PreSync, poisoning, or accelerator run was executed. Origins: `IGM-12`, `IGM2-13`.

- **Confidence**: High for access-contract mismatch; downstream driver failures remain conditional.

- **Desired postcondition and fix direction**: Declare actual accumulator reads, outputs, and regions, or introduce a genuinely velocity-only consumer converter. Another occurrence's complete declaration does not repair this one.

- **Verification**: Inspect each corrected schedule item under the intended validity checker and communication/transfer mode. Require existing accumulator contributions to survive and actual converter inputs/outputs to be tracked.

</details>

<details>
<summary><strong>Issue 13 (medium, optional initialization path):</strong> Initial potential perturbation lacks ordering before magnetic reconstruction.</summary>

<br>

- **Description**: Initial perturbation changes Ai, but perturbation and the A-to-B curl both follow import and precede Prim2Con without a relative dependency. Modified ghost potentials also have no new explicit synchronization after import.

- **Production trigger and likelihood**: `perturb_initial_data=yes` with nonzero amplitude. Default perturbation-off evolution avoids the path; favorable realized ordering can also avoid the stale B state.

- **Impact**: The declared graph permits curl-before-perturbation, leaving initial B/conservatives tied to older A while evolution uses perturbed A. Ghost disagreement additionally depends on ownership/communication requirements. No realized bad schedule or MPI constraint error was measured.

- **Locations**: `IllinoisGRMHD/schedule.ccl:63`; `IllinoisGRMHD/src/Hybrid/perturb_primitives.c:23` and corresponding family functions; `IllinoisGRMHD/src/compute_B_and_Bstagger_from_A.c:50`

- **Repository evidence and checks**: Dependency and write tracing establishes the missing declared perturbation-before-curl guarantee. The later conservative-only perturbation option does not alter Ai and is not included in this claim. RNG reproducibility is not counted as another proven defect. Origins: `IGM-13`, `IGM2-11`.

- **Confidence**: High for the missing ordering guarantee; scheduled/MPI consequences remain conditional.

- **Desired postcondition and fix direction**: Order initial potential perturbation before curl and synchronize modified owner values when required. Derive B and initial conservatives from the final initialized A.

- **Verification**: Inspect the generated initialization schedule and run nonzero perturbations across relevant patches/processes. Compare stored B with a fresh curl of final A and check required owner/ghost consistency.

</details>

<details>
<summary><strong>Issue 14 (medium, field convention):</strong> Legacy magnetic import also forces noncanonical HydroBase B output.</summary>

<br>

- **Description**: One `rescale_magnetics` switch divides legacy Avec by sqrt(4pi) on import and multiplies already normalized internal B by sqrt(4pi) on HydroBase export. It cannot select legacy import together with canonical export.

- **Production trigger and likelihood**: `rescale_magnetics=yes`, enabled export, and a consumer expecting canonical HydroBase Bvec. A legacy diagnostic consumer may intentionally expect the older output convention; cases using the switch off do not exercise this combination.

- **Impact**: Canonical consumers receive B too large by sqrt(4pi), with quadratic diagnostics potentially too large by 4pi. Consumer execution and numerical impact were not measured.

- **Locations**: `IllinoisGRMHD/src/convert_HydroBase_to_IllinoisGRMHD.c:7`; `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c:19`; `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c:96`; `IllinoisGRMHD/param.ccl:18`; `IllinoisGRMHD/doc/documentation.tex:155`

- **Repository evidence and checks**: Scaling/dataflow was compared with the [official normalized HydroBase magnetic-field definition](https://einsteintoolkit.org/thornguide/EinsteinBase/HydroBase/documentation.html) and current library convention. The multiplicative discrepancy follows exactly from the export factor. No downstream consumer was executed. Origin: `IGM2-09`.

- **Confidence**: High for the stated canonical-consumer contract mismatch.

- **Desired postcondition and fix direction**: Separate legacy input conversion from explicitly optional legacy diagnostic output. Publish canonical HydroBase Bvec independently of the import compatibility setting.

- **Verification**: Import a known legacy potential and export to a canonical consumer, checking B normalization and a quadratic diagnostic. Test any retained legacy-output mode explicitly and document both directions.

</details>

<details>
<summary><strong>Issue 15 (medium, field convention):</strong> Hybrid/Simple entropy proxy is published as physical specific entropy.</summary>

<br>

- **Description**: The library proxy `S=P/rho^(Gamma-1)` is assigned to HydroBase's specific-entropy-per-particle field, although it is a different quantity. Internal proxy transport validity does not establish that external field contract.

- **Production trigger and likelihood**: Hybrid/Simple entropy mode with a HydroBase entropy consumer. The mismatch exists even in a constant equal-Gamma ideal-fluid case where Issue 5's internal adiabatic defect is absent.

- **Impact**: A consumer can interpret adiabatic compression as a change in physical entropy or use values with the wrong physical meaning. No coupled-consumer failure was executed.

- **Locations**: `IllinoisGRMHD/src/HybridEntropy/prims_to_conservs.c:54`; `IllinoisGRMHD/src/HybridEntropy/conservs_to_prims.c:241`; `IllinoisGRMHD/src/HybridEntropy/hydro_outer_boundaries.c:224`; `GRHayL-main/GRHayL/EOS/Hybrid/NRPyEOS_compute_entropy_function.c:27`

- **Repository evidence and checks**: Publication sites and proxy definition were compared with the [HydroBase entropy definition](https://einsteintoolkit.org/thornguide/EinsteinBase/HydroBase/documentation.html). For Gamma=2 and the ideal-fluid isentrope `P=rho^2`, the proxy scales with rho while physical specific entropy stays constant. Evidence is source/analytic; no external consumer was run. Origin: `IGM2-10`.

- **Confidence**: High for the external variable-contract mismatch and analytic counterexample.

- **Desired postcondition and fix direction**: Own the proxy in a thorn variable and publish a correctly defined physical entropy, or explicitly make that diagnostic unavailable where unsupported. Preserve the matching internal recovery contract.

- **Verification**: Compare published entropy along a known ideal-fluid isentrope and verify constancy and documented units. Check that internal entropy recovery still round-trips its intended proxy independently.

</details>

<details>
<summary><strong>Issue 16 (medium, schedule interface):</strong> Recovery synchronization depends on the function name instead of its AS alias.</summary>

<br>

- **Description**: The synchronization item is scheduled as `IllinoisGRMHD_sync_conservatives`, but A boundaries specify `after IllinoisGRMHD_sync`. The latter is the source function name rather than the registered scheduling name.

- **Production trigger and likelihood**: Recovery traversal relying on this explicit sync-before-boundary/curl dependency. Incidental order or other driver/MoL/PreSync communication may supply the needed synchronization; no late-sync MPI execution was demonstrated.

- **Impact**: The declared consumer chain lacks its intended dependency on conservative/potential synchronization. Missing communication can affect consumers if no other mechanism supplies it.

- **Locations**: `IllinoisGRMHD/schedule.ccl:124`; `IllinoisGRMHD/schedule.ccl:130`; official Cactus `lib/sbin/CreateScheduleBindings.pl`; official Cactus `src/schedule/ScheduleCreater.c`

- **Repository evidence and checks**: Official [binding generation](https://api.bitbucket.org/2.0/repositories/cactuscode/cactus/src/master/lib/sbin/CreateScheduleBindings.pl) registers AS and preserves dependency strings. Official [schedule creation](https://api.bitbucket.org/2.0/repositories/cactuscode/cactus/src/master/src/schedule/ScheduleCreater.c) matches literal item names and adds no edge for an absent target. This resolves the first audit's open alias question at the inspected implementation level; no runtime schedule was executed. Origin: `IGM2-12`.

- **Confidence**: High for the missing edge under the inspected standard Cactus implementation; runtime effects remain conditional.

- **Desired postcondition and fix direction**: Reference `IllinoisGRMHD_sync_conservatives` consistently, or remove the alias consistently, so the intended ordering edge exists.

- **Verification**: Inspect the generated/resolved recovery schedule and exercise communication across patches/processes. Require synchronization to precede the boundary/curl/C2P consumer chain.

</details>

<details>
<summary><strong>Issue 17 (medium, shifted physical boundary):</strong> Speed limiting can restore inflow after an outflow clamp.</summary>

<br>

- **Description**: The boundary first clips inward coordinate normal velocity to zero, then the general limiter rescales `v+beta` and subtracts beta. With correction c below one, the new clipped component becomes `beta_normal*(c-1)` and can be inward again.

- **Production trigger and likelihood**: Outflow mode at an active coarsest-level physical boundary after initialization, nonzero normal shift, and a state whose clamp activates speed limiting. The incoming interior state need not initially violate the Lorentz bound; zero normal shift avoids this example.

- **Impact**: Final ghost velocity can violate the requested coordinate outflow sign. Net HLLE boundary flux or evolution error was not measured; the claim uses coordinate transport velocity rather than confusing it with Eulerian velocity.

- **Locations**: `IllinoisGRMHD/src/Hybrid/hydro_outer_boundaries.c:54`; enforcement/storage at line 203 and corresponding sequence in all four families; `GRHayL-main/GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:21`

- **Repository evidence and checks**: An actual-library probe used a flat spatial metric, alpha=1, beta=(0.5,0,0), Wmax=10, and admissible interior v=(-0.3,0.9,0), with W approximately 2.582. After the thorn-style vx clamp, the limiter returned vx approximately -0.0167914 and W=10. No full boundary-flux run was performed. Origin: `IGM2-17`.

- **Confidence**: High for the bounded ghost-state sign violation and caller sequence.

- **Desired postcondition and fix direction**: Enforce the normal coordinate-sign constraint and Lorentz bound jointly, then recompute u0/conservatives. Simple re-clipping can violate the speed bound again; define handling for infeasible constraints.

- **Verification**: Exercise upper/lower faces, both shift signs, and tangential velocities that activate limiting. Require final finite states to satisfy both documented constraints and check the resulting boundary flux.

</details>

<details>
<summary><strong>Issue 18 (low, diagnostic):</strong> The terminal atmosphere-reset digit is overwritten.</summary>

<br>

- **Description**: Total recovery exhaustion adds 100 directly to the grid `failure_checker`, then a later assignment overwrites the field from a local accumulator that never received the increment.

- **Production trigger and likelihood**: All available recovery attempts fail and the cell is reset to atmosphere. Frequency was not measured; successful recovery does not take this terminal branch.

- **Impact**: The spatial diagnostic loses the promised reset digit, obscuring where recovery discarded a state. Aggregate failure counts and the atmosphere reset itself still occur.

- **Locations**: `IllinoisGRMHD/src/Hybrid/conservs_to_prims.c:202` and line 239; `IllinoisGRMHD/src/HybridEntropy/conservs_to_prims.c:212` and line 251; `IllinoisGRMHD/src/Tabulated/conservs_to_prims.c:180` and line 219; `IllinoisGRMHD/src/TabulatedEntropy/conservs_to_prims.c:186` and line 226

- **Repository evidence and checks**: Deterministic assignment tracing shows the local accumulator lacks the 100 and the grid field is overwritten. Disposable control-flow/arithmetic checks corroborate the overwrite; no injected full recovery failure was run. Origins: `IGM-15`, `IGM2-14`.

- **Confidence**: High for the diagnostic overwrite under terminal exhaustion.

- **Desired postcondition and fix direction**: Add the terminal digit to the local accumulator and publish the grid diagnostic once. Keep aggregate counts and family-specific decoder wording consistent.

- **Verification**: Force terminal reset in all four families and require the spatial digit and aggregate count to identify the same event. Confirm successful/backup-only recovery retains its distinct markers.

</details>

<details>
<summary><strong>Issue 19 (low, diagnostic):</strong> The threshold-selected failure denominator counts only failed points.</summary>

<br>

- **Description**: `pointcount_inhoriz` increments only alongside `failures_inhoriz` in the total-failure branch when `sqrt_detgamma > psi6threshold`. Successful points passing the same threshold never enter the denominator.

- **Production trigger and likelihood**: Verbose recovery diagnostics for a threshold-selected region. The defect is in reported population accounting rather than evolution; no runtime frequency/log was collected.

- **Impact**: The displayed failure/population statistic degenerates to F/F instead of failures over all selected points. The threshold itself is a proxy, not independent apparent-horizon detection.

- **Locations**: `IllinoisGRMHD/src/Hybrid/conservs_to_prims.c:205` and line 344; `IllinoisGRMHD/src/HybridEntropy/conservs_to_prims.c:217`; `IllinoisGRMHD/src/Tabulated/conservs_to_prims.c:185`; `IllinoisGRMHD/src/TabulatedEntropy/conservs_to_prims.c:191`

- **Repository evidence and checks**: All counter assignments were traced; denominator and numerator increment together only on terminal failures. Evidence is static and no log reproduction is claimed. Origin: `IGM2-15`.

- **Confidence**: High for the counter-accounting defect and diagnostic-only scope.

- **Desired postcondition and fix direction**: Count every threshold-selected point independently of recovery success, and label the selected region accurately.

- **Verification**: Use a selected population containing successes and failures and require F/N, including zero-failure and all-failure cases. Confirm the label describes the threshold criterion without implying independent horizon detection.

</details>

<details>
<summary><strong>Issue 20 (low, example configuration):</strong> The modern magnetized-TOV example requests unstored legacy density.</summary>

<br>

- **Description**: The modern example omits `ID_converter_ILGRMHD` but requests `IllinoisGRMHD::rho_b` and `grmhd_primitives_allbutBi` in output lists. Storage and population of those fields belong to the omitted legacy path.

- **Production trigger and likelihood**: Running the shipped `par/magnetizedTOV.par` configuration with its listed output requests. The separate regression parfile uses HydroBase density and is not accused of this stale request.

- **Impact**: Requested legacy fields cannot provide the intended current density output. The exact output-thorn warning, skip, or failure policy was not exercised.

- **Locations**: `IllinoisGRMHD/par/magnetizedTOV.par:10`; output requests at lines 197, 203, and 217; `IllinoisGRMHD/interface.ccl:126`; `IllinoisGRMHD/schedule.ccl:565`; `IllinoisGRMHD/src/backward_compatible_data.c:13`

- **Repository evidence and checks**: Thorn activation, conditional storage, compatibility writer, and output strings were compared directly. No output-thorn run was performed. Origin: `IGM2-18`.

- **Confidence**: High for the example/storage mismatch; output-driver consequences remain unexecuted.

- **Desired postcondition and fix direction**: Request current HydroBase primitive fields and appropriate current Illinois velocity groups so modern output lists reference populated storage.

- **Verification**: Run the modern example's output setup and require density/primitive requests to resolve to stored current values. Compare with the existing regression configuration and confirm no obsolete converter is needed for these outputs.

</details>

<details>
<summary><strong>Issue 21 (low, mathematical documentation):</strong> Velocity-conversion comments invert the Lorentz-factor identity.</summary>

<br>

- **Description**: Comments equate `alpha*u^0` with `1/sqrt(1+gamma^{ij}u_i*u_j)`. Four-velocity normalization instead gives `W=alpha*u^0=sqrt(1+gamma^{ij}u_i*u_j)`; the printed reciprocal is 1/W.

- **Production trigger and likelihood**: A reader uses the duplicated comments when maintaining velocity conversions or diagnostics. The executable conversion and W calculation are correct; this is not a production computed-W defect.

- **Impact**: Misleading mathematical documentation can propagate into later implementation or numerical debugging.

- **Locations**: `IllinoisGRMHD/src/convert_HydroBase_to_IllinoisGRMHD.c:41`; `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c:40`; `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c:61`

- **Repository evidence and checks**: Analytic normalization, duplicated comment strings, and actual executable formulas were compared. No runtime calculation error is alleged. Origin: `IGM2-19`.

- **Confidence**: High for the documentation contradiction and bounded maintenance impact.

- **Desired postcondition and fix direction**: Correct all duplicated identities while preserving the coordinate/Eulerian velocity distinction and current executable formulas.

- **Verification**: Compare every comment occurrence with the normalization identity. Use a finite nonzero covariant spatial four-velocity to distinguish W from 1/W and confirm the documented conversion agrees with the existing calculation.

</details>

### Review boundaries and excluded findings

This IllinoisGRMHD list synthesizes [IGM_issues_final.md](/work/IGM_issues_final.md) and [IGM2_issues_final.md](/work/IGM2_issues_final.md), incorporating [IGM_fp.md](/work/IGM_fp.md). Their 34 actionable entries reduce to 21 unique findings after merging eleven duplicate pairs and excluding two false-positive claims. Issue numbers preserve the former SYN ordering, and each evidence field records original identifiers. The original audit reports are historical inputs. The issue descriptions below retain their pre-repair evidence; the disposition table above records the current source repairs and validation limits.

The audits covered all four EOS/entropy families and relevant common numerical/magnetic code, CCL/build interfaces, current local GRHayL dependencies, documentation, examples, and test declarations. Executed evidence includes exact retry arithmetic, scalar modulo-zero UBSan reproduction, synthetic-table inverse/enthalpy calls, Simple thermodynamic/characteristic/HLLE calls, checked initializer boundaries, derivative-metric floating-point exceptions, and shifted-outflow limiting. These bounded checks establish the stated mechanisms rather than a full IllinoisGRMHD/Cactus build or regression pass. No MPI evolution, production-table simulation, convergence study, full-thorn sanitizer run, or astrophysical validation is claimed. Verification bullets above specify proposed follow-up checks; they were not newly executed by this formatting update and do not add permanent tests or CI changes.
