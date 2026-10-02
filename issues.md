## Issues Found by Agentic LLM Review

This meta-ticket collects **13 GRHayLHD findings: four high, six medium, and three low severity**. Reviewer examined the thorn source, current supplied GRHayL/GRHayLib interfaces and dependencies, Einstein Toolkit declarations, documentation, parameter files, tests, and historical observations. Three fresh seats then reviewed the findings for false positives. The list incorporates that analysis: the atmosphere-reset diagnostic is Low severity, and the Simple-EOS Font1D finding concerns a recovery-policy bypass and thermal reset, without establishing algebraically invalid final primitives.

Severity applies under the stated trigger. Source/API defects, conditional framework or precision effects, and diagnostic/documentation defects are distinguished from demonstrated production failures. Original audit IDs and historical line numbers are retained in Locations. The findings below describe the pre-fix source; current dispositions are recorded next.


## Resolution status (10-02-2026)

All 13 reported mechanisms are valid and have source fixes. The dispositions
below concern the current source and focused checks; they do not certify a
complete Cactus build, scheduled execution, numerical validity, or a real-table
simulation. No configured Cactus environment or EOS table was supplied.

| Issue | Implemented disposition | Verification and remaining limits |
| --- | --- | --- |
| 1 | Required converter runs unconditionally after hydro flux RHS on every `MoL_CalcRHS` invocation. A separate initial/analysis diagnostic wrapper checks zero before modulo. | Actual converter probe passes intervals 0, 1, 2, 3 and diagnostic cadence checks under UBSan. Local leakage consumer is declared in `MoL_PostRHS`; full substage scheduling remains unrun. |
| 2 | Each active conservative uses `w*(neighbor_sum/N) + (1-w)*center` with fixed count; an empty neighborhood skips averaging. | Forced-recovery probes pass all retries in four families for constant and nonconstant states, clipped edges, and one-point grids. Solvers are stubbed; simulation convergence is unmeasured. |
| 3 | Each recovery entry rejects `ghl_params->calc_prim_guess == false` before solver calls. Previous-state guessing remains unsupported and is documented. | Forced enabled/disabled tests pass for every family. Guard reads the live library configuration at every entry; real restart/configuration recovery remains unrun. |
| 4 | Explicit Font1D emergency fallback is limited to Hybrid EOS. Simple exhausts configured methods/retries and resets to its configured atmosphere, followed by limiting and reprojection. | Forced-failure probes verify Simple never calls Font1D in either entropy mode and Hybrid retains the fallback. Real EOS numerical behavior remains unmeasured. |
| 5 | Both tabulated face inversions send their return status to `ghl_abort_if_error` before characteristic speeds or fluxes. | Successful inversions and injected left/right errors pass in both modes; errors stop before flux evaluation. No real-table inversion was run. |
| 6 | Complete primitive and conservative API objects are initialized before active fields are loaded, including retry, atmosphere/reprojection, and boundary sites. | Header checks and forced recovery input checks pass, including inactive entropy/Ye. Complete instrumented Cactus boundary/evolution execution remains unrun. |
| 7 | Both nonentropy Hybrid speed-limit accumulators start at false. | Actual flux units pass an accumulation/clipping stub probe under UBSan; the original Hybrid flux path reproduced an invalid boolean load. Source flag lifetime is also inspected. |
| 8 | Derivative storage is initialized and raw lapse/shift/symmetric metric components are packed directly without physical inversion. | Actual helper passes constant geometry and independent polynomial derivatives in all directions with `FE_DIVBYZERO`/`FE_INVALID` traps enabled. Full source-term evolution remains unrun. |
| 9 | Tmunu declares reads of its additive targets; leakage conversion declares metric reads and Lorentz-factor writes. | Source/CCL field sets reconciled. Deployed-driver validity/poisoning checks remain unavailable. |
| 10 | API-facing flux stencils, reconstruction outputs, wave speeds, and callback signatures use `double`; grid pointers remain `CCTK_REAL`. | All 33 inspected source units pass syntax checks with double, float, and long-double Cactus shims against upstream headers. This does not establish complete-stack precision support or generated bindings/linking. |
| 11 | Final reset marker accumulates locally and survives the single grid assignment. Regional denominator counts all points above `psi6threshold`; log label is `AbovePsi6Threshold`. | Four-family forced-reset/mixed-success probes pass marker and independently calculated regional counts. Aggregate reset and reprojection paths remain intact. |
| 12 | A conservative sync group always synchronizes core state and conditionally synchronizes Ye/entropy using the same predicates as storage. | All six Simple/Hybrid/Tabulated entropy combinations were statically inspected. Generated schedules and deployed-driver communication remain unrun. |
| 13 | Lorentz identities, HydroBase copy guidance, EOS/guess restrictions, and hydro boundary banners are corrected. Affected KB claims and overlapping KB issues are reconciled. | Normalization/assignments and complete guide context inspected; canonical KB lint passes. |

Validation used disposable probes outside the repository, actual thorn source
units, downloaded upstream GRHayL headers and selected helpers, and independent
Cactus argument shims. Recovery/EOS/flux stubs force exceptional paths and check
caller contracts; they do not validate external solvers, tables, or numerical
flux accuracy. Existing Cactus tests and stored observations were not changed.

<details>
<summary><strong>Issue 1 (high):</strong> Leakage conversion can evaluate an integer remainder with a zero divisor.</summary>

<br>

- **Description**: The converter begins with `cctk_iteration % Convert_to_HydroBase_every`. Zero is the legal default for disabling diagnostic copying, but a separate RHS instance runs whenever NRPyLeakageET is active, independently of that interval. The guards on initial-data and analysis instances do not protect this RHS instance.

- **Production trigger and likelihood**: Leakage enabled with the default conversion interval zero. The parameter and schedule establish a supported configuration reaching the expression; complete leakage-enabled execution was not run. With intervals greater than one, the same gate also skips RHS refreshes on most iterations.

- **Impact**: Undefined integer remainder and a possible runtime failure. Required physical refresh shares a diagnostic cadence gate, but no particular leakage-consumer numerical error or required cadence was established.

- **Locations**: `GRHayLHD/src/convert_GRHayLHD_to_HydroBase.c:8`; `GRHayLHD/param.ccl:9`; `GRHayLHD/schedule.ccl:137`. Original audit ID: `F01`.

- **Repository evidence and checks**: An isolated C99/O0 probe included the actual converter with stubbed Cactus arguments and `-fsanitize=undefined -fno-sanitize-recover=all`. It exited nonzero with a division-by-zero report at converter line 8. This was a function-level reproduction, without an external leakage consumer or Cactus scheduler.

- **Confidence**: High for the zero-divisor mechanism and declared trigger; consumer-specific cadence consequences remain unverified.

- **Desired postcondition and fix direction**: Required physical conversion must not depend on an optional diagnostic interval. Separate the two paths, define the required RHS/substage refresh policy, and handle zero safely. A zero guard alone does not settle the physical cadence question.

- **Verification**: Exercise intervals zero, one, and greater than one with leakage active in the actual schedule, including MoL substages. Confirm defined arithmetic and refresh timing before each identified consumer.

</details>

<details>
<summary><strong>Issue 2 (high):</strong> Recovery neighbor blending fails to preserve a constant conservative state.</summary>

<br>

- **Description**: All four recovery families compute `(w*S + (1-w)*C)/n_avg`, where `S` is a neighbor sum and `C` is the central value. The divisor also grows across the first three retries. This incorrectly divides the central contribution by a neighbor count and changes the normalization between attempts.

- **Production trigger and likelihood**: Initial primitive recovery fails and enters neighbor averaging in any EOS/entropy family. The failure path is explicit; its frequency in actual simulations was not measured. No admissible uniform state was shown to naturally cause the initial solver failure.

- **Impact**: Recovery candidates do not preserve constant data. If a candidate succeeds, subsequent conservative reprojection can carry the erroneous reduction into evolved density, energy, momentum, and active scalar fields.

- **Locations**: `GRHayLHD/src/Hybrid/conservs_to_prims.c:147`; `GRHayLHD/src/HybridEntropy/conservs_to_prims.c:153`; `GRHayLHD/src/Tabulated/conservs_to_prims.c:140`; `GRHayLHD/src/TabulatedEntropy/conservs_to_prims.c:144`. Original audit ID: `F02`.

- **Repository evidence and checks**: An arithmetic probe using the source weights and cumulative divisors gave approximately 0.268519, 0.482143, 0.681034, and 0.896552 for a center and 26 neighbors all equal to one. A constant-preserving convex blend would give one on every retry. No solver-failure simulation was performed by this probe.

- **Confidence**: High for the normalization error and bounded failure-path scope.

- **Desired postcondition and fix direction**: Normalize the neighbor sum once using the actual fixed count, then use `w*(S/N) + (1-w)*C`. Apply the same rule to all active entropy/Ye terms and define edge-neighborhood handling.

- **Verification**: Force entry into each retry in all four families. Check constant preservation and independently calculated blends, including neighborhoods clipped by grid edges.

</details>

<details>
<summary><strong>Issue 3 (high):</strong> Disabling primitive guessing supplies uninitialized active inputs to guess-dependent solvers.</summary>

<br>

- **Description**: Each recovery family creates an automatic primitive structure and initializes only its magnetic components. Loading the previous primitive state is unimplemented. When `calc_primitive_guess=no`, current upstream dispatch skips guess construction and passes this structure to methods that read primitive guesses.

- **Production trigger and likelihood**: The permitted, recovery-steerable option is disabled and a selected recovery method consumes the supplied guess. Default enabled guessing does not trigger this specific omission. The configuration is allowed, but complete runtime failure and deployment frequency were not measured.

- **Impact**: Active velocity, energy, pressure, temperature, or related solver inputs can be indeterminate, permitting invalid reads, nondeterministic convergence, or incorrect recovery behavior.

- **Locations**: `GRHayLHD/src/Hybrid/conservs_to_prims.c:66` and equivalent primitive setup in the other three families; `GRHayL-main/implementations/GRHayLib/param.ccl:58`; `GRHayL-main/GRHayL/Con2Prim/con2prim_multi_method.c:98`. Original audit ID: `F03`.

- **Repository evidence and checks**: Source tracing establishes the caller omission, conditional upstream guess construction, and active reads in guess-dependent Hybrid/Tabulated Noble paths. This is more than copying unspecified padding or discarding an inactive output. No complete Cactus recovery run with the option disabled was performed.

- **Confidence**: High for the input-contract violation under the stated option/method combination.

- **Desired postcondition and fix direction**: Supply a complete valid previous primitive state and declare its scheduled reads, or reject disabled guessing until supported. Merely zeroing the structure does not implement a previous-state guess.

- **Verification**: Capture primary and backup solver inputs with guessing enabled and disabled in each EOS/entropy mode. Include supported recovery-time steering and restart paths; require valid documented guesses or an explicit rejection.

</details>

<details>
<summary><strong>Issue 4 (high, conditional recovery path):</strong> Simple EOS reaches Font1D despite the adapter's explicit restriction.</summary>

<br>

- **Description**: Simple and Hybrid select the same thorn recovery families. After configured methods and neighbor retries fail, both Hybrid families directly call Font1D regardless of EOS selection. Current GRHayLib explicitly rejects configured Font1D for Simple EOS, so the thorn's extra call bypasses that recovery-policy restriction.

- **Production trigger and likelihood**: Simple EOS, with or without entropy, exhausts the preceding recovery attempts and reaches the direct fallback. The source route exists; actual occurrence and numerical effect in a Simple simulation were not measured.

- **Impact**: A successful fallback selects a cold thermal state from Simple's artificial `K=1` representation instead of recovering the independent thermal state from conserved energy, before conservative reprojection. This can materially reset thermodynamics. It does not establish algebraically invalid final ideal-gas primitives: the subsequent Simple limiter recomputes epsilon consistently from limited pressure and density.

- **Locations**: `GRHayLHD/schedule.ccl:165`; `GRHayLHD/src/Hybrid/conservs_to_prims.c:189`; `GRHayLHD/src/HybridEntropy/conservs_to_prims.c:198`; `GRHayL-main/implementations/GRHayLib/src/initialize_and_shutdown.c:213`; `GRHayL-main/GRHayL/Con2Prim/Hybrid/Font1D/hybrid_Font1D.c:197`; `GRHayL-main/GRHayL/GRHayL_Core/initialize_eos.c:123`; `GRHayL-main/GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c:83`. Original audit ID: `F04`.

- **Repository evidence and checks**: Source tracing confirms the explicit restriction and the bypassing call. Simple supplies one piece with `K=1`, the ideal-fluid Gamma, and zero integration constant. The cold helper bounds density and supplies cold pressure/epsilon; the later limiter enforces `eps=P/[rho*(Gamma-1)]`. This counterevidence narrows the claim to policy and thermal-state behavior. No forced Simple fallback execution was performed.

- **Confidence**: High for the restriction bypass and traced thermal prescription; runtime occurrence and magnitude remain unverified.

- **Desired postcondition and fix direction**: Select an explicit EOS-compatible emergency recovery/reset policy. Restrict this cold fallback to compatible Hybrid EOS unless Simple support and its thermal-reset semantics are deliberately defined and validated.

- **Verification**: Force preceding methods to fail in both Simple entropy modes. Inspect the chosen fallback, pressure/energy relation, retained or reset thermal state, and conservative reprojection against the documented policy.

</details>

<details>
<summary><strong>Issue 5 (medium):</strong> Tabulated face reconstruction discards pressure-inversion errors.</summary>

<br>

- **Description**: Both tabulated flux families ignore the status from left/right `ghl_tabulated_compute_eps_T_from_P` calls. Global rho/Ye/pressure bounds do not guarantee attainable pressure at a particular rho/Ye within the permitted temperature range. The inversion can return before assigning epsilon.

- **Production trigger and likelihood**: A reconstructed face pressure cannot be inverted at the selected rho/Ye. Status loss is explicit in both families; no failing face state with an actual supported table was executed, so its frequency is unestablished.

- **Impact**: The failed inverse is not handled at the call site. Later characteristic-speed/HLLE work recomputes pressure and epsilon from temperature, allowing reconstructed pressure to be replaced or a later abort. Inevitable consumption of unset epsilon by the final flux is not established.

- **Locations**: `GRHayLHD/src/Tabulated/evaluate_fluxes_rhs.c:116`; `GRHayLHD/src/TabulatedEntropy/evaluate_fluxes_rhs.c:118`; `GRHayL-main/GRHayL/include/ghl_eos_functions.h:193`; `GRHayL-main/GRHayL/EOS/Tabulated/interpolators/NRPyEOS_eps_and_T_from_rho_Ye_P.c:21`. Original audit ID: `F05`.

- **Repository evidence and checks**: Source tracing confirms discarded return values, early error return, and subsequent temperature-based thermodynamic recomputation. Existing scalar clamps do not establish fixed-rho/Ye inverse attainability. No real-table failure reproduction or numerical flux comparison was performed.

- **Confidence**: High for discarded status and the bounded dataflow; realistic failing-table incidence remains unverified.

- **Desired postcondition and fix direction**: Handle each inverse error before treating the reconstructed face as valid. Apply an explicit abort or defined reconstruction/fallback policy and preserve its status and thermodynamic consistency.

- **Verification**: Use a supported table with successful and unattainable-pressure faces. Inspect final pressure, temperature, epsilon, status handling, and flux for both entropy modes.

</details>

<details>
<summary><strong>Issue 6 (medium, API input initialization):</strong> Generic conservative helpers evaluate omitted inactive scalar members.</summary>

<br>

- **Description**: Generic conservative computation and undensitization unconditionally evaluate entropy and Ye members. Hybrid/Simple callers omit Ye; nonentropy callers omit entropy; corresponding conservative averages omit inactive members. Automatic structures are not fully initialized, and the primitive limiter does not fill every inactive scalar.

- **Production trigger and likelihood**: Generic conversions, recovery/reprojection, averaging, or boundary work in Hybrid/Simple, HybridEntropy, and nonentropy Tabulated modes. These call sites are present in supported modes. TabulatedEntropy supplies both active scalars and is not implicated solely by this omission.

- **Impact**: Indeterminate floating-point inputs are evaluated at the API boundary, with possible undefined/trapping behavior. Their inactive outputs are generally discarded, so this finding alone does not prove evolved mass/energy corruption. It is distinct from the actively consumed solver guesses in Issue 3.

- **Locations**: `GRHayLHD/src/Hybrid/prims_to_conservs.c:28`; `GRHayLHD/src/Hybrid/conservs_to_prims.c:74` and analogous family/boundary/reprojection sites; `GRHayL-main/GRHayL/Con2Prim/compute_conservs.c:97`; `GRHayL-main/GRHayL/Con2Prim/undensitize_conservatives.c:41`. Original audit ID: `F06`.

- **Repository evidence and checks**: Direct tracing confirms omitted member assignments and unconditional helper evaluations, rather than merely unspecified structure padding. No complete instrumented Cactus run or active-output corruption reproduction was performed.

- **Confidence**: High for the omitted-input mechanism; observable downstream effects are bounded by discarded inactive outputs and execution settings.

- **Desired postcondition and fix direction**: Initialize complete API structures to defined values before loading active state, including averaged conservatives and atmosphere/boundary/reprojection paths.

- **Verification**: Inspect complete inputs at every generic conversion/undensitization call for each EOS/entropy mode. Exercise normal, atmosphere, averaging, and boundary paths with suitable uninitialized-read instrumentation.

</details>

<details>
<summary><strong>Issue 7 (medium, API accumulator initialization):</strong> Two Hybrid RHS calls leave the speed-limit flag uninitialized.</summary>

<br>

- **Description**: Nonentropy Hybrid source and flux routines declare `bool speed_limited` without initialization. Current upstream limiting uses `*speed_limited |= true` when clipping activates, so the API treats the incoming flag as an accumulator. Other families initialize it to false.

- **Production trigger and likelihood**: Nonentropy Hybrid/Simple RHS work activates velocity limiting. The uninitialized call sites are established; clipping frequency and observable consequences were not measured.

- **Impact**: An indeterminate boolean participates in the source-level compound assignment. The local result is unused, and optimization may eliminate the machine read. Numerical RHS damage is not established from this flag alone; the upstream accumulation contract itself is not a defect.

- **Locations**: `GRHayLHD/src/Hybrid/evaluate_sources_rhs.c:48`; `GRHayLHD/src/Hybrid/evaluate_fluxes_rhs.c:112`; `GRHayL-main/GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:38`. Original audit ID: `F07`.

- **Repository evidence and checks**: Source inspection confirms uninitialized caller storage and the accumulator expression. C11 compound-assignment semantics support the input obligation while leaving compiler optimization and observed damage separate. No numerical clipping-path failure was reproduced.

- **Confidence**: High for the caller initialization omission and bounded API contract.

- **Desired postcondition and fix direction**: Initialize the flag to false before the first call in each intended accumulation lifetime. Preserve the library's documented accumulation behavior.

- **Verification**: Inspect both call sites and flag lifetimes; activate clipping in a focused input check and verify initialized accumulator inputs.

</details>

<details>
<summary><strong>Issue 8 (medium, floating-point exception dependent):</strong> Metric derivative packing invokes a physical-metric inversion.</summary>

<br>

- **Description**: Derivatives of lapse, shift, and spatial metric are passed to `ghl_initialize_metric`, which computes reciprocal lapse and inverse metric as though its inputs formed a nonsingular physical metric. Derivative tensors have no such guarantee. Constant flat geometry supplies zero derivatives.

- **Production trigger and likelihood**: Source evaluation in any family, including constant Minkowski geometry. Exceptional arithmetic is reachable on valid geometry; an abort requires enabled floating-point traps or equivalent execution policy.

- **Impact**: Zero derivative lapse/determinant produces infinite reciprocal lapse, NaN inverse entries, and floating-point exceptions. Current source kernels use the raw derivative components, so corrupted RHS under masked IEEE exceptions is not established. Trap-enabled execution can abort.

- **Locations**: `GRHayLHD/src/compute_metric_derivs.c:41`; `GRHayL-main/GRHayL/GRHayL_Core/initialize_metric.c:60`; calls from all source families. Original audit ID: `F08`.

- **Repository evidence and checks**: An isolated C11/O0 probe included the actual thorn helper and upstream constructor. Zero derivatives produced `lapseinv=inf`, NaN inverse entries, and raised `FE_DIVBYZERO`/`FE_INVALID`. The declaration/index shim did not include a scheduler or full evolution run.

- **Confidence**: High for constructor misuse and reproduced exceptional arithmetic; default masked-exception RHS corruption is not claimed.

- **Desired postcondition and fix direction**: Pack the derivative components required by source kernels directly into defined storage without applying physical-metric inversion to them.

- **Verification**: Check constant and varying geometry with exception flags and traps. Compare source terms against independent derivative inputs and require unchanged valid results without spurious inversion exceptions.

</details>

<details>
<summary><strong>Issue 9 (medium, schedule metadata):</strong> CCL access declarations omit actual Tmunu and leakage-conversion accesses.</summary>

<br>

- **Description**: Tmunu accumulation uses `+=` for all ten components, but its scheduled block declares writes without the corresponding reads. The leakage converter instance also omits ADM metric reads and HydroBase Lorentz-factor writes present in the implementation and declared by other converter instances.

- **Production trigger and likelihood**: Tmunu addition or the leakage converter executes in a driver/PreSync configuration using access declarations for validity checking or data preparation. The declaration omissions are established; a concrete default-driver failure was not reproduced.

- **Impact**: Incomplete framework contracts can impair access checking, preparation/copies, or validity bookkeeping. Additive Tmunu accumulation itself is correct under the contributor contract.

- **Locations**: `GRHayLHD/schedule.ccl:109`; `GRHayLHD/src/compute_Tmunu.c:40`; `GRHayLHD/schedule.ccl:137`; `GRHayLHD/src/convert_GRHayLHD_to_HydroBase.c:58` and `:81`. Original audit ID: `F10`.

- **Repository evidence and checks**: Source/CCL comparison identifies the missing accesses. Primary [PreSync documentation](https://docs.einsteintoolkit.org/et-docs/PreSync) requires truthful accesses; the [TmunuBase guide](https://einsteintoolkit.org/thornguide/EinsteinBase/TmunuBase/documentation.html) supports addition by contributors. No deployed-driver validity/poisoning run was performed.

- **Confidence**: High for incomplete declarations; operational effects depend on the actual driver and enforcement mode.

- **Desired postcondition and fix direction**: Declare every actual read/write at each scheduled instance, or use a dedicated converter with a narrower truthful access contract.

- **Verification**: Inspect generated access declarations and run the deployed driver's validity/poisoning checks with Tmunu and leakage enabled.

</details>

<details>
<summary><strong>Issue 10 (medium, non-REAL8 conditional):</strong> Flux buffers and callback signatures assume CCTK_REAL is double.</summary>

<br>

- **Description**: All four flux families declare API-facing buffers, callback types, and characteristic outputs using `CCTK_REAL`. Current GRHayL interfaces use `double` and double pointers. Different pointer or function types cannot be repaired by ordinary scalar conversion.

- **Production trigger and likelihood**: A non-REAL8 Cactus configuration makes these interface types differ. Cactus offers precision choices where architecture support exists, but full GRHayLib/selected-thorn support for such builds was not established. The supplied REAL8 configuration is not implicated.

- **Impact**: Compilation can reject incompatible pointers/signatures; if mismatched calls proceed, buffer/ABI misuse is possible. Medium refers to this conditional API correctness concern, rather than a demonstrated default failure.

- **Locations**: All four `GRHayLHD/src/*/evaluate_fluxes_rhs.c` files; representative `GRHayLHD/src/Hybrid/evaluate_fluxes_rhs.c:23`, `:74`, and K/Gamma output addresses near line 6; `GRHayL-main/GRHayL/include/ghl_flux_source.h` and EOS/reconstruction headers. Original audit ID: `F11`.

- **Repository evidence and checks**: A syntax-only shim using actual thorn units/upstream headers accepted double for all units and rejected the four flux units for float and long double with incompatible-pointer diagnostics. It used independent Cactus arguments and did not validate CCL-generated bindings, linking, or an actual Cactus precision configuration. The primary [Cactus compilation guide](https://www.cactuscode.org/documentation/usersguide/UsersGuidech6.html) does not establish support by every thorn.

- **Confidence**: High for conditional type incompatibility; actual supported-stack build outcomes remain unverified.

- **Desired postcondition and fix direction**: Use `double` for API-facing local buffers/callbacks and convert at grid boundaries, or explicitly require REAL8 and clearly reject unsupported precision configurations.

- **Verification**: Build the precision configurations the complete stack claims to support against the supplied upstream version. Check generated bindings and actual linking, with clear rejection of unsupported configurations.

</details>

<details>
<summary><strong>Issue 11 (low, diagnostic):</strong> Final atmosphere-reset attribution is overwritten.</summary>

<br>

- **Description**: Every recovery family adds 100 directly to `failure_checker[index]` after all attempts fail, then overwrites the field from `local_failure_checker`, which never received that contribution. Separately, both logged `InHoriz` counts increment within the final-failure branch, so their denominator is not the size of the entire threshold region.

- **Production trigger and likelihood**: Final recovery failure invokes the atmosphere reset. The write/overwrite sequence is established in all families; reset frequency was not measured.

- **Impact**: The per-point diagnostic loses final-reset attribution, and a whole-region failure fraction cannot be inferred from the logged denominator. Atmosphere resetting, aggregate failure counting, primitive limiting, and conservative reprojection still execute. No inspected evolution decision consumes the grid diagnostic; the demonstrated loss is diagnostic-only.

- **Locations**: `GRHayLHD/interface.ccl:60`; reset/final-assignment pairs in `GRHayLHD/src/Hybrid/conservs_to_prims.c:200` and `:237`, HybridEntropy `:210` and `:248`, Tabulated `:178` and `:217`, TabulatedEntropy `:184` and `:224`. Original audit ID: `F09`.

- **Repository evidence and checks**: Direct source tracing and a scoped search of `failure_checker` uses confirm the overwrite, diagnostic declaration, preserved reset/count/reprojection, and absence of an inspected numerical evolution consumer. No forced-reset runtime probe was performed. The false-positive review corrected the original Medium ranking to Low.

- **Confidence**: High for the diagnostic mechanism and corrected bounded severity.

- **Desired postcondition and fix direction**: Accumulate the reset marker locally and make one final grid assignment, avoiding the unnecessary read of previous scratch contents. Label or recompute the horizon denominator according to its intended meaning.

- **Verification**: Force a final atmosphere reset in each family. Check the spatial marker, unchanged aggregate/reset behavior, and the independently calculated regional logging denominator.

</details>

<details>
<summary><strong>Issue 12 (low, conditional synchronization):</strong> SYNC declarations name optional groups without storage.</summary>

<br>

- **Description**: `GRHayLHD_sync_conservatives` always names mass/energy/momentum, Ye, and entropy conservative groups. Ye storage exists only for Tabulated EOS and entropy storage only when entropy evolves, so supported modes request synchronization of unavailable optional groups.

- **Production trigger and likelihood**: An EOS/entropy mode lacks one optional group and the deployed driver processes this explicit SYNC list. Actual behavior depends on driver/version/mode; warning visibility is not guaranteed.

- **Impact**: An invalid synchronization request can produce warning/failure status. Current public Carpet skips unavailable groups while handling valid ones; fatal default evolution failure is not established.

- **Locations**: `GRHayLHD/schedule.ccl:83`; conditional storage at `GRHayLHD/schedule.ccl:7`. Original audit ID: `F12`.

- **Repository evidence and checks**: Source declarations establish storage/SYNC disagreement. Primary [Carpet communication source](https://raw.githubusercontent.com/EinsteinToolkit/Carpet/master/Carpet/src/Comm.cc) issues a level-4 warning, reports negative status, and retains available groups; [scheduled-call source](https://raw.githubusercontent.com/EinsteinToolkit/Carpet/master/Carpet/src/CallFunction.cc) bounds explicit SYNC processing by mode. The deployed driver was unavailable, and no synchronization run was performed.

- **Confidence**: High for the invalid request; operational consequences are bounded to the inspected implementation and actual mode.

- **Desired postcondition and fix direction**: Condition optional SYNC entries on the same EOS/entropy predicates as storage, leaving required conservative synchronization intact.

- **Verification**: Inspect generated schedules for each supported EOS/entropy combination. Verify synchronization of allocated groups and actual missing-storage warning/status behavior in the deployed driver.

</details>

<details>
<summary><strong>Issue 13 (low, documentation):</strong> Conversion comments, the public guide, and boundary banners misstate behavior.</summary>

<br>

- **Description**: Both converter comments give `alpha*u0=1/sqrt(1+gamma^ij*u_i*u_j)` instead of the correct positive square root. The guide initially says default behavior prevents HydroBase diagnostic use by never copying data back, although thermodynamic scalars are already shared HydroBase fields. Boundary banners describe magnetic A/B operations absent from this HD thorn.

- **Production trigger and likelihood**: A maintainer or user follows these comments or documentation. The guide later clarifies directly shared scalars and distinct velocity fields, narrowing that particular defect to its overbroad opening language.

- **Impact**: Incorrect derivation and descriptions can misdirect implementation or configuration. Executable Valencia/transport conversion uses the correct Lorentz relation; no executable sign error or unavailable default density/pressure diagnostics is established.

- **Locations**: `GRHayLHD/doc/documentation.tex:62`; `GRHayLHD/src/convert_HydroBase_to_GRHayLHD.c:31`; `GRHayLHD/src/convert_GRHayLHD_to_HydroBase.c:50`; opening banners in all four `GRHayLHD/src/*/outer_boundaries.c` files. Original audit ID: `F14`.

- **Repository evidence and checks**: Four-velocity normalization gives `alpha*u0=sqrt(1+gamma^ij*u_i*u_j)`. Comments, executable assignments, the complete guide paragraph, and owned HD fields were compared directly. The [HydroBase guide](https://einsteintoolkit.org/thornguide/EinsteinBase/HydroBase/documentation.html) supplies the velocity contract. No numerical conversion failure is claimed from these prose errors.

- **Confidence**: High for the localized math/documentation errors and bounded impact.

- **Desired postcondition and fix direction**: Correct the identity, name the velocity/Lorentz-factor fields needing optional refresh and their consumers, and replace magnetic banners with the actual hydro boundary sequence. Describe actual EOS/entropy solver restrictions accurately.

- **Verification**: Compare the revised statements with normalization, assignments, shared-field ownership, CCL selection, and current upstream restrictions. Check the complete guide context for consistent diagnostic guidance.

</details>
