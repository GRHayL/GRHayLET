# GRHayLID Issues

> Page status: reviewed · Last reviewed: 07-19-2026

Rows below are evidence-backed candidates. Affected Page IDs may precede their
pages during unpublished construction. Once an affected page exists, it must
link the exact issue anchor. Safe wording remains narrower than any runtime
conclusion.

| ID | Kind | Status | Claim / ambiguity | Locator A | Locator B | Affected Page IDs | Safe wording | Impact | Owner / trigger | Resolution test | Resolution locator | Opened | Resolved |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `GID-0001` | mismatch | open | ThornGuide says entropy is controlled by a `compute_entropy` parameter, but that local parameter is absent; declarations instead expose `initial_entropy="GRHayLID"` and `EOS_type` dispatch. | `doc:GRHayLID/doc/documentation.tex#section=Parameters` | `ccl:GRHayLID/param.ccl#parameter=initial_entropy` | `grhaylid.initial-data.entropy-computation`; `grhaylid.integration.hydrobase-keyword-extensions`; `grhaylid.integration.parameters-and-configurations` | Documentation wording and checked-in control declarations differ; framework acceptance and execution are unverified. | Can misdirect configuration attempts. | Documentation or control declarations change | Reconcile documented control with every local declaration and schedule guard. | - | 07-19-2026 | - |
| `GID-0002` | mismatch | open | ThornGuide Introduction says only Avec is set and B variables are not set, while the magnetic body visibly assigns both Avec and Bvec. | `doc:GRHayLID/doc/documentation.tex#section=Introduction` | `c:GRHayLID/src/1D_tests_magnetic_data.c#symbol=GRHayLID_1D_tests_magnetic_data` | `grhaylid.initial-data.one-d-tests-magnetic`; `grhaylid.architecture.purpose-and-build-surface` | Documentation and visible local writes differ; scheduled execution is not asserted. | Can misstate magnetic output surface. | Documentation, schedule, or magnetic body changes | Reconcile documentation with declared and visible magnetic writes. | - | 07-19-2026 | - |
| `GID-0003` | mismatch | open | ThornGuide abstract advertises a cylindrical explosion, while its locatable setup enumeration, parameter values, schedules, and source units provide no matching local setup. | `doc:GRHayLID/doc/documentation.tex#section=Introduction` | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` | `grhaylid.architecture.purpose-and-build-surface`; `grhaylid.initial-data.one-d-tests-hydro` | The abstract phrase has no direct heading locator; Locator A anchors the locatable setup enumeration that omits it. | Can imply an unsupported checked-in selection. | Documentation or setup surface changes | Add a complete local setup surface or remove the abstract claim, then re-audit all selectors. | - | 07-19-2026 | - |
| `GID-0004` | mismatch | open | `stagger_A_fields` declares half-cell staggering intent and defaults to `yes`, but no checked-in C unit reads it; `x_stag`, `y_stag`, and `z_stag` are assigned directly from `x[index]`, `y[index]`, and `z[index]` with no visible offset arithmetic, while a comment asserts a staggered convention. | `ccl:GRHayLID/param.ccl#parameter=stagger_A_fields` | `c:GRHayLID/src/1D_tests_magnetic_data.c#symbol=GRHayLID_1D_tests_magnetic_data` | `grhaylid.initial-data.one-d-tests-magnetic`; `grhaylid.integration.parameters-and-configurations`; `grhaylid.validation.coverage-gaps` | Declaration, checked-in consumer absence, visible expressions, and comment are recorded separately; numerical placement is unverified. | Leaves staggering intent without implementation evidence. | Parameter or magnetic implementation changes | Add a visible consumer and checked evidence of intended coordinate placement, or revise the declaration/comment. | - | 07-19-2026 | - |
| `GID-0005` | mismatch | open | `initial_data_1D` permits `sound wave` and the point loop contains a sine branch, but its pre-loop dispatch arm is commented out; the visible live chain reaches the invalid-name error arm before the loop. | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` | `c:GRHayLID/src/1D_tests_hydro_data.c#symbol=GRHayLID_1D_tests_hydro_data` | `grhaylid.initial-data.one-d-tests-hydro`; `grhaylid.validation.coverage-gaps` | Visible source ordering only; no claim depends on external error-macro termination semantics. | Selection declaration and visible dispatch disagree. | Parameter or hydro dispatch changes | Restore a live matching dispatch arm or remove the selection and loop branch, then inspect the whole chain. | - | 07-19-2026 | - |
| `GID-0006` | hazard | open | Sound-wave loop branch assigns `press = 1.0` beside the comment `should add kinetic energy here`; this loop-side branch is textually ordered after the invalid-name `CCTK_VERROR` arm described by GID-0005, and its reachability is unverified. | `c:GRHayLID/src/1D_tests_hydro_data.c#symbol=GRHayLID_1D_tests_hydro_data` | `ccl:GRHayLID/param.ccl#parameter=wave_amplitude` | `grhaylid.initial-data.one-d-tests-hydro`; `grhaylid.validation.coverage-gaps` | Source contains an explicit incompleteness note; no numerical defect or runtime path is asserted. | Marks an incomplete visible formula if dispatch is restored. | Sound-wave implementation changes | Define the intended local expression and add checked evidence for the selectable branch. | - | 07-19-2026 | - |
| `GID-0007` | hazard | open | In `GRHayLID_BetaEquilibrium`, the table-bounds comparison using `beq_temperature` visibly precedes `CHECK_PARAMETER(beq_temperature)`, so the unset sentinel reaches the comparison first in source order. | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_BetaEquilibrium` | `macro:GRHayLID/src/GRHayLID.h#name=CHECK_PARAMETER` | `grhaylid.initial-data.beta-equilibrium`; `grhaylid.integration.parameters-and-configurations` | Visible order only; reachability and the comparison result depend on external values and error semantics. | Sentinel validation may be later than readers expect. | Beta-equilibrium or macro use changes | Move or otherwise establish the sentinel check before local use, then inspect source order. | - | 07-19-2026 | - |
| `GID-0008` | mismatch | open | Beta-equilibrium schedule description names entropy among outputs, but declared WRITES and the visible body omit entropy while writing press, eps, Y_e, and temperature. | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_BetaEquilibrium` | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_BetaEquilibrium` | `grhaylid.initial-data.beta-equilibrium`; `grhaylid.architecture.schedule-lifecycle` | Schedule description differs from declared and visible write sets; execution is unverified. | Can misstate initialized state. | Schedule description/declaration or body changes | Make description, WRITES, and visible assignments agree, then inspect all outputs. | - | 07-19-2026 | - |
| `GID-0009` | mismatch | open | README says entropy can be computed for any EOS, while checked-in schedule declarations dispatch only Hybrid and Tabulated selections. | `doc:GRHayLID/README#section=1. Purpose` | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_hybrid` | `grhaylid.initial-data.entropy-computation`; `grhaylid.architecture.purpose-and-build-surface` | Documentation breadth exceeds the two locally declared dispatch arms; external EOS support is unverified. | Can overstate selectable entropy coverage. | README or entropy dispatch changes | Reconcile wording with all declared EOS selections or add a complete local dispatch arm. | - | 07-19-2026 | - |
| `GID-0010` | lifecycle-ambiguity | open | With `initial_entropy="GRHayLID"`, schedule declarations contain Hybrid and Tabulated arms but no arm or error declaration for another `EOS_type`. | `ccl:GRHayLID/param.ccl#parameter=initial_entropy` | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_compute_entropy_tabulated` | `grhaylid.initial-data.entropy-computation`; `grhaylid.architecture.schedule-lifecycle`; `grhaylid.validation.coverage-gaps` | The visible schedule text selects no entropy routine outside two named EOS values; runtime behavior is not asserted. | Leaves unsupported-selection handling undeclared. | Entropy schedule or shared keyword surface changes | Add a declared supported arm or explicit local handling, then inspect the complete conditional block. | - | 07-19-2026 | - |
| `GID-0011` | mismatch | open | IsotropicGas and ConstantDensitySphere schedule descriptions call the setups one-dimensional tests. | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_IsotropicGas` | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_ConstantDensitySphere` | `grhaylid.architecture.schedule-lifecycle`; `grhaylid.initial-data.isotropic-gas`; `grhaylid.initial-data.constant-density-sphere` | Two local schedule descriptions differ from the visible three-coordinate setup bodies. | Can misclassify both setup families. | Schedule descriptions or setup bodies change | Reconcile both descriptions with their visible coordinate and state construction. | - | 07-19-2026 | - |
| `GID-0012` | hazard | open | The magnetic routine's opening condition accepts either `initial_Avec="GRHayLID"` or `initial_Bvec="GRHayLID"`, while its point loop visibly writes both Avec and Bvec. | `c:GRHayLID/src/1D_tests_magnetic_data.c#symbol=GRHayLID_1D_tests_magnetic_data` | `ccl:GRHayLID/param.ccl#parameter=initial_Avec` | `grhaylid.initial-data.one-d-tests-magnetic`; `grhaylid.integration.hydrobase-keyword-extensions` | In-body selector guard and visible unconditional write surface differ; runtime scheduling is unverified. | One selected destination may imply writes to both. | Function guard or magnetic body changes | Reconcile the in-body guard, declared WRITES, and conditional assignments for both selectors. | - | 07-19-2026 | - |
| `GID-0013` | mismatch | open | README labels one-dimensional tests simple/hybrid only; the hydro guard excludes Tabulated and its error text says hybrid or ideal fluid EOS. | `doc:GRHayLID/README#section=1. Purpose` | `c:GRHayLID/src/1D_tests_hydro_data.c#symbol=GRHayLID_1D_tests_hydro_data` | `grhaylid.initial-data.one-d-tests-hydro` | Three local terms differ; no equivalence among EOS names or external types is inferred. | Can confuse selection requirements. | Documentation, guard, or error text changes | Establish one locally declared terminology and reconcile all three sites. | - | 07-19-2026 | - |
| `GID-0014` | mismatch | open | ThornGuide's two three-dimensional setup tables print `EOS_type` as `"tabulated"`, while both C bodies compare it with `"Tabulated"`. | `doc:GRHayLID/doc/documentation.tex#section=Parameters` | `c:GRHayLID/src/IsotropicGas.c#symbol=GRHayLID_IsotropicGas` | `grhaylid.initial-data.isotropic-gas`; `grhaylid.initial-data.constant-density-sphere` | Literal spellings differ; keyword normalization, acceptance, and execution are external. | Can misdirect EOS selection. | Documentation or both setup guards change | Reconcile the documented and compared literals, then inspect both setup bodies. | - | 07-19-2026 | - |
| `GID-0015` | mismatch | open | The `initial_hydro="HydroTest1D"` description directs users to `test_1D_initial_data`, but the declared local selector is `initial_data_1D`. | `ccl:GRHayLID/param.ccl#parameter=initial_hydro` | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` | `grhaylid.integration.hydrobase-keyword-extensions`; `grhaylid.integration.parameters-and-configurations` | Two checked-in parameter names differ; parser behavior and external aliases are unverified. | Can misdirect one-dimensional-test configuration. | Parameter descriptions or selector declarations change | Reconcile the referenced and declared parameter names, then inspect the complete selector surface. | - | 07-19-2026 | - |

### GID-0001

Entropy control documentation versus local keyword and schedule declarations.

### GID-0002

Documented magnetic output versus visible Avec and Bvec assignments.

### GID-0003

Abstract cylindrical-explosion wording versus checked-in setup surface.

### GID-0004

Declared staggering intent versus visible consumer and coordinate evidence.

### GID-0005

Sound-wave selection versus visible dispatch ordering.

### GID-0006

Sound-wave pressure incompleteness note.

### GID-0007

Beta-equilibrium temperature use versus sentinel-check order.

### GID-0008

Beta-equilibrium schedule description versus declared and visible writes.

### GID-0009

README any-EOS wording versus two-arm entropy schedule.

### GID-0010

Entropy selection outside the two declared EOS arms.

### GID-0011

Three-dimensional setup bodies with one-dimensional schedule descriptions.

### GID-0012

Either-keyword magnetic guard versus both-field visible writes.

### GID-0013

One-dimensional EOS terminology drift.

### GID-0014

ThornGuide EOS spelling versus both three-dimensional setup guards.

### GID-0015

HydroTest1D parameter-reference name versus the declared selector.
