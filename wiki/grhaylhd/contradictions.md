# GRHayLHD Issues

> Page status: reviewed · Last reviewed: 07-17-2026

Rows below are evidence-backed candidates. Affected Page IDs may precede their
pages during unpublished construction. Once an affected page exists, it must
link the exact issue anchor. Safe wording remains narrower than any runtime
conclusion.

| ID | Kind | Status | Claim / ambiguity | Locator A | Locator B | Affected Page IDs | Safe wording | Impact | Owner / trigger | Resolution test | Resolution locator | Opened | Resolved |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `GRH-0001` | mismatch | open | Parameter permits only `none`, while initialization contains an equatorial branch. | `ccl:GRHayLHD/param.ccl#parameter=Symmetry` | `c:GRHayLHD/src/InitSymBound.c#symbol=GRHayLHD_InitSymBound` | `grhaylhd.evolution.matter-boundaries-and-symmetry`; `grhaylhd.integration.parameters-and-configurations` | Only `none` is locally selectable; equatorial code is dormant. | Prevents unsupported symmetry claims. | Parameter or implementation changes | Add matching local selectable declaration and statically reconcile every branch. | - | 07-17-2026 | - |
| `GRH-0002` | mismatch | open | Some recovery variants request `grhd_conservatives`; interface declares `grmhd_conservatives`. | `ccl:GRHayLHD/interface.ccl#group=grmhd_conservatives` | `c:GRHayLHD/src/Tabulated/conservs_to_prims.c#symbol=GRHayLHD_tabulated_conservs_to_prims` | `grhaylhd.architecture.variables-and-storage`; `grhaylhd.evolution.conservative-recovery`; `grhaylhd.evolution.matter-boundaries-and-symmetry`; `grhaylhd.integration.adm-mol-tmunu-contracts` | Three visible variant routines contain a group-name mismatch; runtime effect is unverified. | May affect symmetry registration paths. | Interface or variant recovery changes | Make local names consistent, then inspect all four routines and declarations. | - | 07-17-2026 | - |
| `GRH-0003` | mismatch | open | Initialization names `Stilde_z`; interface declares `Stildez`. | `ccl:GRHayLHD/interface.ccl#group=grmhd_conservatives` | `c:GRHayLHD/src/InitSymBound.c#symbol=GRHayLHD_InitSymBound` | `grhaylhd.evolution.matter-boundaries-and-symmetry`; `grhaylhd.architecture.variables-and-storage` | Visible strings differ; no runtime outcome is asserted. | May affect dormant equatorial registration. | Interface or initialization changes | Reconcile variable spelling and inspect generated/local declarations. | - | 07-17-2026 | - |
| `GRH-0004` | mismatch | open | All four terminal recovery branches increment `failure_checker`, then each visible point writeback assigns a separate local value. | `c:GRHayLHD/src/Hybrid/conservs_to_prims.c#symbol=GRHayLHD_hybrid_conservs_to_prims` | `ccl:GRHayLHD/interface.ccl#group=failure_checker` | `grhaylhd.architecture.variables-and-storage`; `grhaylhd.evolution.conservative-recovery`; `grhaylhd.evolution.perturbations-and-diagnostics` | Intended legend and final visible assignment are documented separately. | Diagnostic interpretation may be ambiguous. | Recovery writeback or legend changes | Inspect all four recovery routines; preserve terminal contribution through final assignment or revise legend with local evidence. | - | 07-17-2026 | - |
| `GRH-0005` | mismatch | open | Boundary headers describe magnetic stages absent from declared local state and lifecycle. | `doc:GRHayLHD/README#section=1. Purpose` | `c:GRHayLHD/src/Hybrid/outer_boundaries.c#symbol=GRHayLHD_hybrid_outer_boundaries` | `grhaylhd.architecture.purpose-build-surface`; `grhaylhd.evolution.matter-boundaries-and-symmetry` | Header text is stale or unexplained; no magnetic branch is inferred. | Protects non-magnetic scope. | Boundary documentation or state changes | Reconcile all four headers with local declarations and scheduled operations. | - | 07-17-2026 | - |
| `GRH-0006` | hazard | open | Leakage scheduling is outside the positive cadence guard, while converter uses modulo cadence. | `ccl:GRHayLHD/schedule.ccl#schedule=convert_GRHayLHD_to_HydroBase?context=GRHayLHD_RHS` | `c:GRHayLHD/src/convert_GRHayLHD_to_HydroBase.c#symbol=convert_GRHayLHD_to_HydroBase` | `grhaylhd.architecture.schedule-lifecycle`; `grhaylhd.integration.hydrobase-velocity-conversion`; `grhaylhd.integration.parameters-and-configurations` | Static precondition mismatch is visible; no runtime crash is claimed. | Zero default cadence may reach an unverified path. | Schedule guard or converter guard changes | Establish a positive local precondition on every call path or handle zero before modulo. | - | 07-17-2026 | - |
| `GRH-0007` | lifecycle-ambiguity | open | `update_Tmunu` is always steerable, while scheduling and registration are setup-conditional. | `ccl:GRHayLHD/param.ccl#parameter=update_Tmunu` | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_compute_Tmunu` | `grhaylhd.integration.adm-mol-tmunu-contracts`; `grhaylhd.integration.parameters-and-configurations` | Runtime steering effect is externally unresolved. | Registration and scheduled contribution may not share lifecycle. | Parameter or schedule/registration changes | Provide local lifecycle semantics or a test showing supported steering transitions. | - | 07-17-2026 | - |
| `GRH-0008` | hazard | open | Three Tmunu constrained-registration return values are not visibly accumulated into `ierr` like surrounding calls. | `c:GRHayLHD/src/MoL_registration.c#symbol=GRHayLHD_RegisterVars` | `ccl:GRHayLHD/param.ccl#parameter=update_Tmunu` | `grhaylhd.integration.adm-mol-tmunu-contracts` | Three unchecked-status review candidates are visible; no observed failure is claimed. | Registration errors may be less visible. | Registration implementation changes | Accumulate or explicitly handle all three return values, then statically inspect final error handling. | - | 07-17-2026 | - |
| `GRH-0009` | provenance-ambiguity | open | Balsara oracle header and companion generation header make distinct path assertions without proving a shared chain. | `oracle:GRHayLHD/test/Balsara0/rho.x.asc#file` | `par:GRHayLHD/test/Balsara0/Balsara0.par#file` | `grhaylhd.integration.parameters-and-configurations`; `grhaylhd.validation.test-inventory-and-oracles`; `grhaylhd.validation.coverage-gaps` | Assertions do not prove identity with current authored input, companion-to-oracle production, chronology, or a shared chain. | Limits reproducibility claims. | New provenance evidence | Supply an admissible local generation record tying exact current input, companion, oracle, and artifact chronology. | - | 07-17-2026 | - |
| `GRH-0010` | provenance-ambiguity | open | TOV companion asserts generation and original path while oracle headers lack comparable detail. | `oracle:GRHayLHD/test/TOV/hydrobase-rho.x.asc#file` | `par:GRHayLHD/test/TOV/TOV.par#file` | `grhaylhd.integration.parameters-and-configurations`; `grhaylhd.validation.test-inventory-and-oracles`; `grhaylhd.validation.coverage-gaps` | Checked-in files do not prove identity with current authored input, companion-to-oracle production, chronology, or a shared chain. | Limits reproducibility claims. | New provenance evidence | Supply an admissible local generation record tying exact current input, companion, oracle, and artifact chronology. | - | 07-17-2026 | - |
| `GRH-0011` | mismatch | open | Leakage schedule omits metric reads and the `w_lorentz` write used by the shared converter body. | `ccl:GRHayLHD/schedule.ccl#schedule=convert_GRHayLHD_to_HydroBase?context=GRHayLHD_RHS` | `c:GRHayLHD/src/convert_GRHayLHD_to_HydroBase.c#symbol=convert_GRHayLHD_to_HydroBase` | `grhaylhd.architecture.schedule-lifecycle`; `grhaylhd.integration.hydrobase-velocity-conversion` | Declared and visible field sets differ; runtime scheduler effect is unverified. | May hide dependencies or writes from framework tooling. | Schedule declaration or converter dataflow changes | Reconcile leakage `READS`/`WRITES` with the shared body, then inspect all three call contexts. | - | 07-17-2026 | - |
| `GRH-0012` | mismatch | open | Primitive perturbation schedules declare `HydroBase::eps` reads and writes, while four visible primitive perturbation bodies neither read nor assign `eps`. | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_hybrid_perturb_primitives` | `c:GRHayLHD/src/Hybrid/perturb_primitives.c#symbol=GRHayLHD_hybrid_perturb_primitives` | `grhaylhd.evolution.perturbations-and-diagnostics` | Declared and visible read/write sets differ; no runtime scheduler effect is asserted. | May overstate dependencies and writes to framework tooling. | Schedule declaration or perturbation body changes | Reconcile all four primitive `READS`/`WRITES` sets with their bodies, then inspect every family. | - | 07-17-2026 | - |
| `GRH-0013` | mismatch | open | All four Prim2Con schedules declare `HydroBase::eps` reads, while none of four visible bodies loads `eps[index]` into `prims`; every body later writes `prims.eps`. | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_hybrid_prims_to_conservs` | `c:GRHayLHD/src/Hybrid/prims_to_conservs.c#symbol=GRHayLHD_hybrid_prims_to_conservs` | `grhaylhd.evolution.primitive-conservative-conversion` | Declared and visible read sets differ; no runtime scheduler or numerical effect is asserted. | May overstate a framework dependency while leaving limited-output provenance ambiguous. | Prim2Con declaration or body changes | Reconcile `eps` input handling and inspect all four schedules and bodies. | - | 07-17-2026 | - |
| `GRH-0014` | mismatch | open | Tabulated recovery comments define code 100 as both C2P and Font Fix failure, although neither Tabulated routine contains an explicit local Font1D call; Hybrid families do. | `c:GRHayLHD/src/Hybrid/conservs_to_prims.c#symbol=GRHayLHD_hybrid_conservs_to_prims` | `c:GRHayLHD/src/Tabulated/conservs_to_prims.c#symbol=GRHayLHD_tabulated_conservs_to_prims` | `grhaylhd.evolution.conservative-recovery` | Source legend and visible variant fallback structure differ; no external-library or runtime diagnostic behavior is asserted. | May mislead interpretation of code 100 for Tabulated variants. | Recovery legend or fallback structure changes | Reconcile legend with all four visible fallback paths and inspect every decoder use. | - | 07-17-2026 | - |

### GRH-0001

Selectable symmetry surface versus dormant branch.

### GRH-0002

Conservative-group spelling across recovery variants.

### GRH-0003

Momentum-variable spelling in initialization.

### GRH-0004

Recovery diagnostic legend versus visible writeback.

### GRH-0005

Boundary header scope versus local hydrodynamic state.

### GRH-0006

Conversion cadence precondition.

### GRH-0007

Tmunu steering and setup lifecycle.

### GRH-0008

Tmunu constrained-registration status handling.

### GRH-0009

Balsara companion and oracle provenance.

### GRH-0010

TOV companion and oracle provenance.

### GRH-0011

Leakage converter declaration versus visible field use.

### GRH-0012

Primitive perturbation declared `eps` reads/writes versus visible bodies.

### GRH-0013

Prim2Con declared `eps` reads versus visible body input setup.

### GRH-0014

Recovery code-100 legend versus visible variant fallback structure.
