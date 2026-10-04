# GRHayLHD Issues

> Page status: reviewed · Last reviewed: 10-04-2026

Rows below are evidence-backed candidates. Affected Page IDs may precede their
pages during unpublished construction. Once an affected page exists, it must
link the exact issue anchor. Safe wording remains narrower than any runtime
conclusion.

| ID | Kind | Status | Claim / ambiguity | Locator A | Locator B | Affected Page IDs | Safe wording | Impact | Owner / trigger | Resolution test | Resolution locator | Opened | Resolved |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `GRH-0001` | mismatch | open | Parameter permits only `none`, while initialization contains an equatorial branch. | `ccl:GRHayLHD/param.ccl#parameter=Symmetry` | `c:GRHayLHD/src/InitSymBound.c#symbol=GRHayLHD_InitSymBound` | `grhaylhd.evolution.matter-boundaries-and-symmetry`; `grhaylhd.integration.parameters-and-configurations` | Only `none` is locally selectable; equatorial code is dormant. | Prevents unsupported symmetry claims. | Parameter or implementation changes | Add matching local selectable declaration and statically reconcile every branch. | - | 07-17-2026 | - |
| `GRH-0007` | lifecycle-ambiguity | open | `update_Tmunu` is always steerable, while scheduling and registration are setup-conditional. | `ccl:GRHayLHD/param.ccl#parameter=update_Tmunu` | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_compute_Tmunu` | `grhaylhd.integration.adm-mol-tmunu-contracts`; `grhaylhd.integration.parameters-and-configurations` | Runtime steering effect is externally unresolved. | Registration and scheduled contribution may not share lifecycle. | Parameter or schedule/registration changes | Provide local lifecycle semantics or a test showing supported steering transitions. | - | 07-17-2026 | - |
| `GRH-0009` | provenance-ambiguity | open | Balsara oracle header and companion generation header make distinct path assertions without proving a shared chain. | `oracle:GRHayLHD/test/Balsara0/rho.x.asc#file` | `par:GRHayLHD/test/Balsara0/Balsara0.par#file` | `grhaylhd.integration.parameters-and-configurations`; `grhaylhd.validation.test-inventory-and-oracles`; `grhaylhd.validation.coverage-gaps` | Assertions do not prove identity with current authored input, companion-to-oracle production, chronology, or a shared chain. | Limits reproducibility claims. | New provenance evidence | Supply an admissible local generation record tying exact current input, companion, oracle, and artifact chronology. | - | 07-17-2026 | - |
| `GRH-0010` | provenance-ambiguity | open | TOV companion asserts generation and original path while oracle headers lack comparable detail. | `oracle:GRHayLHD/test/TOV/hydrobase-rho.x.asc#file` | `par:GRHayLHD/test/TOV/TOV.par#file` | `grhaylhd.integration.parameters-and-configurations`; `grhaylhd.validation.test-inventory-and-oracles`; `grhaylhd.validation.coverage-gaps` | Checked-in files do not prove identity with current authored input, companion-to-oracle production, chronology, or a shared chain. | Limits reproducibility claims. | New provenance evidence | Supply an admissible local generation record tying exact current input, companion, oracle, and artifact chronology. | - | 07-17-2026 | - |
| `GRH-0012` | mismatch | open | Primitive perturbation schedules declare `HydroBase::eps` reads and writes, while four visible primitive perturbation bodies neither read nor assign `eps`. | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_hybrid_perturb_primitives` | `c:GRHayLHD/src/Hybrid/perturb_primitives.c#symbol=GRHayLHD_hybrid_perturb_primitives` | `grhaylhd.evolution.perturbations-and-diagnostics` | Declared and visible read/write sets differ; no runtime scheduler effect is asserted. | May overstate dependencies and writes to framework tooling. | Schedule declaration or perturbation body changes | Reconcile all four primitive `READS`/`WRITES` sets with their bodies, then inspect every family. | - | 07-17-2026 | - |
| `GRH-0013` | mismatch | open | All four Prim2Con schedules declare `HydroBase::eps` reads, while none of four visible bodies loads `eps[index]` into `prims`; every body later writes `prims.eps`. | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_hybrid_prims_to_conservs` | `c:GRHayLHD/src/Hybrid/prims_to_conservs.c#symbol=GRHayLHD_hybrid_prims_to_conservs` | `grhaylhd.evolution.primitive-conservative-conversion` | Declared and visible read sets differ; no runtime scheduler or numerical effect is asserted. | May overstate a framework dependency while leaving limited-output provenance ambiguous. | Prim2Con declaration or body changes | Reconcile `eps` input handling and inspect all four schedules and bodies. | - | 07-17-2026 | - |

### GRH-0001

Selectable symmetry surface versus dormant branch.

### GRH-0007

Tmunu steering and setup lifecycle.

### GRH-0009

Balsara companion and oracle provenance.

### GRH-0010

TOV companion and oracle provenance.

### GRH-0012

Primitive perturbation declared `eps` reads/writes versus visible bodies.

### GRH-0013

Prim2Con declared `eps` reads versus visible body input setup.

