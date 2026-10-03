# IllinoisGRMHD Contradictions

> Active source disagreements and bounded authority decisions. · Status: confirmed · Last reconciled: 10-02-2026

## Register

| ID | Claim | Claim status | Source A | Source B | Authority decision | Affected pages | Page-status rationale | Owner/trigger | Resolution test | Opened | Resolved | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `CONTR-0001` | Documentation says backward-compatibility support ends in `ET_2024_11`, while current tree retains build, declaration, parameter, schedule, and C implementation surfaces. | stale | [`documentation.tex`](../IllinoisGRMHD/doc/documentation.tex), `Updating Old Parfiles` | [`interface.ccl`](../IllinoisGRMHD/interface.ccl), backward-compatibility groups; [`param.ccl`](../IllinoisGRMHD/param.ccl), deprecated parameters/options; [`schedule.ccl`](../IllinoisGRMHD/schedule.ccl), `ID_converter_ILGRMHD` branch; [`src/make.code.defn`](../IllinoisGRMHD/src/make.code.defn), `SRCS`; [`backward_compatible_data.c`](../IllinoisGRMHD/src/backward_compatible_data.c) and [`backward_compatible_initialize.c`](../IllinoisGRMHD/src/backward_compatible_initialize.c), named functions | Current code/build/config decide surfaces retained by this tree. They do not decide support policy for an external Einstein Toolkit release. | [`integration/migration-and-backward-compatibility.md`](integration/migration-and-backward-compatibility.md) | Page may remain confirmed because current-tree presence is central and directly supported; stale sunset wording is bounded, marked, and carries no external support guarantee. | Migration owner; trigger on migration prose, deprecated CCL surface, compatibility schedules/functions, or common build list change. | Obtain maintainer/release decision, then deterministically reconcile ThornGuide, CCL, build list, both C functions, reverse dependents, and page marker. | 07-17-2026 | - | Do not infer whether any external release supports these surfaces. |
| `CONTR-0002` | Terminal recovery executes `failure_checker[index] += 100`, but a later unconditional assignment omits that increment in all four variants; decoder text promises a hundreds marker. | resolved | [`Hybrid`](../IllinoisGRMHD/src/Hybrid/conservs_to_prims.c), [`HybridEntropy`](../IllinoisGRMHD/src/HybridEntropy/conservs_to_prims.c), [`Tabulated`](../IllinoisGRMHD/src/Tabulated/conservs_to_prims.c), and [`TabulatedEntropy`](../IllinoisGRMHD/src/TabulatedEntropy/conservs_to_prims.c) terminal fallback increments and decoder comments | Same four functions, later `failure_checker[index] = local_failure_checker + 1000*diagnostics.backup[0] + 10000*diagnostics.tau_fix + 100000*diagnostics.Stilde_fix` assignments | All four variants add 100 to the local accumulator before final publication; decoder wording describes exhausted supported recovery. | [`evolution/con2prim-recovery-and-diagnostics.md`](evolution/con2prim-recovery-and-diagnostics.md) | Recovery page is confirmed after code/decoder reconciliation and disposable forced-failure execution. | Recovery owner; trigger on any variant recovery ladder, `failure_checker` writes, decoder, or diagnostic schedule/interface change. | Inspect all four variants and require terminal failure in final assigned value plus path-accurate decoder wording; then run an authorized targeted case forcing terminal fallback and observe the hundreds digit. | 07-17-2026 | 10-02-2026 | Controlled solver-failure harness executed actual four recovery bodies; final hundreds marker retained. No Cactus evolution or external solver validation. |

## CONTR-0001

The ThornGuide's `Updating Old Parfiles` section gives an `ET_2024_11` sunset.
Current `interface.ccl`, `param.ccl`, and `schedule.ccl` retain deprecated groups,
controls, and an `ID_converter_ILGRMHD` branch. Common `SRCS` still includes
`backward_compatible_data.c` and `backward_compatible_initialize.c`; those files
define their named data-copy and GRHayL-initialization functions. Thus current
tree presence is confirmed, while external release policy remains unsupported.

Claim evidence:
- Claim: This tree retains compiled and conditionally scheduled compatibility surfaces despite documented sunset wording; this does not guarantee support in any external release.
- Role: descriptive behavior
- Deciding authority: registered `IllinoisGRMHD/src/make.code.defn` `SRCS`; `IllinoisGRMHD/schedule.ccl` `ID_converter_ILGRMHD` branch; `IllinoisGRMHD/interface.ccl` backward-compatibility groups
- Corroboration: registered `IllinoisGRMHD/src/backward_compatible_initialize.c` and `backward_compatible_data.c`, named functions
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-run; options=ID_converter_ILGRMHD branch inspected; date=07-17-2026`

## CONTR-0002

Resolved on 10-02-2026: all four recovery bodies now increment
`local_failure_checker` by 100 on terminal atmosphere reset and publish that
accumulator once. Decoder wording names exhausted supported methods, including
Simple's explicit no-Font policy. A disposable controlled-failure caller harness
executed the actual four function bodies with full/truncated/empty neighborhoods
and observed 100 at every terminal-reset point. It used real GRHayL metric,
velocity, diagnostics, and undensitization helpers but controlled inversion/EOS
callbacks; this is not a Cactus schedule, solver, or evolution validation.

Claim evidence:
- Claim: Terminal recovery retains the hundreds marker in all four local bodies; bounded caller execution does not certify external inversion or numerical evolution.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/src/Hybrid/conservs_to_prims.c`, `IllinoisGRMHD_hybrid_conservs_to_prims`; registered `IllinoisGRMHD/src/HybridEntropy/conservs_to_prims.c`, `IllinoisGRMHD_hybrid_entropy_conservs_to_prims`; registered `IllinoisGRMHD/src/Tabulated/conservs_to_prims.c`, `IllinoisGRMHD_tabulated_conservs_to_prims`; registered `IllinoisGRMHD/src/TabulatedEntropy/conservs_to_prims.c`, `IllinoisGRMHD_tabulated_entropy_conservs_to_prims`, `local_failure_checker` accumulation and final publication
- Corroboration: disposable actual-body caller check described below; no checked-in runtime oracle available
- Validation: `inspected=pass; generated=not-run; built=pass; run=pass; result_checked=pass` (bounded C caller only; Cactus build and scheduling not run)
- Dimensions: `platform=Linux x86_64 CPU; tool_version=GCC Ubuntu 13.3.0; backend=disposable C caller with UBSan and controlled inversion/EOS callbacks; precision=double; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=forced all-method failure and terminal atmosphere reset; options=all four actual recovery bodies, 1-point and 27-point grids, empty/truncated/full neighborhoods; date=10-02-2026`

Executed check: `bash /tmp/tri-issues/run_recovery.sh`, working directory
`/tmp/tri-issues`, candidate recovery inputs directly compared with the
frozen repair snapshot. GCC compiled the actual recovery bodies against
inspected GRHayL headers with substitute CCTK declarations. The caller used
real metric/velocity/diagnostic/undensitization helpers and controlled
inversion/EOS callbacks. All assertions passed with exit 0: the final
published marker was 100 at every forced terminal-reset point in each
family and neighborhood case. Additional mixed-success checks observed
threshold-selected diagnostic populations of 0/27 and 9/27. This receipt
does not claim a generated Cactus build, driver entry validity, actual
inversion success, schedule execution, MPI/restart evolution, or
production-table/numerical validation.

## Known Gaps, Not Contradictions

- Balsara4 has example and top-level test parfiles, while its test block is
  commented and no per-case oracle directory is visible. This is a coverage
  gap.
- Public `Symmetry` permits only `none` and calls equatorial support in progress;
  dormant local equatorial branches contain parity setup. Unsupported public
  selection and partial dormant code are limitations, not competing claims.
- Visible cases select Simple or Hybrid, none explicitly Tabulated. They do not
  set shared `evolve_entropy`; its out-of-scope default leaves effective entropy
  family execution unknown.
- No initial static ingest demonstrates current build/test success, numerical
  convergence, divergence behavior, AMR/restart behavior, or external-library
  semantics.

## Rules

Follow [Schema](SCHEMA.md#contradiction-contract). Resolve only after the row's
test and all affected-page, reverse-dependency, catalog, alias, and typed-neighbor
reconciliation steps pass.
