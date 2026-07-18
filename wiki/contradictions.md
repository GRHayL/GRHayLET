# IllinoisGRMHD Contradictions

> Active source disagreements and bounded authority decisions. · Status: confirmed · Last reconciled: 07-17-2026

## Register

| ID | Claim | Claim status | Source A | Source B | Authority decision | Affected pages | Page-status rationale | Owner/trigger | Resolution test | Opened | Resolved | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `CONTR-0001` | Documentation says backward-compatibility support ends in `ET_2024_11`, while current tree retains build, declaration, parameter, schedule, and C implementation surfaces. | stale | [`documentation.tex`](../IllinoisGRMHD/doc/documentation.tex), `Updating Old Parfiles` | [`interface.ccl`](../IllinoisGRMHD/interface.ccl), backward-compatibility groups; [`param.ccl`](../IllinoisGRMHD/param.ccl), deprecated parameters/options; [`schedule.ccl`](../IllinoisGRMHD/schedule.ccl), `ID_converter_ILGRMHD` branch; [`src/make.code.defn`](../IllinoisGRMHD/src/make.code.defn), `SRCS`; [`backward_compatible_data.c`](../IllinoisGRMHD/src/backward_compatible_data.c) and [`backward_compatible_initialize.c`](../IllinoisGRMHD/src/backward_compatible_initialize.c), named functions | Current code/build/config decide surfaces retained by this tree. They do not decide support policy for an external Einstein Toolkit release. | [`integration/migration-and-backward-compatibility.md`](integration/migration-and-backward-compatibility.md) | Page may remain confirmed because current-tree presence is central and directly supported; stale sunset wording is bounded, marked, and carries no external support guarantee. | Migration owner; trigger on migration prose, deprecated CCL surface, compatibility schedules/functions, or common build list change. | Obtain maintainer/release decision, then deterministically reconcile ThornGuide, CCL, build list, both C functions, reverse dependents, and page marker. | 07-17-2026 | - | Do not infer whether any external release supports these surfaces. |
| `CONTR-0002` | Terminal recovery executes `failure_checker[index] += 100`, but a later unconditional assignment omits that increment in all four variants; decoder text promises a hundreds marker. | contested | [`Hybrid`](../IllinoisGRMHD/src/Hybrid/conservs_to_prims.c), [`HybridEntropy`](../IllinoisGRMHD/src/HybridEntropy/conservs_to_prims.c), [`Tabulated`](../IllinoisGRMHD/src/Tabulated/conservs_to_prims.c), and [`TabulatedEntropy`](../IllinoisGRMHD/src/TabulatedEntropy/conservs_to_prims.c) terminal fallback increments and decoder comments | Same four functions, later `failure_checker[index] = local_failure_checker + 1000*diagnostics.backup[0] + 10000*diagnostics.tau_fix + 100000*diagnostics.Stilde_fix` assignments | Final unconditional write decides current observable local dataflow: atmosphere reset and counters remain, but this code does not retain terminal `100`. Tabulated variants also have no local Font1D call despite decoder wording. | [`evolution/con2prim-recovery-and-diagnostics.md`](evolution/con2prim-recovery-and-diagnostics.md) | Contested issue changes page's central diagnostic answer; page stays contested until code, decoder, and targeted observation agree. | Recovery owner; trigger on any variant recovery ladder, `failure_checker` writes, decoder, or diagnostic schedule/interface change. | Inspect all four variants and require terminal failure in final assigned value plus path-accurate decoder wording; then run an authorized targeted case forcing terminal fallback and observe the hundreds digit. | 07-17-2026 | - | Static inspection only; no runtime result claimed. |

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

Each Hybrid, HybridEntropy, Tabulated, and TabulatedEntropy
`conservs_to_prims.c` initializes `local_failure_checker` to zero. Terminal
failure increments the grid function directly by `100`; later in the same loop
iteration an unconditional assignment rebuilds that grid value from
`local_failure_checker` and diagnostic terms without the terminal increment.
Hybrid variants locally try `ghl_hybrid_Font1D`; Tabulated variants go directly
from exhausted neighbor averaging to atmosphere, although all decoder comments
say “C2P and Font Fix failed.”

Claim evidence:
- Claim: In all four local recovery variants, the final unconditional assignment overwrites the earlier terminal `failure_checker += 100`; static inspection does not establish an observed runtime value.
- Role: descriptive behavior
- Deciding authority: registered four `IllinoisGRMHD/src/*/conservs_to_prims.c` files, terminal fallback branches and final `failure_checker[index]` assignments
- Corroboration: same functions' decoder comments expose the intended hundreds marker and Tabulated wording mismatch
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=inspected-not-run; options=Hybrid,HybridEntropy,Tabulated,TabulatedEntropy; date=07-17-2026`

## Known Gaps, Not Contradictions

- `TEST TOV` in `test/test.ccl` and `test/magnetizedTOV` coexist, but local files
  do not define Cactus test-discovery mapping. No pass/failure inference.
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
