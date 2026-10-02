# GRHayLET Evidence And Validation

Use applicable sections with [the protocol](review-protocol.md). Start at
[AGENTS.md](../AGENTS.md) and follow the owning branch; this reference adds no
workflow or budget.

## Repository And KB Scope

[GRHayLET](../README.md) contains Einstein Toolkit thorns that depend on external
GRHayL/GRHayLib: IllinoisGRMHD (GRMHD), GRHayLHD/GRHayLHDX (hydrodynamics),
GRHayLID/GRHayLIDX (initial data), and NRPyLeakageET (neutrino leakage). The `X`
variants target CarpetX; distinguish them from Carpet variants during review.
Library implementation, other Toolkit thorns, and configured Cactus environments
are external evidence, not assumed contents of this checkout.

Each KB branch admits only its matching source tree as domain evidence. A sibling
thorn or its KB is not authority for another branch. Use these owning routes:

| Branch / exclusive evidence | Contracts and maintenance | Impact and validation |
| --- | --- | --- |
| IllinoisGRMHD / `IllinoisGRMHD/**` | [Schema](../wiki/SCHEMA.md), [Workflows](../wiki/workflows.md), [Sources](../raw/SOURCES.md) | [Source Map](../wiki/source-map.md), [Issues](../wiki/contradictions.md), [Validation](../wiki/validation/index.md), [Checks](../wiki/lint/CHECKS.md) |
| GRHayLHD / `GRHayLHD/**` | [Index](../wiki/grhaylhd/index.md), [Schema](../wiki/grhaylhd/SCHEMA.md), [Workflows](../wiki/grhaylhd/workflows.md), [Sources](../raw/grhaylhd/SOURCES.md) | [Source Map](../wiki/grhaylhd/source-map.md), [Issues](../wiki/grhaylhd/contradictions.md), [Validation](../wiki/grhaylhd/validation/index.md), [Checks](../wiki/grhaylhd/lint/CHECKS.md) |
| GRHayLID / `GRHayLID/**` | [Index](../wiki/grhaylid/index.md), [Schema](../wiki/grhaylid/SCHEMA.md), [Workflows](../wiki/grhaylid/workflows.md), [Sources](../raw/grhaylid/SOURCES.md) | [Source Map](../wiki/grhaylid/source-map.md), [Issues](../wiki/grhaylid/contradictions.md), [Validation](../wiki/grhaylid/validation/index.md), [Checks](../wiki/grhaylid/lint/CHECKS.md) |

For GRHayLHDX, GRHayLIDX, NRPyLeakageET, and `bns_parfiles/`, inspect the requested
local sources directly. The existing KB branches do not admit those paths as
domain evidence. Cross-thorn code reviews may inspect relevant integration
dependencies, but must preserve branch isolation when compiling KB claims.

## Authority And Interfaces

Inspect owning C/C++ sources and headers, `interface.ccl`, `param.ccl`,
`configuration.ccl`, `schedule.ccl`, and recursive `src/**/make.code.defn`.
Thorn `README` and `doc/documentation.tex`, where present, supply documented intent;
current source/CCL decides current-tree declarations and visible behavior. KB pages
route and synthesize evidence; reopen exact sources for changed claims. Votes
cannot choose maintainer intent.

Trace parameter/mode selection, GRHayLib calls, HydroBase/ADM/MoL/Tmunu interfaces,
owned/shared variables, storage/time levels, schedule bins and ordering,
synchronization, boundaries, and symmetry where affected. Check driver-specific
loops, indexing, staggering, and CCTK interfaces from the actual variant.
Declarations and visible dataflow do not prove successful build, scheduled
execution, current test pass, numerical validity, or external-library semantics.

## Builds And Numerical Tests

This is a thorn collection, not the standalone GRHayL library build. Derive any
build or run command from the authorized, configured Cactus environment and its
selected thorns/driver/dependencies. Inspect CCL requirements and recursive source
manifests; do not import library `make tests`/`make datagen` or configure options.
Report unavailable build/runtime evidence explicitly.

Trace `test/test.ccl`, case parfiles, and checked-in `.asc`/`.tsv` observations where
present. Keep test declaration, shipped configuration, stored observation,
compile/link, execution, and oracle validity distinct. Check case names, active
versus commented declarations, process counts, tolerances, and mode coverage;
absence is a gap, not a passing test. Process success or self-consistency is not
numerical proof. Justify oracles/tolerances; do not loosen tolerances or regenerate
references just to pass.

KB maintenance permits static inspection and canonical lint, not Cactus builds,
test execution, external dependency execution, or fixture regeneration. Runtime
reproduction needs explicit authority already in the request/session, a configured
isolated environment, bounded resources, and recorded command, inputs, platform,
result, and scope. Preserve checked-in observations and parfiles.

## Documentation And External Dependencies

Review thorn documentation and KB claims against local code/CCL. Establish an
available, supported route before any authorized documentation compilation or
source regeneration; keep generated output in disposable storage and include
required companions. Do not assume standalone-library Doxygen or NRPy generator
infrastructure exists here.

GRHayLib headers/calls establish the local integration boundary, not external
implementation behavior or runtime compatibility. Scope any dependency change
separately under existing authority. Thorn declarations and checked-in observations
do not establish an Einstein Toolkit build, scheduler run, or current regression
pass. Inspect any applicable CI configuration actually present; configuration
presence proves selection, not historical execution.

## KB And Agent Instructions

Prepare one candidate. Follow the owning schema, workflows, registry, source map,
issues, and check contract. Reopen local ground truth for changed claims; update
affected routes, owner pages, catalog, glossary, reverse dependencies, locators,
and issue backlinks together. KB-only scope does not authorize source or external
dependency repairs.

Run `python tools/kb_lint.py` from the repository root. `--all` is an identical
compatibility alias, not stronger coverage; run it when the owning workflow requires
both commands. The canonical checker is deterministic and does not mutate source
trees. It proves structural contracts, not scientific truth or runtime success.
Checked-in KB-linter tests, regression harnesses, and fixture suites are forbidden
regardless of filename or location. Reproduce checker regressions only with
disposable manual fixtures outside the repository.

Never compute or store source fingerprints, digest values, or modification times.
Reconcile through paths, stable locators, and source-map reverse dependencies.
Retained KB dates use `MM-DD-YYYY`; routers use `n/a`. Only `plan1.md`, `plan2.md`,
`plan3.md`, `plan_synth.md`, `tasks1.md`, `tasks2.md`, `tasks3.md`, and `tasks4.md`
are coordination exemptions; keep them outside KB routes, catalog, source map,
and source registration. Operational checkpoints are task state, not KB evidence.

For agent instruction edits, check links/anchors, activation, precedence, coverage,
budgets, and exit/delivery paths; use scoped whitespace/link checks and canonical
lint, not unrelated numerical suites. Preserve the standing preference to skip
cross-model review unless explicitly requested and never ask whether to run it.
Mentioning, auditing, or editing either review skill does not activate its seats.
