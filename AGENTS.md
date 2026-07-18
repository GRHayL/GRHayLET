# GRHayLET Knowledge Base

This root selects between evidence-isolated knowledge-base branches. Each
branch admits only its matching source tree as domain evidence; neither branch
is authority for the other.

## Branch Selector

| Branch | Exclusive evidence scope | Start here |
| --- | --- | --- |
| IllinoisGRMHD | `IllinoisGRMHD/**` | [IllinoisGRMHD routes](#illinoisgrmhd-knowledge-base) |
| GRHayLHD | `GRHayLHD/**` | [GRHayLHD routes](wiki/grhaylhd/index.md) |

## GRHayLHD Router

| Go to | Use it for |
| --- | --- |
| [Architecture](wiki/grhaylhd/architecture/index.md) | Purpose, build surface, variables, storage, and declared schedule lifecycle. |
| [Evolution](wiki/grhaylhd/evolution/index.md) | EOS/entropy modes, conversion, recovery, RHS, matter boundaries, symmetry, perturbations, and diagnostics. |
| [Integration](wiki/grhaylhd/integration/index.md) | HydroBase, GRHayLib, ADM, MoL, Tmunu, parameters, and configurations. |
| [Validation](wiki/grhaylhd/validation/index.md) | Test declarations, parfiles, checked-in observations, provenance limits, and coverage gaps. |
| [Catalog](wiki/grhaylhd/catalog.md) | Exact namespaced page inventory and query aliases. |
| [Glossary](wiki/grhaylhd/glossary.md) | Canonical recurring terms and their owner pages. |
| [Workflows](wiki/grhaylhd/workflows.md) | Source registration, reconciliation, verification, and publication. |
| [Schema](wiki/grhaylhd/SCHEMA.md) | Page, claim, issue, locator, registry, and publication contracts. |
| [Sources](raw/grhaylhd/SOURCES.md) | Source identity, ownership, provenance class, lifecycle, and ingest state. |
| [Source Map](wiki/grhaylhd/source-map.md) | Canonical source-to-page evidence edges and next actions. |
| [Issues](wiki/grhaylhd/contradictions.md) | Open mismatches, hazards, and ambiguities. |
| [Lint Checks](wiki/grhaylhd/lint/CHECKS.md) | Deterministic checks and the manual proof boundary. |

## GRHayLHD Task Shortcuts

| Task | Read first |
| --- | --- |
| Understand purpose, build shape, variables, or declared scheduling | [Architecture](wiki/grhaylhd/architecture/index.md) |
| Change EOS/entropy modes, conversion, recovery, RHS, matter BCs, or perturbations | [Evolution](wiki/grhaylhd/evolution/index.md) |
| Use HydroBase/GRHayLib/ADM/MoL/Tmunu or inspect parameters | [Integration](wiki/grhaylhd/integration/index.md) |
| Inspect tests, tolerances, parfiles, or checked-in observations | [Validation](wiki/grhaylhd/validation/index.md) |
| Maintain the GRHayLHD branch | [Workflows](wiki/grhaylhd/workflows.md) |

---

## IllinoisGRMHD Knowledge Base

This plain-Markdown KB covers only `IllinoisGRMHD/**`. Local IllinoisGRMHD
code, CCL, build inputs, documentation, configurations, tests, and checked-in
oracles are evidence; sibling workspace trees are outside scope.

## Router

| Go to | Use it for |
| --- | --- |
| [Architecture](wiki/architecture/index.md) | Purpose, repository/build surface, and declared schedule lifecycle. |
| [Evolution](wiki/evolution/index.md) | EOS/entropy state, conversion, recovery, hydrodynamic RHS, boundaries, and perturbations. |
| [Magnetics](wiki/magnetics/index.md) | A/B staggering and reconstruction, induction/gauge RHS, EM boundaries, and symmetry. |
| [Integration](wiki/integration/index.md) | HydroBase/GRHayLib/Tmunu boundaries, parameters, and migration. |
| [Validation](wiki/validation/index.md) | Test declarations, oracles, shipped cases, and visible coverage gaps. |
| [Catalog](wiki/catalog.md) | Exact page inventory, aliases, and cross-branch query routing. |
| [Glossary](wiki/glossary.md) | Canonical recurring terms and their owner pages. |
| [Workflows](wiki/workflows.md) | Source registration, query, reconciliation, and safe verification. |
| [Schema](wiki/SCHEMA.md) | Page, authority, citation, status, and scope contracts. |
| [Sources](raw/SOURCES.md) | Source identity, provenance, lifecycle, and ingest state. |
| [Source Map](wiki/source-map.md) | Authority, reverse dependencies, coverage gaps, and next actions. |
| [Contradictions](wiki/contradictions.md) | Active source disagreements and resolution tests. |
| [Lint Checks](wiki/lint/CHECKS.md) | Deterministic checker contract and manual review boundary. |

## Task Shortcuts

| Task | Read first |
| --- | --- |
| Understand purpose, source/build shape, variables, or declared schedule | [Architecture](wiki/architecture/index.md) |
| Change EOS/entropy selection, Prim2Con, recovery, RHS, matter BCs, or perturbations | [Evolution](wiki/evolution/index.md) |
| Work on A/B staggering, magnetic reconstruction, induction/gauge RHS, EM BCs, or symmetry | [Magnetics](wiki/magnetics/index.md) |
| Use HydroBase/GRHayLib/Tmunu, inspect parameters, or migrate old parfiles | [Integration](wiki/integration/index.md) |
| Inspect tests, tolerances, oracles, shipped cases, or coverage | [Validation](wiki/validation/index.md) |
| Resolve an ambiguous or cross-branch term | [Catalog](wiki/catalog.md), then its branch router |
| Maintain this KB | [Workflows](wiki/workflows.md) |

## Source Tracking and Dates

- Retained KB dates use `MM-DD-YYYY`; routers use `n/a`.
- Never compute or store source fingerprints, digest values, or modification
  times. Reconcile drift through paths, stable locators, and reverse
  dependencies in the source map.
- Static inspection proves declarations and visible dataflow only. It does not
  prove a successful build, scheduled execution, current test pass, numerical
  validity, or external-library behavior.

## Checker and Coordination Scope

Run `python tools/kb_lint.py`. `--all` is an identical compatibility alias,
not stronger coverage. Checker is deterministic and never mutates
`IllinoisGRMHD/`.

Exact coordination exemptions are `plan1.md`, `plan2.md`, `plan3.md`,
`plan_synth.md`, `tasks1.md`, `tasks2.md`, `tasks3.md`, and `tasks4.md`. They
stay outside KB routes, catalog, source map, and source registration. No other
plan/task filename is implicitly exempt.
