# IllinoisGRMHD KB Workflows

> Maintenance procedures. · Status: confirmed · Last reconciled: 07-17-2026

## Summary

Navigate root-first. Register local sources before compiling claims. Reconcile
living-source changes through source-map dependencies. Initial maintenance is
static and never runs Cactus or rewrites checked-in fixtures.

## Search and Query

1. Read [`AGENTS.md`](../AGENTS.md).
2. For ambiguous or cross-branch terms, inspect [`catalog.md`](catalog.md).
3. Enter one of five branch routers.
4. Read the owning leaf and typed neighbors.
5. Use scoped `rg` over `wiki/` if routing does not answer the question.
6. Use scoped `rg` over `IllinoisGRMHD/` only to verify or update facts.

Answer from a current compiled leaf when evidence suffices. If a normal query
needs broad source search, repair the owner leaf, alias, catalog row, or route.
Do not route through coordination artifacts or sibling source trees.

## Register and Ingest a Source

1. Add one unique row to [`../raw/SOURCES.md`](../raw/SOURCES.md), using one
   valid IllinoisGRMHD path/glob, lifecycle status, and ingest state.
2. Classify authority and add/update exactly one
   [`source-map.md`](source-map.md) row; do not copy ingest state there.
3. Identify one owning branch and leaf.
4. Inspect source fully, compile claims in prose, and cite stable locators.
5. Follow touched glossary entities and typed neighbors one hop.
6. Update reverse dependencies, catalog metadata, contradictions, and glossary
   owners where affected.
7. Promote `registered` to `partial` or `ingested` only after corresponding
   dependents are reconciled.
8. Run `python tools/kb_lint.py` from repository root and another CWD when path
   behavior changed. `--all` is an identical compatibility alias.

External web context is optional, labeled `background`, and non-deciding.
Record its URL in the work report. Every IllinoisGRMHD claim still needs a
registered local source.

## Durable Answer Filing

Do not file one-off answers. Update an existing leaf for recurring ambiguity,
a gotcha, or a durable single-owner answer. Add a narrow leaf only when content
has independent query value and an unambiguous branch. Add a synthesis branch
only when a recurring cross-branch answer has no clean owner. Update router,
catalog, source map, glossary, and typed neighbors with any structural change.

## Living-Source Reconciliation

1. Preserve and compare scoped workspace state; never erase user changes.
2. Identify changed IllinoisGRMHD paths and matching manifest IDs.
3. Follow source-map reverse dependents and stable locators.
4. Re-inspect affected code/config/docs/tests; do not rely on timestamps or
   stored fingerprints.
5. Reconcile owning leaves and one-hop neighbors. Mark `stale` only when work
   cannot finish immediately.
6. Update dates, catalog metadata, source-map check/action, and ingest state in
   its sole owner.
7. Run lint and scoped searches for old symbols/claims.

## Contradictions and Gaps

When registered authorities disagree, add a fixed 13-column row to
[`contradictions.md`](contradictions.md), add its evidence block, link every
affected page, and add exact page markers before presenting adjudication.
Current code/config decides current-tree behavior; this does not establish
external policy or runtime success.

Keep absent evidence as a gap. Do not invent contradictions from uncertain TOV
test naming, commented Balsara4 coverage, dormant equatorial code, missing
explicit Tabulated selectors, unknown entropy defaults, or unexecuted fixtures.

Before closure, execute the row's resolution test under explicit authority,
review reverse dependents and all pages citing either source, search aliases and
key phrases, reconcile affected page/catalog status, then remove markers and
record `resolved` plus date together.

## Page Add, Move, or Delete

- Add: confirm real fan-out and owner; register sources; create contracted
  leaf/router; update parent, neighbor links, catalog, source map, and glossary.
- Move: preserve content and ownership; update all inbound/outbound links,
  catalog route, source dependents, glossary owner, and contradiction targets.
- Delete: first search inbound links and unique claims; transfer retained facts;
  then remove all routing, catalog, source-map, glossary, and contradiction
  references in the same change.

Routers remain navigation-only. Never add padding routers or bypass branch
routers from root shortcuts.

## Safe Verification

Initial construction permits source inspection, read-only manifest expansion,
Markdown checks, and `python tools/kb_lint.py`. It does not permit Cactus builds,
Cactus test execution, external dependency execution, or rewriting `.asc` or
per-case `.par` fixtures.

A future runtime reproduction needs explicit authority, configured isolated
environment, bounded resources, and recorded revision, command, inputs,
platform, result, and validation scope. Preserve fixture files. A schedule CCL
entry proves declaration, not execution; checked-in output proves only its
stored fixture context, not current pass state.

## Coordination Exemptions

`plan1.md`, `plan2.md`, `plan3.md`, `plan_synth.md`, and commissioned
`tasks1.md` through `tasks4.md` stay outside root routes, catalog, source map,
and checker inventory. No other planning filename is exempt by pattern.
