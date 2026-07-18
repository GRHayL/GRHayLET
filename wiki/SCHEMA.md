# IllinoisGRMHD Knowledge Base Schema

> Governance contract. · Status: confirmed · Last reconciled: 07-17-2026

## Summary

This KB has three ownership layers: registered IllinoisGRMHD sources, compiled
wiki pages, and governance pages. Navigation starts at `AGENTS.md`; routers
only navigate, while sourced leaves own durable implementation facts.

## Scope and Ownership

- Scope is `/work/IllinoisGRMHD/**` plus KB files named by the root inventory.
  Other source trees are neither evidence nor routing targets.
- [`../raw/SOURCES.md`](../raw/SOURCES.md) solely owns source identity,
  provenance, `living`/`frozen` status, and
  `registered`/`partial`/`ingested` state.
- [`source-map.md`](source-map.md) solely owns authority tiers, reverse page
  dependencies, dependency coverage, gaps, checks, and next actions. Its
  dependency coverage states are `unmapped`, `partial`, and `complete`.
- Compiled leaves own IllinoisGRMHD explanations. Governance owns contracts and
  procedures. Routers own no durable implementation fact.
- `plan1.md`, `plan2.md`, `plan3.md`, `plan_synth.md`, `tasks1.md`, `tasks2.md`,
  `tasks3.md`, and `tasks4.md` are coordination artifacts explicitly exempt
  from root routing, catalog inventory, source registration, and KB checks.
  No other plan or task filename is implicitly exempt.

## Authority

Authority is role-sensitive:

- `primary-code`: C, headers, CCL, build inputs, and checked-in runtime/example
  configuration. These decide current locally visible behavior and surfaces.
- `primary-test`: `test/test.ccl` and test parfiles. They decide configured or
  declared fixture conditions only.
- `generated-evidence`: checked-in `.asc` outputs. They are evidence for their
  checked-in fixture context, not proof of a current pass.
- `primary-doc`: README and ThornGuide. They decide documented intent,
  attribution, and migration guidance unless code/config contradict them.
- `background`: orientation only; it never decides an IllinoisGRMHD claim.

For current behavior, local code/config wins over tests, generated evidence,
and prose. Intended public or scientific contracts require their stable owning
declaration and targeted test support where available. Synthesis agreement is
never independent authority. External web material may appear only as clearly
labeled, non-deciding background. Every IllinoisGRMHD claim must still cite
registered local evidence, and external material cannot override it.

## Source Tracking and Dates

All retained dates use `MM-DD-YYYY`; routers use `n/a` where no reconciliation
date applies. Do not compute or store source fingerprints, hash values, digest
columns, or modification-time fields. Reconcile living-source drift by path,
stable locator, reverse dependencies, and affected claims.

Every current IllinoisGRMHD file starts `living` and `registered`, including
checked-in `.asc` files. Generated-evidence authority is distinct from source
lifecycle. Promote ingest only after dependent leaves have been reconciled.

## Leaf Contract

Leaf order is exact:

1. one H1 title;
2. exactly two blockquote lines, first with `Status`, an allowed status,
   `Last reconciled`, and an `MM-DD-YYYY` value; second with `Up:` and one
   parent-router link;
3. `## Summary`;
4. `## Detail`;
5. `## Sources`;
6. `## See Also`.

Leaf status is `confirmed`, `provisional`, `contested`, or `stale`. Use
`confirmed` when registered sources directly support the principal answer;
`provisional` when useful structure lacks full evidence; `contested` when an
active conflict changes a central answer; and `stale` when living-source drift
makes central content unreliable. A bounded, directly sourced contradiction
may live on a confirmed page if its marker and non-guarantee are explicit.
Catalog status/date must match the header.

Each leaf cites at least one registered source using a stable locator. Code
locators are functions, symbols, macros, or header roles; CCL locators are
group, parameter, function-alias, schedule, or storage names; build locators
are `SRCS`/`SUBDIRS` roles; docs use headings; parfiles/oracles use fixture
roles. Generated line numbers are not authorities. `See Also` contains a
`Parent:` link plus at least one typed neighbor: `Depends on:`, `Implements:`,
`Example:`, `Contrasts with:`, or `See also:`. Reserve `Validated by:` for an
executed comparison recording command, environment, revision, measured result,
and scope.

KB routing and citations use ordinary single-line inline Markdown links, the
same authoring and validation boundary as the NRPy template. Images,
reference-style links, and multiline or nested CommonMark link forms do not
serve as KB edges or citations; the checker is not a general Markdown parser.

## Router Contract

Routers contain one H1, one-line blockquote, a
`Page | Go here when...` table, immediate children, parent/root navigation,
and global support links. Their catalog type/status are `router`; date is
`n/a`. They have no `Detail`, source dump, ordinary `Sources`, or unique
implementation claim. Routers link only to immediate children plus parent,
root, and global navigation. Root task shortcuts enter a branch router, never
bypass it to a deep leaf. Add depth only for genuine fan-out; do not create
single-child padding.

## Page Status Decision Matrix

- `contested`: conflict changes Summary, normal routed answer, interface,
  command, guarantee, or multiple required Detail claims.
- `stale`: unreconciled living-source drift undermines principal content.
- `provisional`: evidence/reconciliation is incomplete without known conflict.
- `confirmed`: principal answer remains directly supported; any active defect
  is bounded, marked, and accompanied by an explicit non-guarantee.

Claim status and page status differ. Catalog metadata follows page status, not
claim status.

## Claim-Evidence Contract

Add this block immediately after every new or materially changed high-risk
claim (public/scientific contract, user command/interface, generated-evidence
boundary, authority/contradiction decision, CI guarantee, or easily misread
claim). Active contradiction blocks live in their `CONTR-*` subsection, not
the fixed register row.

```text
Claim evidence:
- Claim: exact qualified text, including conditions and non-guarantees
- Role: descriptive behavior, normative rule, public/scientific contract, CI behavior, or generated evidence
- Deciding authority: registered source plus stable symbol or heading
- Corroboration: separate source plus stable locator, or `none available` plus a reason
```

Behavioral claims also add:

```text
- Validation: `inspected=<pass|fail|not-run>; generated=<pass|fail|not-run>; built=<pass|fail|not-run>; run=<pass|fail|not-run>; result_checked=<pass|fail|not-run>`
- Dimensions: `platform=<value|not-run|not-applicable>; tool_version=<value|not-run|not-applicable>; backend=<value|not-run|not-applicable>; precision=<value|not-run|not-applicable>; GPU=<value|not-run|not-applicable>; restart=<value|not-run|not-applicable>; distributed=<value|not-run|not-applicable>; error_path=<value|not-run|not-applicable>; options=<value|not-run|not-applicable>; date=<MM-DD-YYYY|not-run|not-applicable>`
```

Static inspection is not a build, run, observed schedule execution, numerical
validation, or external-library guarantee.

## Contradiction Contract

[`contradictions.md`](contradictions.md) uses exactly:

`ID | Claim | Claim status | Source A | Source B | Authority decision | Affected pages | Page-status rationale | Owner/trigger | Resolution test | Opened | Resolved | Notes`

IDs are immutable `CONTR-0001` forms. Active statuses are `contested` or
`stale`; closed rows use `resolved` and a resolution date. Every active affected
page carries `Claim status: <status>; contradiction: CONTR-####.` with a link
back to its register heading. Open a row before treating disagreement as a
resolved fact. Resolution requires its test, reverse-dependent review, page and
catalog reconciliation, neighbor/alias search, and marker removal together.
Unknown behavior, missing coverage, unsupported configuration, or absent
external semantics is a gap until two authoritative claims actually disagree.

## Catalog, Glossary, and Source Map

Catalog lists every `wiki/**/*.md` page once using exact columns
`Page | Type | One-line answer | Route | Tags | Aliases / Query terms | Status |
Last reconciled`. Glossary defines recurring local routing terms, with exactly
one valid owner page each. Source map includes every manifest ID exactly once
and uses exact columns `Source / aggregate | Authority tier | Dependency
coverage | Dependent pages | Covered paths / stable locators | Known gaps |
Last check | Next action`. It never stores ingest vocabulary.

## Query Filing and Page Changes

File durable single-topic answers in their existing owner. Create a new leaf
only for a durable narrow topic with a clear branch owner. Create a synthesis
branch only when a recurring cross-branch answer has no clean existing owner.
One-off answers need no page.

Adding, moving, or deleting a page updates its router, typed neighbors,
catalog, source-map dependencies, glossary owners, and contradiction affected
pages in one change. Before deletion, search inbound links and transfer unique
sourced content. Normal questions should resolve root, branch router, owner
leaf; broad source grep signals routing or coverage drift.

## Scope and Safety Guard

Initial KB work may inspect IllinoisGRMHD sources and run deterministic KB
lint only. It must not build/run Cactus, execute the test suite, regenerate or
rewrite `.asc`/`.par` fixtures, mutate IllinoisGRMHD, or cite sibling trees.
Runtime reproduction requires explicit authority, a configured isolated
environment, resource limits, and a recorded command/context/result. Never
upgrade declared, scheduled, configured, or checked-in evidence to observed
execution.
