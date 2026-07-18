# GRHayLHD KB Schema

> Page status: reviewed · Last reviewed: 07-17-2026

## Scope

This contract governs only `wiki/grhaylhd/**`, `raw/grhaylhd/SOURCES.md`, and
the additive GRHayLHD profile in `tools/kb_lint.py`. Domain claims may use only
evidence below `GRHayLHD/**`. External framework documentation may explain
generic semantics, but cannot establish local build, execution, test, numeric,
library, or provenance claims.

## State Vocabularies

- Page status: `draft`, `reviewed`, or `router`.
- Claim status: `declared`, `visible-implementation`,
  `checked-in-observation`, `unresolved`, `coverage-gap`, or `out-of-scope`.
- Issue kind: `contradiction`, `mismatch`, `hazard`,
  `lifecycle-ambiguity`, or `provenance-ambiguity`.
- Issue status: `open`, `accepted`, or `resolved`.
- Branch publication: `unpublished` or `ready`.

New leaves start `draft`. A leaf becomes `reviewed` only after section, source
edge, typed-locator, affected-issue backlink, and static claim review passes.
`reviewed` never means built, run, tested, or numerically validated. Routers
always use `router` and review date `n/a`.

Issues start `open`. `accepted` keeps an intentionally unresolved issue
visible. `resolved` requires changed admissible local evidence or an explicit
documentation-only disposition with a resolution locator. Agreement alone
cannot resolve an issue.

Issue table includes `Resolution locator` before `Opened` and `Resolved`.
Active issues use `-`. A resolved issue uses either an admitted typed
`GRHayLHD/**` source locator or
`kb:wiki/grhaylhd/path#section=Heading` pointing to a checked local
documentation-only disposition.

Branch stays `unpublished` until every target page exists, every routed leaf
is `reviewed`, routes/catalog/edges agree, and all checks pass. Structural
drift returns it to `unpublished` before root routing changes.

## Page Contract

Leaf preamble is H1, blank line, exact metadata lines, blank line:

```text
# Title

> Page status: draft|reviewed · Last reviewed: MM-DD-YYYY
> Up: [Router](index.md)
```

Leaf H2 sections occur once in this order and remain nonempty:

1. `Scope and Non-Scope`
2. `Summary`
3. `Variant Applicability`
4. `Claim-Evidence`
5. `Details`
6. `Caveats`
7. `Sources`
8. `Related Pages`

`Claim-Evidence` uses columns `Claim ID | Claim | Status | Evidence | Typed
locator`. Variant applicability uses only `Common`, `Hybrid/Simple`,
`Hybrid/Simple+Entropy`, `Tabulated`, or `Tabulated+Entropy`.
All eight sections are nonempty; Claim IDs and fields are nonempty, and Claim
IDs are unique within each leaf.

Router preambles use page status `router`, review date `n/a`, an `Up` link,
and a `Page | Use it for` routing table. Governance pages may organize their
contract-specific material freely, but must expose reviewed metadata and a
catalog row.

## Claim Authority

| Claim kind | Local authority | Required wording | Claim status |
| --- | --- | --- | --- |
| Stated purpose or intent | README or ThornGuide | documentation states/describes | declared |
| Cactus interface, parameter, schedule intent | named CCL declaration | declares/permits/schedules | declared |
| Checked-in compilation surface | `make.code.defn` | lists/includes source unit | declared |
| Visible dataflow, expression, or call order | C/header symbol | reads/writes/calls/computes expression | visible-implementation |
| Shipped configuration | exact parfile assignment | sets/activates/omits explicit setting | checked-in-observation |
| Regression case or tolerance | `test/test.ccl` | declares case/tolerance | declared |
| Checked-in numeric content | exact ASCII file/header/data | contains checked-in output | checked-in-observation |
| External behavior | no local authority | delegated/unverified externally | out-of-scope |

Multiple sources may support one claim, but status follows wording actually
asserted. Conflicting local sources require narrow wording plus an issue; no
global "code wins" rule applies.

Canonical edge claim kinds are `stated-purpose`, `cactus-interface`,
`parameter`, `schedule-intent`, `build-surface`, `visible-dataflow`,
`visible-formula`, `visible-call-order`, `shipped-configuration`,
`test-declaration`, `numeric-observation`, and `external-behavior`.

## Typed Locators

Every locator is `TYPE:GRHayLHD/path#QUALIFIER`, with optional uniqueness
query. Admitted forms:

- `doc:path#section=Heading` with optional `?occurrence=N`;
- `ccl:path#group=Name` or `#parameter=Name`, with optional `?occurrence=N`;
- `ccl:path#schedule=Name`, with optional
  `?context=<schedule-bin-or-group>` or `?occurrence=N`;
- `ccl:path#storage=Name`, with optional `?occurrence=N`;
  configuration/interface anchors may use
  `#implementation=Name` or `#requirement=Name` with optional occurrence;
- `build:path#field=SRCS`, with optional `?occurrence=N`;
- `c:path#symbol=Name`, or `#call=Name?function=ContainingFunction`;
- `macro:path#name=Name` or `#include=Header`;
- `par:path#parameter=Thorn::name`, with optional `?occurrence=N`, or
  `par:path#file` for a whole-file absence or multi-line header assertion;
- `test:path#case=Name`, with optional `?occurrence=N`;
- `oracle:path#dataset=Name`, or `oracle:path#file`.

Paths must resolve to registered files within `GRHayLHD/**`. A qualifier must
resolve exactly once; repeated names need `context` or one-based `occurrence`.
Query keys are qualifier-specific; irrelevant keys and combined CCL
`context`/`occurrence` are invalid.
Heading matching extracts Markdown ATX/setext, numbered README headings, and
TeX section arguments, then applies Unicode NFKC, trims outer whitespace,
collapses internal ASCII whitespace, preserves case, and demands one match.
Generated CCTK macro existence and semantic truth remain manual review.

## Registry and Edge Contract

`raw/grhaylhd/SOURCES.md` owns source identity, one exact path or one
non-overlapping pattern, provenance class, lifecycle, ingest state, and family.
Pattern expansion is repository-root-relative, sorted, file-only, nonescaping,
and must assign every `rg --files GRHayLHD` file to exactly one record.
Families group records but add no coverage.

`wiki/grhaylhd/catalog.md` owns page IDs and aliases.
Page IDs are mechanical: lowercase the namespace-relative path, remove `.md`,
and replace `/` with `.`. Branch `index.md` uses its branch name, root index
uses `grhaylhd.index`, `SCHEMA.md` uses `grhaylhd.schema`, and
`lint/CHECKS.md` uses `grhaylhd.lint`.
`wiki/grhaylhd/source-map.md` owns the sole edge table:

`Source ID | Page ID | Claim kind | Typed locator | Claim status | Next action`

Every edge uses known IDs, a locator owned by its source record, and an allowed
claim status. Reviewed leaves need at least one edge. Sources marked `ingested`
need at least one consumer. Forward and reverse views are derived, never copied
into separate hand-written maps.

## Proof Boundary

Static inspection proves declarations and visible dataflow only. It does not
prove successful build, scheduled execution, current test pass, numerical
validity, convergence, thread safety, external-library behavior, or oracle
provenance beyond checked-in assertions.
