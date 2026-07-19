# GRHayLID KB Schema

> Page status: reviewed · Last reviewed: 07-19-2026

## Scope

This contract governs only `wiki/grhaylid/**`, `raw/grhaylid/SOURCES.md`, and
the additive GRHayLID profile in `tools/kb_lint.py`. Domain claims may use only
evidence below `GRHayLID/**`. External framework documentation may explain
generic semantics, but cannot establish local build, execution, test, numeric,
library, or provenance claims.

Navigation stays within `wiki/grhaylid/**` and `GRHayLID/**`, except for the
branch index's `AGENTS.md` Up link and links to `raw/grhaylid/SOURCES.md`,
`tools/kb_lint.py`, or `wiki/lint/CHECKS.md`. Local declaration prose may name
another thorn, but sibling source trees are never domain evidence.

## State Vocabularies

- Page status: `draft`, `reviewed`, or `router`.
- Claim status: `declared`, `visible-implementation`, `unresolved`,
  `coverage-gap`, or `out-of-scope`.
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
`GRHayLID/**` source locator or
`kb:wiki/grhaylid/path#section=Heading` pointing to a checked local
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
3. `Mode Applicability`
4. `Claim-Evidence`
5. `Details`
6. `Caveats`
7. `Sources`
8. `Related Pages`

`Claim-Evidence` uses columns `Claim ID | Claim | Status | Evidence | Typed
locator`. Mode applicability uses first-column header `Applicability` and only
`Common`, `HydroTest1D`, `HydroTest1D+Magnetic`, `IsotropicGas`,
`ConstantDensitySphere`, `BetaEquilibrium`, `Entropy/Hybrid`, or
`Entropy/Tabulated`. All eight sections are nonempty; Claim IDs and fields are
nonempty, and Claim IDs are unique within each leaf.

Router preambles contain exactly one nonempty ATX H1, use page status
`router`, review date `n/a`, an `Up` link, and a `Page | Use it for` routing
table. Governance pages may organize their contract-specific material freely,
but must expose reviewed metadata and a catalog row.

## Claim Authority

| Claim kind | Local authority | Required wording | Claim status |
| --- | --- | --- | --- |
| Stated purpose or intent | README or ThornGuide | documentation states/describes | declared |
| Cactus interface, parameter, schedule intent | named CCL declaration | declares/permits/schedules | declared |
| Checked-in compilation surface | `make.code.defn` | lists/includes source unit | declared |
| Visible dataflow, expression, or call order | C/header symbol | reads/writes/calls/computes expression | visible-implementation |
| Shipped configuration, test declaration, or numeric oracle | No such local artifact exists | absence claims cite existing CCL/build/C locators | coverage-gap |
| External behavior | no local authority | delegated/unverified externally | out-of-scope |

Multiple sources may support one claim, but status follows wording actually
asserted. Conflicting local sources require narrow wording plus an issue; no
global "code wins" rule applies.

Canonical edge claim kinds are `stated-purpose`, `cactus-interface`,
`parameter`, `schedule-intent`, `build-surface`, `visible-dataflow`,
`visible-formula`, `visible-call-order`, and `external-behavior`.

The first checked-in parfile, test, or oracle below `GRHayLID/**` requires a
schema-and-lint change that reinstates `checked-in-observation`,
`shipped-configuration`, `test-declaration`, `numeric-observation`, and the
`par:`, `test:`, and `oracle:` locator types before that artifact is cited.

## Typed Locators

Every locator is `TYPE:GRHayLID/path#QUALIFIER`, with optional uniqueness
query. Admitted forms are:

- `doc:path#section=Heading` with optional `?occurrence=N`;
- `ccl:path#parameter=Name`, with optional `?occurrence=N`;
- `ccl:path#schedule=Name`, with optional
  `?context=<schedule-bin-or-group>` or `?occurrence=N`;
- configuration/interface anchors may use `#implementation=Name` or
  `#requirement=Name` with optional occurrence;
- `build:path#field=SRCS`, with optional `?occurrence=N`;
- `c:path#symbol=Name`, or `#call=Name?function=ContainingFunction`;
- `macro:path#name=Name` or `#include=Header`.

Paths must resolve to registered files within `GRHayLID/**`. A qualifier must
resolve exactly once; repeated names need `context` or one-based `occurrence`.
Query keys are qualifier-specific; irrelevant keys and combined CCL
`context`/`occurrence` are invalid. The `ccl:#parameter=` matcher accepts
plain declarations plus `USES KEYWORD <name>` and `EXTENDS KEYWORD <name>`
forms.

Heading matching extracts Markdown ATX/setext, numbered README headings, and
TeX section-family arguments, then applies Unicode NFKC, trims outer
whitespace, collapses internal ASCII whitespace, preserves case, and demands
one match. Generated CCTK macro existence and semantic truth remain manual
review.

`par:`, `test:`, `oracle:`, `ccl:#group=`, and `ccl:#storage=` are not
admitted because the checked-in source tree has no possible target for them.

## Registry and Edge Contract

`raw/grhaylid/SOURCES.md` owns source identity, one exact path or one
non-overlapping pattern, provenance class, lifecycle, ingest state, and family.
Allowed provenance classes are `narrative`, `declaration`, `build input`, and
`implementation`. Pattern expansion is repository-root-relative, sorted,
file-only, nonescaping, and must assign every `rg --files GRHayLID` file to
exactly one record. Families group records but add no coverage.

`wiki/grhaylid/catalog.md` owns page IDs and aliases. Page IDs are mechanical:
prefix `grhaylid.`, lowercase the namespace-relative path, remove `.md`, and
replace `/` with `.`. Specials are `grhaylid.index`, `grhaylid.schema`, and
`grhaylid.lint`.

`wiki/grhaylid/source-map.md` owns the sole edge table:

`Source ID | Page ID | Claim kind | Typed locator | Claim status | Next action`

Every edge uses known IDs, a locator owned by its source record, and an allowed
claim status. Reviewed leaves need at least one edge. Sources marked `ingested`
need at least one consumer. Forward and reverse views are derived, never copied
into separate hand-written maps.

## Proof Boundary

Static inspection proves declarations and visible dataflow only. It does not
prove successful build, scheduled execution, current test pass, numerical
validity, convergence, thread safety, external-library behavior, or provenance
beyond checked-in assertions. It also does not prove HydroBase keyword-
extension acceptance by the framework, `CCTK_VERROR` or `CCTK_ERROR`
termination semantics, EOS-table content, HDF5 discovery, or GRHayLib results.
