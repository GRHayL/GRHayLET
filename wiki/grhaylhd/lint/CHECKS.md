# GRHayLHD KB Lint Checks

> Page status: reviewed · Last reviewed: 07-17-2026

## Canonical Command

Run `python tools/kb_lint.py`. `python tools/kb_lint.py --all` is an identical
compatibility alias. Both run Illinois checks plus the isolated GRHayLHD
profile, sort GRHayLHD diagnostics by check ID and path, and never write to
`GRHayLHD/**`.

## Automated Checks

- Exact namespaced page schema and preamble, dates, statuses, catalog equality, router
  shape and `Up` targets, glossary ownership, issue backlinks, and publication
  gate with an exact root index route. Issue backlinks and affected-page rows
  must agree in both directions; resolved issues require an admissible source
  or namespaced documentation-disposition locator.
- Leaf sections must be nonempty; Claim IDs/fields must be valid and unique;
  variant applicability must use only schema vocabulary.
- Registry syntax, file-only sorted expansion, exact `rg --files GRHayLHD`
  equality, no unmatched/escaping/overlapping ownership, ingested-source
  consumers, and consumers for every source once target is complete.
- Typed-locator syntax and resolution, registered ownership, unique targets,
  exact schedule-container contexts, qualifier-specific query keys,
  source/page edge integrity, exact claim-locator/status/page pairs, direct
  evidence-link/page pairs, and
  reviewed-page evidence.
- Namespaced navigation isolation, no Magnetics branch, no root GRHayLHD link
  while unpublished, no forbidden metadata, and exact coordination exemptions.
- Existing Illinois profile remains enabled with unchanged meaning.

## Manual Boundary

Checker does not decide scientific correctness, generated CCTK macro
existence, actual schedule execution, external semantics, numerical order not
locally stated, oracle validity or provenance, current build/test status,
thread safety, convergence, conservation, restart behavior, or performance.

## Negative Fixtures

Use only `mktemp -d` copies. Required fixtures cover malformed schema state,
source-prefix escape, overlapping registry patterns, unresolved or non-unique
locators, missing issue backlink, and premature root routing while publication
is `unpublished`. Compare default and `--all` exit code, stdout, stderr, and
emitted check IDs directly; never store fingerprints or modification times.

## Related Pages

- [Schema](../SCHEMA.md)
- [Catalog](../catalog.md)
- [Source Map](../source-map.md)
