# GRHayLID KB Lint Checks

> Page status: reviewed · Last reviewed: 07-19-2026

## Canonical Command

Run `python tools/kb_lint.py`. `python tools/kb_lint.py --all` is an identical
compatibility alias. Both run Illinois and GRHayLHD checks plus the isolated
GRHayLID profile, sort GRHayLID diagnostics by check ID and path, and never
write to `GRHayLID/**`.

Checked-in KB-linter tests, regression harnesses, and fixture suites are
forbidden regardless of filename or location. Use disposable manual fixtures
outside the repository when a checker regression needs reproduction.

## Automated Checks

- Exact 24-page target, namespaced page schema and preamble, dates, statuses,
  catalog equality, router shape and `Up` targets, glossary ownership, issue
  backlinks, and publication gate with an exact root route.
- Leaf sections must be nonempty; Claim IDs and fields must be valid and
  unique; mode applicability must use only schema vocabulary.
- Registry syntax, file-only sorted expansion, exact `rg --files GRHayLID`
  equality, no unmatched, escaping, or overlapping ownership, ingested-source
  consumers, and consumers for every source once the target is complete.
- Five-type locator syntax and resolution, registered ownership, unique
  targets, qualifier-specific query keys, source/page edge integrity, exact
  claim-locator/status/page pairs, direct evidence-link/page pairs, and
  reviewed-page evidence.
- Issue table shape, IDs, kinds, statuses, locators, resolution locators,
  anchors, and bidirectional affected-page backlinks.
- Namespaced navigation isolation; slash-anchored sibling-tree mention scans;
  no root GRHayLID link while unpublished; no forbidden metadata; and exact
  coordination exemptions.
- Existing Illinois and GRHayLHD profiles remain enabled with unchanged
  meaning.

## Manual Boundary

Checker does not decide scientific correctness, generated CCTK macro
existence, actual schedule execution, external semantics, error-macro
termination, EOS-table content, HDF5 discovery, numerical validity, current
build/test status, thread safety, convergence, conservation, restart behavior,
or performance.

## Negative Fixtures

Manual disposable fixtures should cover default/`--all` equivalence;
commented and string-literal C symbol decoys; missing, misplaced, and duplicate
ready root routes; Branch Selector multiplicity; wrong leaf `Up`; mixed
file/directory registry globs; `EXTENDS KEYWORD` resolution; duplicate,
malformed, and dangling domain-router rows; missing, empty, and multiple router
H1s; issue backlinks and anchors; and GRHayLHD locator behavior. Use only
`mktemp -d` copies outside the repository. Compare exit code, stdout, stderr,
and check IDs directly; never store fingerprints or modification times.

## Related Pages

- [Schema](../SCHEMA.md)
- [Catalog](../catalog.md)
- [Source Map](../source-map.md)
- [Repository Checker](../../../tools/kb_lint.py)
