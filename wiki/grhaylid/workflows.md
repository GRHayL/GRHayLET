# GRHayLID KB Workflows

> Page status: reviewed · Last reviewed: 07-19-2026

## Register Sources

1. Run `rg --files GRHayLID` from repository root and sort results.
2. Add one exact path or one non-overlapping pattern to the
   [source registry](../../raw/grhaylid/SOURCES.md).
3. Confirm every expansion is file-only, remains below `GRHayLID/**`, and has
   exactly one owner. Do not store fingerprints, digests, or modification
   times.
4. Keep new records `registered` until at least one reviewed page consumes a
   resolving typed locator through the [source map](source-map.md).

## Add or Reconcile Pages

1. Start from [schema](SCHEMA.md); use `draft` during construction.
2. Inspect only local GRHayLID evidence needed for each claim. Phrase CCL and
   build facts as declarations, C/header facts as visible implementation,
   checked-in absences as coverage gaps, and external semantics as out of
   scope.
3. Add the page ID to [catalog](catalog.md) and canonical edge rows to the
   [source map](source-map.md) in the same slice.
4. Review both sides of any disagreement. Add or update a typed record in
   [issues](contradictions.md); every live affected page links its exact anchor.
5. Promote to `reviewed` only after static claim, locator, edge, link, and issue
   review. This never records build, run, or test success.

## Verify Safely

Run both canonical commands:

```text
python tools/kb_lint.py
python tools/kb_lint.py --all
```

Run each twice when checking determinism; compare exit code, stdout, stderr,
and check IDs directly. Use `mktemp -d` copies for negative fixtures. Never
mutate or execute `GRHayLID/**`, generate source artifacts, build the thorn, or
run Cactus tests during KB maintenance.

## Publication

Keep the branch `unpublished` until exact target, routes, reviewed leaves,
catalog, edges, issue backlinks, registry, and deterministic checks close.
Publication is one atomic slice: set branch index to `ready` and add every root
route together. Any later structural drift returns publication to
`unpublished` before root navigation changes.

Before publication, probe these routes from the [branch index](index.md):

- purpose and build through [Architecture](architecture/index.md);
- a Balsara test, magnetic initialization, beta equilibrium, and entropy
  through [Initial Data](initial-data/index.md);
- keyword extensions and `stagger_A_fields` through
  [Integration](integration/index.md); and
- checked-in test absence through [Validation](validation/index.md).

Each domain leaf must resolve within two links from branch index. Navigation
stays in `wiki/grhaylid/**` except the branch index's exact root `Up` link and
the governance allowances in the schema.
