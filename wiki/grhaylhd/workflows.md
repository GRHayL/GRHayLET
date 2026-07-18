# GRHayLHD KB Workflows

> Page status: reviewed · Last reviewed: 07-17-2026

## Register Sources

1. Run `rg --files GRHayLHD` from repository root and sort results.
2. Add one exact path or one non-overlapping pattern to
   [source registry](../../raw/grhaylhd/SOURCES.md).
3. Confirm every expansion is file-only, remains below `GRHayLHD/**`, and has
   exactly one owner. Do not store fingerprints, digests, or modification
   times.
4. Keep new records `registered` until at least one reviewed page consumes a
   resolving typed locator through [source map](source-map.md).

## Add or Reconcile Pages

1. Start from [schema](SCHEMA.md); use `draft` during construction.
2. Inspect only local GRHayLHD evidence needed for each claim. Phrase CCL and
   build facts as declarations, C facts as visible implementation, parfiles and
   oracles as checked-in observations, and external semantics as out of scope.
3. Add page ID to [catalog](catalog.md) and canonical edge rows to
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
mutate or execute `GRHayLHD/**`, regenerate oracles, or run Cactus tests during
KB maintenance.

## Publication

Keep branch `unpublished` until exact target, routes, reviewed leaves, catalog,
edges, issue backlinks, registry, and deterministic checks close. Publication
is one atomic slice: set branch index to `ready` and add root route together.
Any later structural drift returns publication to `unpublished` before root
navigation changes.

Before publication, probe these routes from [branch index](index.md):

- purpose/build through [Architecture](architecture/index.md);
- mode selection and C2P failure through [Evolution](evolution/index.md);
- velocity conversion through [Integration](integration/index.md);
- matter boundaries through [Evolution](evolution/index.md); and
- tests/oracles through [Validation](validation/index.md).

Each domain leaf must resolve within two links from branch index. Navigation
stays in `wiki/grhaylhd/**` except branch index's exact root `Up` link.
