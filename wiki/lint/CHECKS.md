# IllinoisGRMHD KB Lint Checks

> Deterministic checks and manual review boundary. · Status: confirmed · Last reconciled: 07-17-2026

## Canonical Command

Run `python tools/kb_lint.py`. `python tools/kb_lint.py --all` invokes the
same dispatcher, checks, diagnostic ordering, output, and exit code. Success
prints exactly `KB lint passed.`. Checker resolves repository root from its own
file path and never writes to `IllinoisGRMHD/`.

## Hard Failures

- Ordinary single-line inline Markdown links use the same validation boundary
  as NRPy's checker: `(?<!!)\[([^\]\n]+)\]\(([^)\n]+)\)`. Recognized local
  links must stay within repository, resolve to existing files, and name
  existing Markdown anchors. Authors use this form for KB routing and
  citations; images, reference-style links, and multiline/nested CommonMark
  forms are outside deterministic validation. Obsidian links are forbidden.
- Exact 31-node inventory must exist: root, 28 wiki pages, source manifest,
  and checker. Missing required nodes produce stable diagnostics without later
  dereference. All wiki pages must be reachable from root. Routers must have
  exact immediate-child fan-out, no `Detail`/`Sources`, and no illegal edges.
- Every leaf must begin with H1, blank line, exact status/date line, exact
  `Up:` line, blank line, and immediate `Summary`; its H2 order is `Summary`,
  `Detail`, `Sources`, `See Also`. Sections must be nonempty, parent must be
  correct, and at least one typed neighbor must resolve to another live wiki
  page.
- Catalog live-page set, uniqueness, type, status, and date must agree with
  pages. Its Markdown table requires a delimiter row. Router rows use `router`
  and `n/a`; recovery remains `contested` while `CONTR-0002` is active.
- Manifest table must have exact four columns and valid one-literal rows.
  Its Markdown table requires a delimiter row.
  IDs/provenance are unique; paths/globs are normalized IllinoisGRMHD-relative,
  nonescaping, syntactically valid, matched, and file-producing. Union must
  cover every `git ls-files -- IllinoisGRMHD` path. Directory matches do not
  imply descendants.
- Every local leaf link and existing path-like inline literal outside wiki/
  governance navigation must lie under `IllinoisGRMHD/**`. Every local source
  path/glob in `Sources` must be covered by a manifest provenance pattern.
  Every named manifest ID anywhere in leaf prose must exist. Source paths and
  IDs are recognized with or without inline-code formatting.
- Source map must use exact columns and allowed authority/dependency-coverage
  values, contain every manifest ID exactly once with no unknown IDs, link only
  live dependent pages, and give each `complete` row every direct dependent
  derived from registered IDs and IllinoisGRMHD path/glob citations, whether
  bare or inline-code formatted. It must not contain an ingest column or
  source-state vocabulary; its table requires a delimiter row.
- Contradiction table requires a delimiter row, exact columns, unique IDs,
  statuses/dates, linked affected pages, reciprocal exact markers/statuses/
  backlinks, and matching active affected-page sets.
- Glossary table requires a delimiter row; every row has exactly one live wiki
  owner page and no duplicate term.
- Governed text (`AGENTS.md`, `wiki/**/*.md`, `raw/SOURCES.md`, and
  `tools/kb_lint.py`) requires final newline and forbids trailing horizontal
  whitespace, Obsidian links, unfilled template markers, stored source
  fingerprints/modification-time values, every empty or invalid retained date,
  generic formatted or plain `Hash:`/`Digest:` values (including `0x` hex),
  and non-Markdown files under `wiki/` or `raw/`.
- Exact commissioned exemptions are `plan1.md`, `plan2.md`, `plan3.md`,
  `plan_synth.md`, and `tasks1.md` through `tasks4.md`. Other root plan/task
  maintenance artifacts—including natural or dotted names such as `plan4.md`,
  `plan.4.md`, `plan.notes.md`, and `tasks.5.md`—or any commissioned artifact
  routed/cataloged as KB content, fail.

## Manual and Report-Only Review

Checker does not decide scientific truth, generated CCTK macro/symbol
existence, actual schedule execution, external-library semantics, numerical
oracle validity, current build/test pass state, convergence, divergence, AMR,
restart, or performance. Review Summary-to-source claims, stable locators,
four-variant comparisons, route probes, source-map reverse lookups, and scope
preservation manually. Runtime reproduction needs explicit authority and an
isolated configured environment; structural KB completion does not run Cactus
or regenerate fixtures.

## See Also

- [Schema](../SCHEMA.md)
- [Workflows](../workflows.md)
- [Source Map](../source-map.md)
