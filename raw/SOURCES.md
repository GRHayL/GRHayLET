# IllinoisGRMHD Source Manifest

> Source identity and ingest-state registry. Last audited: 07-17-2026.

## Pattern Grammar

Each row owns one stable source ID and exactly one backtick-delimited,
repo-relative POSIX path or glob. A checker resolves repository root as
`Path(__file__).resolve().parents[1]`, never from process working directory.
Expand with `pathlib.Path.glob`, sort by normalized repo-relative path, and
retain files only. `*` matches one path component. `**` is recursive only as a
complete component. A matched directory does not imply its descendants.

Reject empty values or cells, backslashes, absolute paths, URIs, leading
slashes, empty/`.`/`..` components, paths outside `IllinoisGRMHD/`, duplicate
IDs, duplicate provenance values, malformed table rows, missing or extra
backtick spans, invalid globs, unmatched patterns, directory-only patterns,
and results resolving outside the repository or IllinoisGRMHD tree. Aggregate
and exact rows may overlap intentionally; coverage is their file-set union.
Compare that union with `git ls-files -- IllinoisGRMHD`: all tracked files must
be covered, and untracked matches do not become registered automatically.

`Status` is `living` or `frozen`. `Ingest` is `registered`, `partial`, or
`ingested`. Authority and reverse dependencies belong in
[`wiki/source-map.md`](../wiki/source-map.md), not this manifest.

## Registered Sources

| Source | Provenance | Status | Ingest |
| --- | --- | --- | --- |
| `illinoisgrmhd-readme` | `IllinoisGRMHD/README` | living | ingested |
| `illinoisgrmhd-docs` | `IllinoisGRMHD/doc/**/*` | living | ingested |
| `illinoisgrmhd-root-ccl` | `IllinoisGRMHD/*.ccl` | living | ingested |
| `illinoisgrmhd-common-runtime` | `IllinoisGRMHD/src/*` | living | ingested |
| `illinoisgrmhd-hybrid` | `IllinoisGRMHD/src/Hybrid/**/*` | living | ingested |
| `illinoisgrmhd-hybrid-entropy` | `IllinoisGRMHD/src/HybridEntropy/**/*` | living | ingested |
| `illinoisgrmhd-tabulated` | `IllinoisGRMHD/src/Tabulated/**/*` | living | ingested |
| `illinoisgrmhd-tabulated-entropy` | `IllinoisGRMHD/src/TabulatedEntropy/**/*` | living | ingested |
| `illinoisgrmhd-example-parfiles` | `IllinoisGRMHD/par/**/*` | living | ingested |
| `illinoisgrmhd-test-ccl` | `IllinoisGRMHD/test/test.ccl` | living | ingested |
| `illinoisgrmhd-test-parfiles` | `IllinoisGRMHD/test/**/*.par` | living | ingested |
| `illinoisgrmhd-test-oracles` | `IllinoisGRMHD/test/**/*.asc` | living | partial |
| `illinoisgrmhd-thornguide` | `IllinoisGRMHD/doc/documentation.tex` | living | ingested |
| `illinoisgrmhd-configuration-ccl` | `IllinoisGRMHD/configuration.ccl` | living | ingested |
| `illinoisgrmhd-interface-ccl` | `IllinoisGRMHD/interface.ccl` | living | ingested |
| `illinoisgrmhd-param-ccl` | `IllinoisGRMHD/param.ccl` | living | ingested |
| `illinoisgrmhd-schedule-ccl` | `IllinoisGRMHD/schedule.ccl` | living | ingested |
| `illinoisgrmhd-common-build` | `IllinoisGRMHD/src/make.code.defn` | living | ingested |
| `illinoisgrmhd-hybrid-build` | `IllinoisGRMHD/src/Hybrid/make.code.defn` | living | ingested |
| `illinoisgrmhd-hybrid-entropy-build` | `IllinoisGRMHD/src/HybridEntropy/make.code.defn` | living | ingested |
| `illinoisgrmhd-tabulated-build` | `IllinoisGRMHD/src/Tabulated/make.code.defn` | living | ingested |
| `illinoisgrmhd-tabulated-entropy-build` | `IllinoisGRMHD/src/TabulatedEntropy/make.code.defn` | living | ingested |
| `illinoisgrmhd-header` | `IllinoisGRMHD/src/IllinoisGRMHD.h` | living | ingested |
| `illinoisgrmhd-hybrid-con2prim` | `IllinoisGRMHD/src/Hybrid/conservs_to_prims.c` | living | ingested |
| `illinoisgrmhd-hybrid-entropy-con2prim` | `IllinoisGRMHD/src/HybridEntropy/conservs_to_prims.c` | living | ingested |
| `illinoisgrmhd-tabulated-con2prim` | `IllinoisGRMHD/src/Tabulated/conservs_to_prims.c` | living | ingested |
| `illinoisgrmhd-tabulated-entropy-con2prim` | `IllinoisGRMHD/src/TabulatedEntropy/conservs_to_prims.c` | living | ingested |
| `illinoisgrmhd-a-flux-rhs` | `IllinoisGRMHD/src/A_flux_rhs.c` | living | ingested |
| `illinoisgrmhd-magnetic-curl` | `IllinoisGRMHD/src/compute_B_and_Bstagger_from_A.c` | living | ingested |
| `illinoisgrmhd-gauge-rhs` | `IllinoisGRMHD/src/evaluate_phitilde_and_A_gauge_rhs.c` | living | ingested |
| `illinoisgrmhd-hydrobase-ingress` | `IllinoisGRMHD/src/convert_HydroBase_to_IllinoisGRMHD.c` | living | ingested |
| `illinoisgrmhd-hydrobase-egress` | `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c` | living | ingested |
| `illinoisgrmhd-tmunu` | `IllinoisGRMHD/src/compute_Tmunu.c` | living | ingested |
| `illinoisgrmhd-compat-data` | `IllinoisGRMHD/src/backward_compatible_data.c` | living | ingested |
| `illinoisgrmhd-compat-initialize` | `IllinoisGRMHD/src/backward_compatible_initialize.c` | living | ingested |
| `illinoisgrmhd-example-balsara1` | `IllinoisGRMHD/par/Balsara1.par` | living | ingested |
| `illinoisgrmhd-example-balsara2` | `IllinoisGRMHD/par/Balsara2.par` | living | ingested |
| `illinoisgrmhd-example-balsara3` | `IllinoisGRMHD/par/Balsara3.par` | living | ingested |
| `illinoisgrmhd-example-balsara4` | `IllinoisGRMHD/par/Balsara4.par` | living | ingested |
| `illinoisgrmhd-example-balsara5` | `IllinoisGRMHD/par/Balsara5.par` | living | ingested |
| `illinoisgrmhd-example-magnetized-tov` | `IllinoisGRMHD/par/magnetizedTOV.par` | living | ingested |
| `illinoisgrmhd-test-balsara1` | `IllinoisGRMHD/test/Balsara1.par` | living | ingested |
| `illinoisgrmhd-test-balsara2` | `IllinoisGRMHD/test/Balsara2.par` | living | ingested |
| `illinoisgrmhd-test-balsara3` | `IllinoisGRMHD/test/Balsara3.par` | living | ingested |
| `illinoisgrmhd-test-balsara4` | `IllinoisGRMHD/test/Balsara4.par` | living | ingested |
| `illinoisgrmhd-test-balsara5` | `IllinoisGRMHD/test/Balsara5.par` | living | ingested |
| `illinoisgrmhd-test-magnetized-tov` | `IllinoisGRMHD/test/magnetizedTOV.par` | living | ingested |
| `illinoisgrmhd-fixture-balsara1` | `IllinoisGRMHD/test/Balsara1/Balsara1.par` | living | ingested |
| `illinoisgrmhd-fixture-balsara2` | `IllinoisGRMHD/test/Balsara2/Balsara2.par` | living | ingested |
| `illinoisgrmhd-fixture-balsara3` | `IllinoisGRMHD/test/Balsara3/Balsara3.par` | living | ingested |
| `illinoisgrmhd-fixture-balsara5` | `IllinoisGRMHD/test/Balsara5/Balsara5.par` | living | ingested |
| `illinoisgrmhd-fixture-magnetized-tov` | `IllinoisGRMHD/test/magnetizedTOV/magnetizedTOV.par` | living | ingested |
