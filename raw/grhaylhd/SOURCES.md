# GRHayLHD Source Registry

> Page status: reviewed · Last reviewed: 07-17-2026

Each row owns one exact repository path or one non-overlapping pattern.
Expansion is from repository root; Family is grouping metadata only.

| Source ID | Provenance | Provenance class | Lifecycle | Ingest | Family |
| --- | --- | --- | --- | --- | --- |
| `grhaylhd-readme` | `GRHayLHD/README` | narrative | living | registered | documentation |
| `grhaylhd-thornguide` | `GRHayLHD/doc/documentation.tex` | narrative | living | registered | documentation |
| `grhaylhd-ccl` | `GRHayLHD/*.ccl` | declaration | living | registered | Cactus interface |
| `grhaylhd-build-root` | `GRHayLHD/src/make.code.defn` | build input | living | registered | build manifests |
| `grhaylhd-build-variants` | `GRHayLHD/src/*/make.code.defn` | build input | living | registered | build manifests |
| `grhaylhd-common-header` | `GRHayLHD/src/GRHayLHD.h` | implementation | living | registered | common implementation |
| `grhaylhd-common-c` | `GRHayLHD/src/*.c` | implementation | living | registered | common implementation |
| `grhaylhd-variant-c` | `GRHayLHD/src/*/*.c` | implementation | living | registered | four variants |
| `grhaylhd-example-par` | `GRHayLHD/par/*.par` | shipped configuration | living | registered | example input |
| `grhaylhd-test-inputs` | `GRHayLHD/test/*.par` | authored test input | living | registered | tests |
| `grhaylhd-test-companions` | `GRHayLHD/test/*/*.par` | checked-in companion configuration | frozen | registered | tests |
| `grhaylhd-test-declarations` | `GRHayLHD/test/test.ccl` | test declaration | living | registered | tests |
| `grhaylhd-oracles` | `GRHayLHD/test/*/*.asc` | checked-in numeric observation | frozen | registered | tests |

No row records a fingerprint, digest, or modification time. `registered` means
identity and ownership are known; later domain review may promote a record to
`ingested` only after canonical page edges exist.
