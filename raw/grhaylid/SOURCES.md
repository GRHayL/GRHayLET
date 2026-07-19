# GRHayLID Source Registry

> Page status: reviewed · Last reviewed: 07-19-2026

Each row owns one exact repository path or one non-overlapping pattern.
Expansion is from repository root; Family is grouping metadata only.

| Source ID | Provenance | Provenance class | Lifecycle | Ingest | Family |
| --- | --- | --- | --- | --- | --- |
| `grhaylid-readme` | `GRHayLID/README` | narrative | living | ingested | documentation |
| `grhaylid-thornguide` | `GRHayLID/doc/documentation.tex` | narrative | living | ingested | documentation |
| `grhaylid-ccl` | `GRHayLID/*.ccl` | declaration | living | ingested | Cactus interface |
| `grhaylid-build` | `GRHayLID/src/make.code.defn` | build input | living | ingested | build manifest |
| `grhaylid-header` | `GRHayLID/src/GRHayLID.h` | implementation | living | ingested | common implementation |
| `grhaylid-c` | `GRHayLID/src/*.c` | implementation | living | ingested | initial-data implementation |

No row records a fingerprint, digest, or modification time. `registered` means
identity and ownership are known; domain review promotes a record to `ingested`
only after canonical page edges exist.
