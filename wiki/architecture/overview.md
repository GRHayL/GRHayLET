# IllinoisGRMHD Overview

> Status: confirmed · Last reconciled: 07-17-2026
> Up: [Architecture](index.md)

## Summary

IllinoisGRMHD is documented as a Cactus thorn that evolves general-relativistic
magnetohydrodynamics using GRHayL through the inherited GRHayLib thorn. This
repository supplies CCL surfaces, common runtime code, four EOS/entropy source
families, examples, tests, and checked-in oracle outputs; external library
internals remain outside this KB.

## Detail

### Purpose and Provenance

README `Purpose` and ThornGuide `Introduction` describe a high-resolution
shock-capturing GRMHD evolution thorn built on GRHayL through GRHayLib. They
attribute lineage to the original IllinoisGRMHD and Illinois Numerical
Relativity code, list Samuel Cupp, Leonardo Rosa Werneck, Terrence Pierre
Jacques, and Zachariah B. Etienne as current authors, and identify BSD-2 in the
README. README also lists maintainers and required/optional citations.

Claim evidence:
- Claim: IllinoisGRMHD is documented as a GRMHD evolution thorn built on GRHayL via GRHayLib; this KB does not claim GRHayL internals.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/README`, `Purpose`
- Corroboration: registered `IllinoisGRMHD/doc/documentation.tex`, `Introduction`; `IllinoisGRMHD/configuration.ccl`, declared `GRHayL` requirement

### Repository Shape

- Root `configuration.ccl`, `interface.ccl`, `param.ccl`, and `schedule.ccl`
  declare dependencies, variables/functions, controls, and schedule.
- `doc/` and `README` own local documentation and attribution.
- `src/` contains common integration, magnetic, registration, reconstruction,
  compatibility, and coupling code.
- `src/Hybrid`, `HybridEntropy`, `Tabulated`, and `TabulatedEntropy` contain
  parallel conversion, recovery, boundary, perturbation, source, and flux
  implementations.
- `par/` contains six example configurations. `test/` contains the CCL test
  declaration, top-level and checked-in per-case parfiles, and `.asc` outputs.
  File presence alone is not a current test result.

### Responsibility Boundary

Local CCL includes `GRHayLib.h` and inherits `GRHayLib`; local C functions call
`ghl_*` interfaces. This tree decides how IllinoisGRMHD declares, schedules,
assembles, and calls those interfaces. It does not expose or decide GRHayL/
GRHayLib algorithm internals, full supported-method semantics, or transitive
dependency behavior.

### Vector-Potential Intent

README `Purpose` says IllinoisGRMHD evolves staggered vector potential rather
than magnetic fields directly to keep magnetic fields divergence-free across
AMR boundaries, and says uniform-grid results are equivalent to staggered
flux-constrained transport. Treat these as repository-attributed design intent,
not a newly executed divergence/AMR theorem. Exact local placement and curl
construction belong to [Staggered State and Magnetic Reconstruction](../magnetics/staggered-state-and-magnetic-reconstruction.md).

Claim evidence:
- Claim: Repository documentation attributes divergence-control and uniform-grid equivalence intent to staggered vector-potential evolution; initial static ingest does not validate those guarantees.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/README`, `Purpose`
- Corroboration: registered `IllinoisGRMHD/interface.ccl`, `Ax`, `Ay`, `Az`, and `phitilde` declarations establish local staggered state but not numerical guarantees

## Sources

- [`IllinoisGRMHD/README`](../../IllinoisGRMHD/README) — authors, maintainers,
  license, `Purpose`, modeled systems, and citations.
- [`IllinoisGRMHD/doc/documentation.tex`](../../IllinoisGRMHD/doc/documentation.tex)
  — `Introduction`, `Parameters`, `Updating Old Parfiles`, and acknowledgements.
- [`IllinoisGRMHD/configuration.ccl`](../../IllinoisGRMHD/configuration.ccl) —
  declared HDF5 and GRHayL requirements.
- [`IllinoisGRMHD/src/make.code.defn`](../../IllinoisGRMHD/src/make.code.defn)
  — common `SRCS` and four variant `SUBDIRS` repository shape.
- [`IllinoisGRMHD/test/test.ccl`](../../IllinoisGRMHD/test/test.ccl) — local
  test-declaration role; representative checked-in generated-evidence roles are
  [`Balsara1/rho.x.asc`](../../IllinoisGRMHD/test/Balsara1/rho.x.asc) and
  [`magnetizedTOV/rho.maximum.asc`](../../IllinoisGRMHD/test/magnetizedTOV/rho.maximum.asc).

## See Also

- Parent: [Architecture](index.md)
- See also: [Cactus Surface and Build](cactus-surface-and-build.md)
- Depends on: [HydroBase, GRHayLib, and Tmunu](../integration/hydrobase-grhaylib-and-tmunu.md)
- See also: [Staggered State and Magnetic Reconstruction](../magnetics/staggered-state-and-magnetic-reconstruction.md)
