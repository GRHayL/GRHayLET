# Purpose and Build Surface

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Architecture](index.md)

## Scope and Non-Scope

This page records GRHayLID's locally stated purpose, Cactus interface and
configuration declarations, checked-in compilation manifest, and common header
surface. It does not establish build success, resolved dependencies, scheduled
execution, numerical validity, or behavior supplied by Cactus, HDF5, HydroBase,
or GRHayLib.

## Summary

README presents GRHayLID as an initial-data thorn for GR(M)HD systems using
GRHayL. It enumerates six setup families: five Balsara tests, static
equilibrium, shock tube, and sound wave under one-dimensional tests, plus
IsotropicGas and ConstantDensitySphere. Magnetic initialization is described
as optional. Beta-equilibrium imposition and entropy computation are described
as standalone features usable when another thorn supplies initial data.

Checked-in declarations implement `GRHayLID`, inherit `GRHayLib`, `Grid`, and
`HydroBase`, request `GRHayLib.h`, and require HDF5. One `SRCS` manifest lists
six C units. Local interface and schedule CCL declare no GRHayLID gridfunction
group, storage clause, or MoL registration; every declared write names a
HydroBase field.

## Mode Applicability

| Applicability | Locally stated surface |
| --- | --- |
| Common | Purpose, interface, dependency, build-manifest, no-owned-storage, and header declarations span all GRHayLID modes. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ARCH-PURPOSE-01` | README declares six initial-data setup families, optional magnetic initialization, and standalone beta-equilibrium and entropy features. | declared | Purpose section | `doc:GRHayLID/README#section=1. Purpose` |
| `ARCH-PURPOSE-02` | ThornGuide Introduction declares five Balsara tests, equilibrium, sound wave, shock tube, IsotropicGas, ConstantDensitySphere, and optional beta equilibrium. | declared | Introduction | `doc:GRHayLID/doc/documentation.tex#section=Introduction` |
| `ARCH-SURFACE-01` | Interface CCL implements GRHayLID, inherits three named interfaces, and requests `GRHayLib.h`. | declared | Complete interface declaration | `ccl:GRHayLID/interface.ccl#implementation=GRHayLID` |
| `ARCH-SURFACE-02` | Configuration CCL declares an HDF5 requirement. | declared | Requirement declaration | `ccl:GRHayLID/configuration.ccl#requirement=HDF5` |
| `ARCH-BUILD-01` | Build manifest lists six C translation units in `SRCS`. | declared | `SRCS` field | `build:GRHayLID/src/make.code.defn#field=SRCS` |
| `ARCH-HEADER-01` | Common header visibly defines `CHECK_PARAMETER` as an equality check followed by a `CCTK_VERROR` call. | visible-implementation | Macro definition | `macro:GRHayLID/src/GRHayLID.h#name=CHECK_PARAMETER` |

## Details

### Purpose inventory

[README](../../../GRHayLID/README) groups its setup inventory into six numbered
items: all five Balsara tests, static equilibrium, shock tube, sound wave,
IsotropicGas, and ConstantDensitySphere. It labels the first four numbered
items simple/hybrid-only and the final two tabulated-only. It separately states
that magnetic fields are optional and that beta equilibrium and entropy can be
applied when another initial-data thorn supplies the base state. These are
documentation claims; individual implementation and dispatch details belong to
Initial Data pages.

[ThornGuide source](../../../GRHayLID/doc/documentation.tex) divides the
inventory into one-dimensional and three-dimensional sets. Its abstract also
advertises a cylindrical explosion, but no section heading, parameter choice,
schedule declaration, or implementation unit names such a setup. That local
documentation mismatch remains open under GID-0003.

### Interface and build surface

[interface.ccl](../../../GRHayLID/interface.ccl) contains only implementation,
inheritance, and include-use declarations. It declares `implements: GRHayLID`,
`inherits: GRHayLib, Grid, HydroBase`, and `USES INCLUDE: GRHayLib.h`.
[configuration.ccl](../../../GRHayLID/configuration.ccl) consists of the
declaration `requires HDF5`. This requirement does not prove header discovery,
linkage, or any HDF5 behavior.

[make.code.defn](../../../GRHayLID/src/make.code.defn) lists
`1D_tests_hydro_data.c`, `1D_tests_magnetic_data.c`, `BetaEquilibrium.c`,
`ComputeEntropy.c`, `ConstantDensitySphere.c`, and `IsotropicGas.c`. It does not
list `GRHayLID.h` in `SRCS`; the C units include that header textually. The
manifest proves only the checked-in source listing.

### State ownership and common header

Complete local inspection of `interface.ccl` and `schedule.ccl` finds no
GRHayLID gridfunction-group declaration, `STORAGE` clause, or MoL registration.
Schedule `WRITES` clauses name only `HydroBase::rho`, `press`, `eps`, `vel`,
`Y_e`, `temperature`, `entropy`, `Avec`, and `Bvec`. Thus "no owned storage"
means absence from these checked-in declarations, not a framework-level claim
about allocation or alias behavior.

[GRHayLID.h](../../../GRHayLID/src/GRHayLID.h) includes Cactus headers and
`GRHayLib.h`. Its `CHECK_PARAMETER(par)` macro visibly compares `par` with
`-1` and calls `CCTK_VERROR` when equal. Whether that external macro terminates
or how it reports errors is outside local authority.

## Caveats

- Static declarations do not prove a successful build, dependency resolution,
  storage allocation, or execution.
- ThornGuide says only `Avec` is set and that `B` variables are not set, while
  checked-in magnetic declarations and code include `Bvec`; see
  [GID-0002](../contradictions.md#gid-0002).
- Cylindrical explosion appears only in ThornGuide abstract prose; see
  [GID-0003](../contradictions.md#gid-0003).
- README's "any EOS" entropy wording is broader than two locally declared
  entropy schedule arms; see [GID-0009](../contradictions.md#gid-0009).
- HDF5 and external include/library semantics remain out of scope.

## Sources

- [README](../../../GRHayLID/README)
- [ThornGuide source](../../../GRHayLID/doc/documentation.tex)
- [Interface declarations](../../../GRHayLID/interface.ccl)
- [Configuration declaration](../../../GRHayLID/configuration.ccl)
- [Schedule declarations](../../../GRHayLID/schedule.ccl)
- [Build manifest](../../../GRHayLID/src/make.code.defn)
- [Common header](../../../GRHayLID/src/GRHayLID.h)

## Related Pages

- [Declared Schedule Lifecycle](schedule-lifecycle.md)
- [Initial Data](../initial-data/index.md)
- [GID-0002](../contradictions.md#gid-0002)
- [GID-0003](../contradictions.md#gid-0003)
- [GID-0009](../contradictions.md#gid-0009)
