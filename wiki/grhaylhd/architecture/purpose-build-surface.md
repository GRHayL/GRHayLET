# Purpose and Build Surface

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Architecture](index.md)

## Scope and Non-Scope

This page records GRHayLHD's stated purpose, Cactus surface, dependency
declarations, and five checked-in build manifests. It does not establish a
successful build, link result, runtime activation, or external-library
behavior.

## Summary

GRHayLHD documentation describes a general-relativistic hydrodynamics thorn
built on GRHayL through GRHayLib. The ThornGuide calls it a trimmed-down
IllinoisGRMHD form whose magnetic fields are set to zero and says it evolves
systems without magnetic fields. Local manifests list seven common C units and
the same seven unit names in each of four variant directories.

## Variant Applicability

| Applicability | Build family | Locally stated role |
| --- | --- | --- |
| Common | `GRHayLHD/src/make.code.defn` | Lists common support units and four variant subdirectories. |
| Hybrid/Simple | `Hybrid` | Lists no-entropy Hybrid/Simple operation units. |
| Hybrid/Simple+Entropy | `HybridEntropy` | Lists entropy Hybrid/Simple operation units. |
| Tabulated | `Tabulated` | Lists no-entropy tabulated operation units. |
| Tabulated+Entropy | `TabulatedEntropy` | Lists entropy tabulated operation units. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ARCH-PURPOSE-01` | README states GRHD purpose and GRHayLib use. | declared | Purpose section | `doc:GRHayLHD/README#section=1. Purpose` |
| `ARCH-PURPOSE-02` | ThornGuide describes magnetic fields as zero and not evolved. | declared | Introduction | `doc:GRHayLHD/doc/documentation.tex#section=Introduction` |
| `ARCH-SURFACE-01` | Interface CCL declares GRHayLHD implementation plus inherited/include/function boundaries. | declared | Complete interface declaration | `ccl:GRHayLHD/interface.ccl#implementation=GRHayLHD` |
| `ARCH-SURFACE-02` | Configuration CCL declares an HDF5 requirement. | declared | Requirement declaration | `ccl:GRHayLHD/configuration.ccl#requirement=HDF5` |
| `ARCH-BUILD-01` | Root manifest lists seven common units. | declared | `SRCS` field | `build:GRHayLHD/src/make.code.defn#field=SRCS` |
| `ARCH-BUILD-02` | Hybrid manifest lists seven variant units. | declared | `SRCS` field | `build:GRHayLHD/src/Hybrid/make.code.defn#field=SRCS` |
| `ARCH-BUILD-03` | HybridEntropy manifest lists seven variant units. | declared | `SRCS` field | `build:GRHayLHD/src/HybridEntropy/make.code.defn#field=SRCS` |
| `ARCH-BUILD-04` | Tabulated manifest lists seven variant units. | declared | `SRCS` field | `build:GRHayLHD/src/Tabulated/make.code.defn#field=SRCS` |
| `ARCH-BUILD-05` | TabulatedEntropy manifest lists seven variant units. | declared | `SRCS` field | `build:GRHayLHD/src/TabulatedEntropy/make.code.defn#field=SRCS` |

## Details

### Thorn and dependency declarations

[interface.ccl](../../../GRHayLHD/interface.ccl) declares implementation name
`GRHayLHD`; inheritance from `ADMBase`, `Tmunubase`, `HydroBase`, and
`GRHayLib`; includes `Symmetry.h` and `GRHayLib.h`; and four used MoL function
interfaces plus `GetRefinementLevel`. [configuration.ccl](../../../GRHayLHD/configuration.ccl)
declares `requires HDF5`. These are declarations only. Neither file proves
resolved headers, aliases, library linkage, or runtime behavior.

### Checked-in compilation surface

Root `SUBDIRS` names `Hybrid`, `HybridEntropy`, `Tabulated`, and
`TabulatedEntropy`. Root `SRCS` lists metric derivative/interpolation support,
stress-energy computation, both HydroBase velocity converters, symmetry setup,
and MoL registration. Every variant `SRCS` lists
`conservs_to_prims.c`, `evaluate_fluxes_rhs.c`,
`evaluate_sources_rhs.c`, `outer_boundaries.c`,
`perturb_conservatives.c`, `perturb_primitives.c`, and
`prims_to_conservs.c`. Identical filenames do not imply identical bodies; all
four families require independent inspection.

### Non-magnetic contract

Documentation supplies purpose-level authority. Visible C routines strengthen
only a narrower implementation observation: routines that construct local
`ghl_primitive_quantities` repeatedly set `BU[0..2]` to zero. No local
interface group, storage clause, scheduled operation, or build unit introduces
magnetic evolution. Magnetic reconstruction, induction, gauge, and EM
boundaries therefore remain outside this KB branch.

## Caveats

- GRHayLib feature semantics and HDF5 discovery/link behavior are delegated and
  unverified externally.
- README says GRHayLHD supports "most features" and "all Con2Prim routines";
  these are statements of intent, not an exhaustive local API proof.
- Variant boundary headers contain stale magnetic-stage wording; see
  [GRH-0005](../contradictions.md#grh-0005).
- Static manifests prove checked-in source listing only, not compilation.

## Sources

- [README](../../../GRHayLHD/README)
- [ThornGuide source](../../../GRHayLHD/doc/documentation.tex)
- [Interface declarations](../../../GRHayLHD/interface.ccl)
- [Configuration declaration](../../../GRHayLHD/configuration.ccl)
- [Root build manifest](../../../GRHayLHD/src/make.code.defn)
- Variant build manifests: [Hybrid](../../../GRHayLHD/src/Hybrid/make.code.defn),
  [HybridEntropy](../../../GRHayLHD/src/HybridEntropy/make.code.defn),
  [Tabulated](../../../GRHayLHD/src/Tabulated/make.code.defn), and
  [TabulatedEntropy](../../../GRHayLHD/src/TabulatedEntropy/make.code.defn)

## Related Pages

- [Variables and Storage](variables-and-storage.md)
