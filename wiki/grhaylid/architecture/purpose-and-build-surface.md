# Purpose and Build Surface

> Page status: reviewed · Last reviewed: 10-02-2026
> Up: [Architecture](index.md)

## Scope and Non-Scope

This page owns the local declarations and visible dataflow described below.
Only GRHayLID sources are domain evidence. Framework execution, library
semantics, production-table provenance, and numerical validation remain
out of scope.

## Summary

README enumerates Balsara, equilibrium, shock tube, sound-wave velocity
perturbation, gas, and sphere setups, plus standalone beta/entropy features.
Interface CCL inherits GRHayLib, Grid, HydroBase, and ADMBase. The manifest
lists seven C units including ParamCheck.c. No private gridfunction group or
MoL registration is declared; beta scheduling now declares storage for the
HydroBase Ye and temperature groups.

## Mode Applicability

| Applicability | Local surface |
| --- | --- |
| Common | Purpose, seven-unit build surface, double precision, and external state boundaries. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `PURPOSE-AND-BUILD-SURFACE-01` | README states the supported families and explicit normalization/entropy contracts. | declared | Named local source | `doc:GRHayLID/README#section=1. Purpose` |
| `PURPOSE-AND-BUILD-SURFACE-02` | Interface declares four inherited implementations including ADMBase. | declared | Named local source | `ccl:GRHayLID/interface.ccl#implementation=GRHayLID` |
| `PURPOSE-AND-BUILD-SURFACE-03` | Configuration declares HDF5. | declared | Named local source | `ccl:GRHayLID/configuration.ccl#requirement=HDF5` |
| `PURPOSE-AND-BUILD-SURFACE-04` | SRCS lists seven C translation units. | declared | Named local source | `build:GRHayLID/src/make.code.defn#field=SRCS` |
| `PURPOSE-AND-BUILD-SURFACE-05` | Common header rejects builds without CCTK_REAL_PRECISION_8. | visible-implementation | Named local source | `macro:GRHayLID/src/GRHayLID.h#include=GRHayLib.h` |
| `PURPOSE-AND-BUILD-SURFACE-06` | Beta is conditional and writes HydroBase fields. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_BetaEquilibrium` |

## Details

The compilation units are 1D_tests_hydro_data.c,
1D_tests_magnetic_data.c, BetaEquilibrium.c, ComputeEntropy.c,
ConstantDensitySphere.c, IsotropicGas.c, and ParamCheck.c. The header is
included textually and is not a separate translation unit.

ADMBase supplies the metric read by native test routines. The local header
checks finite identity-metric components, table-state effective bounds,
active group storage, and the required eight-byte real configuration.
CHECK_PARAMETER retains the local -1 sentinel convention.

The thorn declares no private state. Conditional beta STORAGE requests act
on shared HydroBase fields; absence of a private group does not imply absence
of local allocation requests. No local body produces w_lorentz.

ThornGuide now describes both magnetic outputs, the actual entropy selector,
and the locally implemented setup inventory. Resolved historical mismatches
remain recorded in GID-0002, GID-0003, and GID-0009.

## Caveats

Storage and schedule declarations do not prove allocation or execution.
Local guards and calls do not establish external error or interpolation
semantics. No checked-in GRHayLID test/parfile/oracle validates this path.

## Sources

- [README](../../../GRHayLID/README)
- [interface.ccl](../../../GRHayLID/interface.ccl)
- [configuration.ccl](../../../GRHayLID/configuration.ccl)
- [make.code.defn](../../../GRHayLID/src/make.code.defn)
- [GRHayLID.h](../../../GRHayLID/src/GRHayLID.h)
- [schedule.ccl](../../../GRHayLID/schedule.ccl)

## Related Pages

- [Schedule Lifecycle](schedule-lifecycle.md)
- [Initial Data](../initial-data/index.md)
- [GID-0002](../contradictions.md#gid-0002)
- [GID-0003](../contradictions.md#gid-0003)
- [GID-0009](../contradictions.md#gid-0009)
