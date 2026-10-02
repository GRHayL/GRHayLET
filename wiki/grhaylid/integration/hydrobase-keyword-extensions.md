# HydroBase Keyword Extensions

> Page status: reviewed · Last reviewed: 10-02-2026
> Up: [Integration](index.md)

## Scope and Non-Scope

This page owns the local declarations and visible dataflow described below.
Only GRHayLID sources are domain evidence. Framework execution, library
semantics, production-table provenance, and numerical validation remain
out of scope.

## Summary

Six HydroBase keywords are extended: initial_hydro gains HydroTest1D,
IsotropicGas, and ConstantDensitySphere; initial_Y_e, initial_temperature,
initial_entropy, initial_Avec, and initial_Bvec each gain GRHayLID.
ParamCheck now validates their complete local producer matrix.

## Mode Applicability

| Applicability | Local surface |
| --- | --- |
| Common | Six keyword extensions with parameter-check and body preconditions. |
| BetaEquilibrium | Ye/T may be supplied by beta after an external hydro owner. |
| HydroTest1D+Magnetic | Both GRHayLID selectors required. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `HYDROBASE-KEYWORD-EXTENSIONS-01` | initial_hydro has a local HydroBase extension. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=initial_hydro` |
| `HYDROBASE-KEYWORD-EXTENSIONS-02` | initial_Y_e has a local HydroBase extension. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=initial_Y_e` |
| `HYDROBASE-KEYWORD-EXTENSIONS-03` | initial_temperature has a local HydroBase extension. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=initial_temperature` |
| `HYDROBASE-KEYWORD-EXTENSIONS-04` | initial_entropy has a local HydroBase extension. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=initial_entropy` |
| `HYDROBASE-KEYWORD-EXTENSIONS-05` | initial_Avec has a local HydroBase extension. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=initial_Avec` |
| `HYDROBASE-KEYWORD-EXTENSIONS-06` | initial_Bvec has a local HydroBase extension. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=initial_Bvec` |
| `HYDROBASE-KEYWORD-EXTENSIONS-07` | ParamCheck diagnoses missing/incompatible local producers. | visible-implementation | Named local source | `c:GRHayLID/src/ParamCheck.c#symbol=GRHayLID_ParamCheck` |

## Details

| Selection | Required local producer / contract |
| --- | --- |
| HydroTest1D | Simple or Hybrid hydro body |
| IsotropicGas / ConstantDensitySphere | Tabulated EOS; both Ye/T selectors GRHayLID |
| GRHayLID Avec/Bvec | HydroTest1D, enabled magnetics, both selectors GRHayLID |
| GRHayLID Ye/T | Gas, sphere, or enabled beta |
| GRHayLID entropy | Non-none hydro producer; supported EOS and initialized inputs |

Standalone beta intentionally replaces Ye/T after HydroBase_Initial and requests
storage even if their selectors are none. Tabulated entropy without beta requires
non-none Ye/T selectors; README requires their producers to initialize fields in
HydroBase_Initial. Native entropy requires compatible-consumer opt-in. No local
body initializes w_lorentz; a reader must establish another producer and order it.

Schedule reads/writes name shared HydroBase groups. Native tests additionally
read ADMBase::metric. Keyword/storage declarations alone do not establish actual
framework acceptance or initialized contents. GID-0001/0012/0015 preserve the
resolved historical selector and control mismatches.

## Caveats

Storage and schedule declarations do not prove allocation or execution.
Local guards and calls do not establish external error or interpolation
semantics. No checked-in GRHayLID test/parfile/oracle validates this path.

## Sources

- [param.ccl](../../../GRHayLID/param.ccl)
- [ParamCheck.c](../../../GRHayLID/src/ParamCheck.c)

## Related Pages

- [GRHayLib Contract](grhaylib-contract.md)
- [Parameters](parameters-and-configurations.md)
- [Entropy](../initial-data/entropy-computation.md)
- [GID-0001](../contradictions.md#gid-0001)
- [GID-0012](../contradictions.md#gid-0012)
- [GID-0015](../contradictions.md#gid-0015)
