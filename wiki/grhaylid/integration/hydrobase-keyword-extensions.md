# HydroBase Keyword Extensions

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Integration](index.md)

## Scope and Non-Scope

This page owns the six HydroBase keyword extensions declared by GRHayLID,
their local selection sites, and HydroBase gridfunctions named in local
schedule declarations. It does not establish HydroBase variable semantics,
keyword acceptance, storage, or execution. `EOS_type` belongs to the
[GRHayLib Contract](grhaylib-contract.md), not this page.

## Summary

Parameter CCL shares HydroBase and extends `initial_hydro`, `initial_Y_e`,
`initial_temperature`, `initial_entropy`, `initial_Avec`, and `initial_Bvec`.
`initial_hydro` gains three GRHayLID setup values; each other keyword gains
`GRHayLID`. Schedule CCL consumes `initial_hydro` and `initial_entropy`
directly. C bodies check the electron-fraction, temperature, and magnetic
selections. Local schedules declare writes only to HydroBase state.

## Mode Applicability

| Applicability | Keyword boundary |
| --- | --- |
| Common | Six shared HydroBase keywords are extended; their declarations do not establish external acceptance or storage. |
| HydroTest1D | `initial_hydro="HydroTest1D"` selects hydro setup; magnetic body checks either `initial_Avec` or `initial_Bvec` for `GRHayLID`. |
| HydroTest1D+Magnetic | Magnetic schedule is additionally guarded by local `initialize_magnetic_quantities`; its body checks the two magnetic keyword values. |
| IsotropicGas | `initial_hydro="IsotropicGas"`; body checks `initial_Y_e` and `initial_temperature`. |
| ConstantDensitySphere | `initial_hydro="ConstantDensitySphere"`; body checks `initial_Y_e` and `initial_temperature`. |
| BetaEquilibrium | Feature has no `initial_hydro` requirement; it operates on HydroBase density and declares four HydroBase writes. |
| Entropy/Hybrid | `initial_entropy="GRHayLID"` participates in schedule dispatch; hybrid branch declares entropy output. |
| Entropy/Tabulated | Same keyword participates; tabulated branch declares six HydroBase writes. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `INT-HB-01` | `initial_hydro` is extended with `HydroTest1D`, `IsotropicGas`, and `ConstantDensitySphere`. | declared | HydroBase shared-keyword extension | `ccl:GRHayLID/param.ccl#parameter=initial_hydro` |
| `INT-HB-02` | `initial_Y_e` is extended with `GRHayLID`. | declared | HydroBase shared-keyword extension | `ccl:GRHayLID/param.ccl#parameter=initial_Y_e` |
| `INT-HB-03` | `initial_temperature` is extended with `GRHayLID`. | declared | HydroBase shared-keyword extension | `ccl:GRHayLID/param.ccl#parameter=initial_temperature` |
| `INT-HB-04` | `initial_entropy` is extended with `GRHayLID`. | declared | HydroBase shared-keyword extension | `ccl:GRHayLID/param.ccl#parameter=initial_entropy` |
| `INT-HB-05` | `initial_Avec` is extended with `GRHayLID`. | declared | HydroBase shared-keyword extension | `ccl:GRHayLID/param.ccl#parameter=initial_Avec` |
| `INT-HB-06` | `initial_Bvec` is extended with `GRHayLID`. | declared | HydroBase shared-keyword extension | `ccl:GRHayLID/param.ccl#parameter=initial_Bvec` |
| `INT-HB-07` | Interface CCL declares the GRHayLID implementation and inheritance from GRHayLib, Grid, and HydroBase. | declared | Local interface declaration | `ccl:GRHayLID/interface.ccl#implementation=GRHayLID` |
| `INT-HB-08` | ThornGuide describes the keyword controls and claims an entropy control that differs from local CCL. | declared | Parameters section | `doc:GRHayLID/doc/documentation.tex#section=Parameters` |
| `INT-HB-09` | IsotropicGas visibly checks the electron-fraction and temperature keyword values before writing HydroBase fields. | visible-implementation | Local guard site | `c:GRHayLID/src/IsotropicGas.c#symbol=GRHayLID_IsotropicGas` |

## Details

### Extension and consumption map

| Extended keyword | Added values | Visible local consumer |
| --- | --- | --- |
| `initial_hydro` | `HydroTest1D`, `IsotropicGas`, `ConstantDensitySphere` | Top-level schedule dispatch; the `HydroTest1D` description names undeclared `test_1D_initial_data`, while the local selector is `initial_data_1D`; see [GID-0015](../contradictions.md#gid-0015) |
| `initial_Y_e` | `GRHayLID` | `IsotropicGas.c` and `ConstantDensitySphere.c` preconditions |
| `initial_temperature` | `GRHayLID` | `IsotropicGas.c` and `ConstantDensitySphere.c` preconditions |
| `initial_entropy` | `GRHayLID` | Hybrid/Tabulated schedule dispatch |
| `initial_Avec` | `GRHayLID` | Magnetic function either/or precondition |
| `initial_Bvec` | `GRHayLID` | Magnetic function either/or precondition |

The magnetic precondition accepts either magnetic keyword, while the body
unconditionally writes both arrays. This declaration/implementation hazard is
tracked as [GID-0012](../contradictions.md#gid-0012). ThornGuide instead names
a nonexistent local `compute_entropy` control; [GID-0001](../contradictions.md#gid-0001)
records the mismatch with `initial_entropy` plus `EOS_type` dispatch.
The `HydroTest1D` extension description instead directs users to
`test_1D_initial_data`, while the restricted declaration names
`initial_data_1D`; [GID-0015](../contradictions.md#gid-0015) records that
literal name drift without inferring an external alias.

### Declared HydroBase writes

Across seven schedule blocks, local CCL names `rho`, `press`, `eps`, `vel`,
`Y_e`, `temperature`, `entropy`, `Avec`, and `Bvec` as HydroBase writes.
Beta equilibrium reads `rho`; Hybrid entropy reads `rho` and `press`;
Tabulated entropy reads `rho`, `Y_e`, and `temperature`. These are declared
accesses only. Interface CCL declares no local gridfunction group, so this
page does not attribute storage ownership to GRHayLID.

### Guard locations

`initial_hydro` and `initial_entropy` are schedule conditions.
`initial_Y_e` and `initial_temperature` are checked in each 3D setup body.
`initial_Avec` and `initial_Bvec` are checked in the magnetic body, while the
separate local boolean `initialize_magnetic_quantities` guards whether that
body is scheduled. External error behavior and keyword validation are not
inferred.

## Caveats

- `EXTENDS KEYWORD` declares accepted text locally; Cactus/HydroBase handling
  is external.
- Schedule READS/WRITES do not establish allocation, storage ownership, or
  execution.
- Error messages and checks expose intent but do not prove termination.
- The `HydroTest1D` description's parameter reference differs from the local
  selector declaration; see [GID-0015](../contradictions.md#gid-0015).
- `EOS_type` is a GRHayLib shared keyword and is intentionally documented on
  another page.

## Sources

- [Parameter and keyword declarations](../../../GRHayLID/param.ccl)
- [Interface declaration](../../../GRHayLID/interface.ccl)
- [Schedule declarations](../../../GRHayLID/schedule.ccl)
- [ThornGuide parameters](../../../GRHayLID/doc/documentation.tex)
- [IsotropicGas guard site](../../../GRHayLID/src/IsotropicGas.c)

## Related Pages

- [GRHayLib Contract](grhaylib-contract.md)
- [Parameters and Configurations](parameters-and-configurations.md)
- [Entropy Computation](../initial-data/entropy-computation.md)
- [One-D Magnetic Tests](../initial-data/one-d-tests-magnetic.md)
- [GID-0001: entropy-control mismatch](../contradictions.md#gid-0001)
- [GID-0012: magnetic selection/write hazard](../contradictions.md#gid-0012)
- [GID-0015: HydroTest1D parameter-name mismatch](../contradictions.md#gid-0015)
