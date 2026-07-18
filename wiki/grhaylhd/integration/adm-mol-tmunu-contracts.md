# ADM, MoL, and Tmunu Contracts

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Integration](index.md)

## Scope and Non-Scope

This page records GRHayLHD's visible ADMBase reads, MoL registration calls,
declared aliases, and optional additive TmunuBase writes. External ADMBase,
MoL, Carpet, TmunuBase, and scheduler semantics remain unverified.

## Summary

GRHayLHD visibly calls MoL APIs with core and conditional evolved/RHS group
indices, HydroBase and optional TmunuBase group indices, and four ADMBase group
indices. Tmunu computation reads ADM fields, base thermodynamics, native
velocity, and `u0`; zeros `BU`; delegates tensor computation; then adds ten
components into TmunuBase fields. Setup-conditional registration versus an
always-steerable parameter and three unchecked registration results remain
open static issues.

## Variant Applicability

| Applicability | Registration and Tmunu boundary |
| --- | --- |
| Common | Calls MoL APIs with core conservative/RHS, `rho`/`press`/`eps`, ADM, and optional Tmunu group indices; Tmunu computation is additive. |
| Hybrid/Simple | No extra group-index call. |
| Hybrid/Simple+Entropy | Adds calls with entropy/RHS and HydroBase entropy indices. |
| Tabulated | Adds calls with electron-fraction/RHS and HydroBase `Y_e`/temperature indices. |
| Tabulated+Entropy | Adds both optional group-index call sets. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `INT-AMT-01` | Function visibly calls MoL APIs with core, conditional, and HydroBase group indices. | visible-implementation | Registration function | `c:GRHayLHD/src/MoL_registration.c#symbol=GRHayLHD_RegisterVars` |
| `INT-AMT-02` | Function visibly passes ADMBase lapse, shift, metric, and curvature indices to save-and-restore API calls. | visible-implementation | Registration function | `c:GRHayLHD/src/MoL_registration.c#symbol=GRHayLHD_RegisterVars` |
| `INT-AMT-03` | Schedule declares optional Tmunu computation in `AddToTmunu` with ADM, primitive, native-velocity, and `u0` reads. | declared | Tmunu schedule | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_compute_Tmunu` |
| `INT-AMT-04` | Tmunu function visibly zeros `BU`, calls `ghl_compute_TDNmunu`, and uses additive writes. | visible-implementation | Tmunu function | `c:GRHayLHD/src/compute_Tmunu.c#symbol=GRHayLHD_compute_Tmunu` |
| `INT-AMT-05` | Interface declares stored `u0` as a separate Tmunu dependency. | declared | `u0` declaration | `ccl:GRHayLHD/interface.ccl#group=u0` |
| `INT-AMT-06` | Runtime steering effects on setup-conditional Tmunu registration and scheduling are unresolved externally. | unresolved | Parameter/schedule lifecycle mismatch | `ccl:GRHayLHD/param.ccl#parameter=update_Tmunu` |

## Details

### ADMBase boundary

Schedule declarations read ADMBase lapse, shift, and metric across conversion,
recovery, boundaries, fluxes, and Tmunu. Source RHS additionally reads
curvature. `GRHayLHD_RegisterVars` visibly passes lapse, shift, metric, and
curvature group indices to `MoLRegisterSaveAndRestoreGroup`; comments state
intent to preserve them across timestep setup. Actual preservation behavior is
delegated to MoL.

### MoL registration classes

Visible API-call matrix; external registration effects remain unverified:

- Evolved with RHS: core `grmhd_conservatives`; conditional `ent_star` when
  `ghl_params->evolve_entropy`; conditional `Ye_star` when
  `ghl_eos->eos_type == ghl_eos_tabulated`.
- Constrained: HydroBase `rho`, `press`, and `eps`; conditional entropy;
  conditional `Y_e` and temperature; three TmunuBase stress-energy groups when
  `update_Tmunu` is true.
- Save-and-restore: ADMBase lapse, shift, metric, and curvature.

`u0` has declared storage and is written by Prim2Con, recovery, and boundary
paths, but is not visibly registered in this function. Interface CCL declares
aliases for `GetRefinementLevel`, `MoLRegisterEvolved`,
`MoLRegisterEvolvedGroup`, `MoLRegisterConstrainedGroup`, and
`MoLRegisterSaveAndRestoreGroup`. Local C calls group registration forms,
constrained registration, save-and-restore registration, and uses
`GetRefinementLevel` in boundary/diagnostic paths; it does not visibly call
`MoLRegisterEvolved`.

Registration uses `GRHayLHD::grmhd_conservatives`. Three variant recovery
symmetry branches visibly request `GRHayLHD::grhd_conservatives`; runtime
effect remains unverified under [GRH-0002](../contradictions.md#grh-0002).

### Tmunu dataflow

When `update_Tmunu` condition is met, CCL schedules function in `AddToTmunu`
and declares reads of ADM metric/lapse/shift, HydroBase `rho`/`press`/`eps`,
native velocity, and `u0`. Function initializes metric and auxiliary objects,
loads those primitive fields, explicitly assigns all `prims.BU` components
zero, and calls `ghl_compute_TDNmunu`.

Ten visible assignments use `+=` for `eTtt`, `eTtx`, `eTty`, `eTtz`, `eTxx`,
`eTxy`, `eTxz`, `eTyy`, `eTyz`, and `eTzz`. This proves additive local writes,
not initialization, tensor semantics, or successful contribution ordering.

### Open lifecycle and status handling

`update_Tmunu` is declared `STEERABLE=ALWAYS`, while Tmunu scheduling and
registration are both setup-conditional. Whether runtime changes reconfigure
either boundary is unresolved; see
[GRH-0007](../contradictions.md#grh-0007).

Most registration return values are accumulated into `ierr` and checked.
Three Tmunu constrained-registration calls visibly omit that accumulation;
they remain unchecked-status review candidates under
[GRH-0008](../contradictions.md#grh-0008), not observed failures.

## Caveats

- CCL `READS`/`WRITES` and source calls establish declared or visible
  boundaries only, not framework execution.
- Additive tensor writes require external lifecycle context to determine prior
  initialization and contribution ordering.
- Zero `BU` is local to Tmunu primitive setup; no magnetic evolution is
  inferred.
- Generated group indices, alias binding, and return-code semantics are
  external.

## Sources

- [Interface declarations](../../../GRHayLHD/interface.ccl)
- [Parameter declarations](../../../GRHayLHD/param.ccl)
- [Schedule declarations](../../../GRHayLHD/schedule.ccl)
- [MoL registration](../../../GRHayLHD/src/MoL_registration.c)
- [Tmunu computation](../../../GRHayLHD/src/compute_Tmunu.c)

## Related Pages

- [Variables and Storage](../architecture/variables-and-storage.md)
- [Declared Schedule Lifecycle](../architecture/schedule-lifecycle.md)
- [GRHayLib Contract](grhaylib-contract.md)
