# Variables and Storage

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Architecture](index.md)

## Scope and Non-Scope

This page maps locally declared GRHayLHD state, storage conditions, interface
tags, and visible MoL registration consumers. HydroBase, ADMBase, and
TmunuBase internals remain external.

## Summary

GRHayLHD owns native coordinate velocities, `u0`, five core conservatives,
their RHS and flux temporary groups, optional entropy and electron-fraction
state, and `failure_checker`. HydroBase owns thermodynamic primitives used by
GRHayLHD; TmunuBase owns stress-energy outputs. `ent_star` storage follows
`evolve_entropy`; `Ye_star` storage follows tabulated `EOS_type`.

## Variant Applicability

| Applicability | Extra evolved state | Extra HydroBase primitives | Extra RHS and flux state |
| --- | --- | --- | --- |
| Common | Five core conservatives | `rho`, `press`, `eps`; native `vx`, `vy`, `vz` remain GRHayLHD-owned | Core conservative RHS and flux temporaries |
| Hybrid/Simple | None | Base thermodynamics | None |
| Hybrid/Simple+Entropy | `ent_star` | `entropy` | `ent_star_rhs`, `ent_star_flux` |
| Tabulated | `Ye_star` | `Y_e`, `temperature` | `Ye_star_rhs`, `Ye_star_flux` |
| Tabulated+Entropy | `ent_star`, `Ye_star` | `entropy`, `Y_e`, `temperature` | Both optional RHS and flux pairs |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ARCH-VAR-01` | Interface declares native `vx`, `vy`, and `vz`. | declared | Velocity group | `ccl:GRHayLHD/interface.ccl#group=grmhd_velocities` |
| `ARCH-VAR-02` | Interface declares five core conservatives with three timelevels. | declared | Conservative group | `ccl:GRHayLHD/interface.ccl#group=grmhd_conservatives` |
| `ARCH-VAR-03` | Interface declares entropy evolved state. | declared | Entropy group | `ccl:GRHayLHD/interface.ccl#group=ent_star` |
| `ARCH-VAR-04` | Interface declares electron-fraction evolved state. | declared | Electron-fraction group | `ccl:GRHayLHD/interface.ccl#group=Ye_star` |
| `ARCH-VAR-05` | Interface declares diagnostic `failure_checker`. | declared | Diagnostic group | `ccl:GRHayLHD/interface.ccl#group=failure_checker` |
| `ARCH-VAR-06` | Schedule CCL conditionally stores electron-fraction state for Tabulated EOS. | declared | `Ye_star` storage clause | `ccl:GRHayLHD/schedule.ccl#storage=Ye_star` |
| `ARCH-VAR-07` | Function visibly calls MoL APIs with optional evolved/RHS and HydroBase group indices. | visible-implementation | Registration function | `c:GRHayLHD/src/MoL_registration.c#symbol=GRHayLHD_RegisterVars` |
| `ARCH-VAR-08` | Schedule CCL conditionally stores entropy state when `evolve_entropy` is true. | declared | `ent_star` storage clause | `ccl:GRHayLHD/schedule.ccl#storage=ent_star` |

## Details

### Owned state and tags

| State | Declaration | Storage metadata | Owner / visible consumers |
| --- | --- | --- | --- |
| `vx`, `vy`, `vz` | `grmhd_velocities` | GF, `InterpNumTimelevels=1`, no prolongation | GRHayLHD; conversion, conversion/recovery, RHS, BC, perturbation, Tmunu |
| `u0` | scalar GF | `InterpNumTimelevels=1`, no prolongation, no checkpoint | GRHayLHD; conversion/recovery/BC write, Tmunu read |
| `rho_star`, `tau`, `Stildex`, `Stildey`, `Stildez` | `grmhd_conservatives` | three timelevels, ENO prolongation | GRHayLHD; MoL-evolved against core RHS |
| `ent_star` | scalar GF | three timelevels, ENO prolongation | GRHayLHD; conditional MoL evolution |
| `Ye_star` | scalar GF | three timelevels, ENO prolongation | GRHayLHD; conditional MoL evolution |
| Core RHS | `grmhd_conservatives_rhs` | no prolongation, no checkpoint | Source/flux writers; MoL RHS |
| Optional RHS | `ent_star_rhs`, `Ye_star_rhs` | no prolongation, no checkpoint | Conditional source/flux writers; MoL RHS |
| Core flux temporary | `grmhd_flux_temps` | no prolongation, no checkpoint | Variant flux routines |
| Optional flux | `ent_star_flux`, `Ye_star_flux` | no prolongation, no checkpoint | Conditional variant flux routines |
| `failure_checker` | scalar diagnostic GF | no prolongation/checkpoint; one interpolation timelevel | Recovery writers; interpretation has open mismatch |

`HydroBase::rho`, `press`, and `eps` indices are visibly passed to constrained-
group registration APIs. Entropy adds `HydroBase::entropy`; tabulated modes add
`HydroBase::Y_e` and `temperature`. When `update_Tmunu` is true, function calls
same API with three TmunuBase group indices. ADMBase lapse, shift, metric, and
curvature indices are passed to save-and-restore registration APIs.

### Conditional storage

Top-level schedule declarations always store core conservatives, native
velocities, `u0`, core RHS/flux temporaries, and `failure_checker`. A tabulated
EOS condition adds `Ye_star[3]`, its RHS, and flux. `evolve_entropy` adds
`ent_star[3]`, its RHS, and flux independently. Conditional registration in
`GRHayLHD_RegisterVars` mirrors these two axes through GRHayLib-owned runtime
objects `ghl_eos` and `ghl_params`; equivalence between CCL conditions and
external object initialization is not proved locally.

### Interface mismatches

`InitSymBound.c` asks for `GRHayLHD::Stilde_z` in a dormant equatorial branch,
while interface state spells the variable `Stildez`. Three recovery variants
also contain `grhd_conservatives` where interface declares
`grmhd_conservatives`. These strings are recorded without a runtime claim.

## Caveats

- Conditional `STORAGE` and visible MoL calls establish declarations/dataflow,
  not allocation or registration success.
- Generated Cactus argument macros and external group semantics are not local
  evidence.
- See [GRH-0003](../contradictions.md#grh-0003) for momentum spelling and
  [GRH-0002](../contradictions.md#grh-0002) for conservative-group spelling,
  and [GRH-0004](../contradictions.md#grh-0004) for diagnostic overwrite
  context.

## Sources

- [Interface declarations](../../../GRHayLHD/interface.ccl)
- [Schedule and storage declarations](../../../GRHayLHD/schedule.ccl)
- [MoL registration](../../../GRHayLHD/src/MoL_registration.c)
- [Symmetry initialization](../../../GRHayLHD/src/InitSymBound.c)

## Related Pages

- [Purpose and Build Surface](purpose-build-surface.md)
