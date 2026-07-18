# Cactus Surface and Build

> Status: confirmed · Last reconciled: 07-17-2026
> Up: [Architecture](index.md)

## Summary

IllinoisGRMHD declares one implemented thorn interface, four inherited thorns,
two used includes, HDF5/GRHayL requirements, CCL grid-function and function-alias
surfaces, and a build split between common sources and four parallel variant
subdirectories. These declarations show local inputs and exposed surfaces, not
a successful or transitively complete build.

## Detail

### Thorn and Alias Surface

`configuration.ccl` declares requirements `HDF5 GRHayL`.
`interface.ccl` declares `implements: IllinoisGRMHD`; inherits `ADMBase`,
`Tmunubase`, `HydroBase`, and `GRHayLib`; and uses `Symmetry.h` and
`GRHayLib.h`. It declares used aliases `Driver_SelectVarForBC`,
`Driver_SelectGroupForBC`, `GetRefinementLevel`, `MoLRegisterEvolved`,
`MoLRegisterEvolvedGroup`, `MoLRegisterConstrainedGroup`, and
`MoLRegisterSaveAndRestoreGroup`.

Claim evidence:
- Claim: These are locally declared requirements, inheritance, includes, and aliases; they do not prove availability or behavior of external dependencies.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/configuration.ccl`, `requires`; registered `IllinoisGRMHD/interface.ccl`, `implements`, `inherits`, `USES INCLUDE`, and function declarations
- Corroboration: registered `IllinoisGRMHD/src/MoL_registration.c`, `IllinoisGRMHD_RegisterVars`, contains local alias call sites

### Declared Groups

All groups are CCTK real grid functions. Declaration names and metadata are:

- `grmhd_velocities` (`vx`, `vy`, `vz`) and scalar `u0`:
  `InterpNumTimelevels=1`, no prolongation; `u0` also disables checkpointing.
- `grmhd_conservatives` (`rho_star`, `tau`, three `Stilde` components),
  `ent_star`, and `Ye_star`: three timelevels, ENO prolongation.
- scalar `Ax`, `Ay`, `Az`, `phitilde`: three timelevels with `STAGGER011`,
  `STAGGER101`, `STAGGER110`, and `STAGGER111` prolongation tags respectively.
- `grmhd_B_stagger` (`Bx_stagger`, `By_stagger`, `Bz_stagger`) and
  `grmhd_B_center` (`Bx_center`, `By_center`, `Bz_center`):
  `InterpNumTimelevels=1`, no prolongation.
- `grmhd_conservatives_rhs` (`rho_star_rhs`, `tau_rhs`, `Stildex_rhs`,
  `Stildey_rhs`, `Stildez_rhs`), scalar `ent_star_rhs`, scalar `Ye_star_rhs`,
  and `EM_rhs` (`Ax_rhs`, `Ay_rhs`, `Az_rhs`, `phitilde_rhs`): no prolongation
  and no checkpoint.
- `grmhd_cmin_cmax_temps` (`cmin_x`, `cmax_x`, `cmin_y`, `cmax_y`, `cmin_z`,
  `cmax_z`); `grmhd_flux_temps` (`rho_star_flux`, `tau_flux`,
  `Stildex_flux`, `Stildey_flux`, `Stildez_flux`); scalar `ent_star_flux` and
  `Ye_star_flux`; and `grmhd_primitives_reconstructed_temps` (`vxr`, `vyr`,
  `vzr`, `vxl`, `vyl`, `vzl`, right/left staggered B components, and the 12
  double-reconstructed velocity components): no prolongation and no checkpoint.
- diagnostic `failure_checker`: no prolongation/checkpoint and one interpolated
  timelevel; its description warns of overwrite at every RK substep.
- deprecated `em_psi6phi` (`psi6phi`) has three timelevels and `STAGGER111`;
  deprecated `grmhd_primitives_allbutBi` (`rho_b`, `P`) and
  `grmhd_primitives_Bi` (`Bx`, `By`, `Bz`) have
  `InterpNumTimelevels=1` and no prolongation.

This page owns names and tags only. Scientific state roles and conditional
extras belong to [State and EOS Modes](../evolution/state-and-eos-modes.md);
magnetic placement belongs to [Staggered State and Magnetic Reconstruction](../magnetics/staggered-state-and-magnetic-reconstruction.md).

### Build Inputs

Common `src/make.code.defn` selects `Hybrid`, `HybridEntropy`, `Tabulated`, and
`TabulatedEntropy` as `SUBDIRS`. Its 18 `SRCS` entries are `A_flux_rhs.c`,
`A_i_outer_boundaries.c`, `compute_B_and_Bstagger_from_A.c`,
`compute_metric_derivs.c`, `compute_Tmunu.c`, both `convert_*HydroBase*.c`
files, `evaluate_phitilde_and_A_gauge_rhs.c`, `InitSymBound.c`,
`interpolate_metric_to_face.c`, `MoL_registration.c`, `reconstruction_loop.c`,
`set_gz_symmetries.c`, `symmetry_gzs_staggered.c`, both
`backward_compatible_*.c` files, `specify_driver_BCs.c`, and `sync.c`. Each
variant's `make.code.defn` lists `calculate_fluxes_rhs.c`,
`conservs_to_prims.c`, `evaluate_fluxes_rhs.c`, `evaluate_sources_rhs.c`,
`hydro_outer_boundaries.c`, `perturb_conservatives.c`,
`perturb_primitives.c`, and `prims_to_conservs.c`.

Claim evidence:
- Claim: Five checked-in build definitions include common sources plus four eight-file variant families; static inclusion does not prove compilation, linking, or transitive completeness.
- Role: descriptive behavior
- Deciding authority: registered `IllinoisGRMHD/src/make.code.defn`, `SUBDIRS` and `SRCS`; four registered variant `make.code.defn` files, `SRCS`
- Corroboration: registered source aggregates contain every named file; no build was run
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=all five make.code.defn files; date=07-17-2026`

### Local Header Surface

`IllinoisGRMHD.h` includes Cactus argument/parameter headers and `GRHayLib.h`.
It defines reconstruction indices, a perturbation macro, interpolation and
derivative coefficients/macros, and prototypes for metric face interpolation,
metric derivatives, staggered symmetry ghost filling, PPM reconstruction loop,
and A-flux RHS. These declarations establish local helper interfaces, not CCTK-
generated symbol availability or external implementation semantics.

## Sources

- [`configuration.ccl`](../../IllinoisGRMHD/configuration.ccl) — `requires`.
- [`interface.ccl`](../../IllinoisGRMHD/interface.ccl) — thorn, include, group,
  tag/timelevel, deprecated-group, and function-alias declarations.
- [`src/make.code.defn`](../../IllinoisGRMHD/src/make.code.defn) — common
  `SUBDIRS` and `SRCS`.
- [`Hybrid/make.code.defn`](../../IllinoisGRMHD/src/Hybrid/make.code.defn),
  [`HybridEntropy/make.code.defn`](../../IllinoisGRMHD/src/HybridEntropy/make.code.defn),
  [`Tabulated/make.code.defn`](../../IllinoisGRMHD/src/Tabulated/make.code.defn),
  and [`TabulatedEntropy/make.code.defn`](../../IllinoisGRMHD/src/TabulatedEntropy/make.code.defn)
  — variant `SRCS`.
- [`IllinoisGRMHD.h`](../../IllinoisGRMHD/src/IllinoisGRMHD.h) — macros, enum,
  includes, and helper prototypes.

## See Also

- Parent: [Architecture](index.md)
- See also: [Schedule Lifecycle](schedule-lifecycle.md)
- Depends on: [State and EOS Modes](../evolution/state-and-eos-modes.md)
- See also: [Parameters and Runtime Controls](../integration/parameters-and-runtime-controls.md)
