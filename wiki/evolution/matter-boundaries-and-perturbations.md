# Matter Boundaries and Perturbations

> Status: confirmed · Last reconciled: 07-17-2026
> Up: [Evolution](index.md)

## Summary

`Matter_BC` selects copy, outflow, or frozen primitive boundaries. Approximate
fills run only after iteration zero on refinement level zero, require equal
ghost widths, and rebuild limited primitives plus conservatives per face.
Separate gates perturb initial primitives or conservatives before every
recovery.

Claim evidence:

- Claim: Local matter boundary code implements the stated frozen/iteration/refinement/equal-ghost guards, copy/outflow face fills, and family limit/recompute steps; perturbation gates and fields are schedule/code declarations, not observed execution.
- Role: public/scientific contract
- Deciding authority: registered aggregates `illinoisgrmhd-hybrid`, `illinoisgrmhd-hybrid-entropy`, `illinoisgrmhd-tabulated`, and `illinoisgrmhd-tabulated-entropy`; four `*_hydro_outer_boundaries` and eight perturbation functions
- Corroboration: registered `IllinoisGRMHD/param.ccl` controls, `schedule.ccl` groups, and `src/InitSymBound.c` shared preconditions
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=inspected-not-run; options=copy,outflow,frozen and four EOS families; date=07-17-2026`

## Detail

### Shared setup constraints

At `BASEGRID`, `IllinoisGRMHD_InitSymBound` enforces two preconditions:

- `Matter_BC=frozen` and `EM_BC=frozen` must be selected together; selecting
  only one is a fatal error.
- Every dimension must provide at least three ghost zones for this thorn's PPM
  implementation.

EM behavior routes to
[Electromagnetic Boundaries and Symmetry](../magnetics/electromagnetic-boundaries-and-symmetry.md).
Parameter declarations and defaults route to
[Parameters and Runtime Controls](../integration/parameters-and-runtime-controls.md).

### Matter boundary control and ordering

`IllinoisGRMHD_hydro_outer_boundaries` is declared after Con2Prim. Each of the
four family functions follows the same guards and face order:

1. `Matter_BC=frozen` returns immediately.
2. Iteration zero or any nonzero refinement level returns before approximate
   filling. Therefore filling continues only after iteration zero on the
   coarsest level.
3. Unequal x/y/z ghost widths cause a fatal error.
4. For each ghost layer, code processes x-max, x-min, y-max, y-min, z-max,
   then z-min. Only faces marked by `cctk_bbox` are filled. z-min is skipped
   unless `Symmetry=none`.

`copy` copies rho, pressure, all velocity components, and family extras from
the adjacent inward point. `outflow` begins with the same copy but zeros only
the face-normal velocity when it points inward: negative at a maximum face,
positive at a minimum face. This is done independently on each axis.

At every filled point, centered B is taken from the destination point. The
family helper initializes/enforces the local metric, applies primitive limits
and computes `u0`, computes conservatives, then writes corrected primitives
and `rho_star`, `tau`, and `Stilde`. HybridEntropy also copies/writes entropy
and `ent_star`; Tabulated copies/writes `Y_e`, temperature, and `Ye_star`;
TabulatedEntropy does both. Boundary-fill limiting and recomputation are owned
here, not by ordinary conversion.

### Driver registration

`IllinoisGRMHD_specify_driver_BCs` registers `none` for HydroBase rho,
pressure, internal energy, entropy, `Y_e`, and temperature, plus IllinoisGRMHD
velocity, conservative, `ent_star`, and `Ye_star` fields. Its source comment
says IllinoisGRMHD handles these manually. Registration is unconditional in
this file even when a field is unused by the active family.

### Perturbation gates and fields

If `perturb_initial_data` is true, the selected family primitive perturbation
is declared after HydroBase ingress and before Prim2Con. Each full-grid loop
calls `srand(random_seed)` once and multiplies each selected value by
`one_plus_pert(random_pert)`, defined locally as
`1 + random_pert * rand() / RAND_MAX`.

All primitive variants perturb rho, pressure, three velocities, `phitilde`,
and `Ax/Ay/Az`; HybridEntropy adds entropy, Tabulated adds `Y_e` and
temperature, and TabulatedEntropy adds all three. The code does not directly
perturb internal energy.

If `perturb_every_con2prim` is true, the selected conservative perturbation is
declared after B reconstruction and before Con2Prim. All variants perturb
`rho_star`, `tau`, and all `Stilde` components; entropy variants add
`ent_star`, tabulated variants add `Ye_star`, and the combined variant adds
both. These routines expose no local perturbation counter or result log.
File presence and schedule declarations do not prove fixture execution.

## Sources

- `IllinoisGRMHD/param.ccl` — `Matter_BC`, `random_seed`, `random_pert`,
  `perturb_initial_data`, and `perturb_every_con2prim` declarations.
- `IllinoisGRMHD/schedule.ccl` — groups
  `IllinoisGRMHD_hydro_outer_boundaries`,
  `IllinoisGRMHD_perturb_primitives`, and
  `IllinoisGRMHD_perturb_conservatives`.
- Registered aggregates `illinoisgrmhd-hybrid`,
  `illinoisgrmhd-hybrid-entropy`, `illinoisgrmhd-tabulated`, and
  `illinoisgrmhd-tabulated-entropy` — four `*_hydro_outer_boundaries`, four
  `*_perturb_primitives`, and four `*_perturb_conservatives` functions.
- `IllinoisGRMHD/src/IllinoisGRMHD.h` — `one_plus_pert` macro.
- `IllinoisGRMHD/src/InitSymBound.c` —
  `IllinoisGRMHD_InitSymBound` frozen pairing and PPM ghost-zone checks.
- `IllinoisGRMHD/src/specify_driver_BCs.c` —
  `IllinoisGRMHD_specify_driver_BCs` matter-side registrations.

## See Also

- Parent: [Evolution](index.md)
- Depends on: [Con2Prim Recovery and Diagnostics](con2prim-recovery-and-diagnostics.md)
- See also: [Electromagnetic Boundaries and Symmetry](../magnetics/electromagnetic-boundaries-and-symmetry.md)
- See also: [Parameters and Runtime Controls](../integration/parameters-and-runtime-controls.md)
