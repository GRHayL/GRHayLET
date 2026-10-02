# Reconstruction, Fluxes, and Sources

> Status: confirmed · Last reconciled: 10-02-2026
> Up: [Evolution](index.md)

## Summary

Declared RHS scheduling initializes sources before adding directional flux
divergences. Each family builds metric/curvature source data, then performs
x/y/z PPM reconstruction, face-metric interpolation, characteristic-speed and
HLLE calls, and conservative flux differences. Reconstructed data and speeds
also feed magnetic induction, whose algorithm belongs to Magnetics.

Claim evidence:

- Claim: Local schedule declares source-before-flux; all four family implementations then use directional reconstruction, local face/derivative helpers, GRHayL speed/HLLE call boundaries, and flux divergence without a convergence guarantee.
- Role: descriptive behavior
- Deciding authority: registered `IllinoisGRMHD/schedule.ccl` RHS groups and four family `evaluate_*_rhs.c`/`calculate_fluxes_rhs.c` implementations
- Corroboration: registered `IllinoisGRMHD/src/reconstruction_loop.c`, `interpolate_metric_to_face.c`, `compute_metric_derivs.c`, and `IllinoisGRMHD.h` helper roles
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-run; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-run; options=Hybrid,HybridEntropy,Tabulated,TabulatedEntropy; date=07-17-2026`

## Detail

### Declared stage order

`IllinoisGRMHD_RHS` is scheduled in `MoL_CalcRHS`.
`IllinoisGRMHD_evaluate_sources_rhs` has no internal predecessor, and
`IllinoisGRMHD_evaluate_fluxes_rhs` is declared after it. Thus
source-before-flux is a schedule declaration, not an observed runtime trace.
All four family source functions set `rho_star_rhs` and enabled advected-extra
RHS values to zero, also zero `phitilde_rhs` and all three A RHS arrays, then
write tau and momentum source values. Flux functions subsequently add
directional differences.

### Metric and source stage

Over interior points, each `*_evaluate_sources_rhs` function initializes ADM
metric and extrinsic-curvature records and assembles family primitives. It
limits velocity/computes `u0`, then calls
`ghl_calculate_source_terms` with three metric-derivative records.

`IllinoisGRMHD_compute_metric_derivs` uses the header's four-point centered
`COMPUTE_DERIV` stencil on lapse, shift, and six spatial-metric components in
each coordinate direction. The local code supplies reciprocal grid spacing.
These implementation facts do not establish a measured convergence order.

### Directional flux stage

Each family `*_evaluate_fluxes_rhs` invokes its directional calculator for x,
y, and z. For a direction, `*_calculate_flux_dir_rhs`:

1. chooses the directional velocity, normal staggered-B component, component
   permutation, characteristic-speed function, and family/direction HLLE
   function;
2. reconstructs all velocity components with six-point PPM stencils, using
   pressure and normal velocity to compute flattening inputs;
3. interpolates lapse, shift, and metric to the face through
   `IllinoisGRMHD_interpolate_metric_to_face`, whose `COMPUTE_FCVAL` uses
   coefficients `-0.0625, 0.5625, 0.5625, -0.0625`;
4. reconstructs rho with the steepening form, plus pressure, two transverse
   centered-B components, and family extras; the normal B comes from the
   densitized staggered value divided by face `sqrt_detgamma`;
5. limits left/right velocity, calls the selected characteristic-speed and
   HLLE routines, and writes face fluxes; and
6. adds `(flux[index] - flux[indp1]) / delta` to rho, tau, momentum, and
   enabled entropy/electron-fraction RHS arrays.

Hybrid families calculate an effective Gamma locally for density steepening.
Tabulated families pass `1.0`; reconstruct `Y_e`, and optionally entropy; seed
face temperature from the current point; enforce rho/`Y_e`/pressure bounds;
then locally call the matching pressure-to-thermodynamics routine. Entropy
families carry entropy flux, and tabulated families carry electron-fraction
flux.

### Shared reconstruction and induction handoff

`IllinoisGRMHD_reconstruction_loop` is EOS-independent by local design: it is
used only for staggered B and already reconstructed velocities. The family
flux callers reuse it for the second, transverse reconstruction required by
the A RHS, retain `cmin/cmax`, and invoke `IllinoisGRMHD_A_flux_rhs` between
directional hydro stages. Details belong to
[Induction and Lorenz-Gauge RHS](../magnetics/induction-and-lorenz-gauge-rhs.md).

This page makes no untested claim about convergence, thread safety, numerical
stability, or external GRHayL algorithms.

### Defined derivative carriers and checked face thermodynamics

`IllinoisGRMHD_compute_metric_derivs` initializes a raw derivative carrier with
lapse, shift, and symmetric covariant spatial-metric derivatives. It never
passes that tensor to physical-metric inversion. Both tabulated face kernels
check each pressure-to-thermal inverse return status before characteristic
speed/flux calls and abort through the existing error policy on failure.
Hybrid entropy reconstruction reads the thorn proxy; the compiled tabulated
entropy kernel remains guarded out at startup.

Claim evidence:
- Claim: Derivative preparation avoids physical-metric inversion and both table face inversions check errors before flux evaluation.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/src/compute_metric_derivs.c`, `IllinoisGRMHD_compute_metric_derivs`; registered `IllinoisGRMHD/src/Tabulated/calculate_fluxes_rhs.c`, `IllinoisGRMHD_tabulated_calculate_flux_dir_rhs`; registered `IllinoisGRMHD/src/TabulatedEntropy/calculate_fluxes_rhs.c`, `IllinoisGRMHD_tabulated_entropy_calculate_flux_dir_rhs`
- Corroboration: none available for the local status-check policy; the two registered face kernels independently implement the same check, and the derivative helper contains the raw assignments
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=inspected-not-run; options=local source and declaration inspection; date=10-02-2026`

## Sources

- `IllinoisGRMHD/schedule.ccl` — groups `IllinoisGRMHD_RHS`,
  `IllinoisGRMHD_evaluate_sources_rhs`, and
  `IllinoisGRMHD_evaluate_fluxes_rhs`.
- Registered aggregates `illinoisgrmhd-hybrid`,
  `illinoisgrmhd-hybrid-entropy`, `illinoisgrmhd-tabulated`, and
  `illinoisgrmhd-tabulated-entropy` — four `*_evaluate_sources_rhs`,
  `*_evaluate_fluxes_rhs`, and `*_calculate_flux_dir_rhs` function sets.
- `IllinoisGRMHD/src/reconstruction_loop.c` —
  `IllinoisGRMHD_reconstruction_loop`.
- `IllinoisGRMHD/src/interpolate_metric_to_face.c` —
  `IllinoisGRMHD_interpolate_metric_to_face`.
- `IllinoisGRMHD/src/compute_metric_derivs.c` —
  `IllinoisGRMHD_compute_metric_derivs`.
- `IllinoisGRMHD/src/IllinoisGRMHD.h` — `COMPUTE_FCVAL` and
  `COMPUTE_DERIV` macros.

## See Also

- Parent: [Evolution](index.md)
- Depends on: [State and EOS Modes](state-and-eos-modes.md)
- Implements: [Schedule Lifecycle](../architecture/schedule-lifecycle.md)
- See also: [Induction and Lorenz-Gauge RHS](../magnetics/induction-and-lorenz-gauge-rhs.md)
