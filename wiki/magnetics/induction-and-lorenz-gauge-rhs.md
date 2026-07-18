# Induction and Lorenz-Gauge RHS

> Status: confirmed · Last reconciled: 07-17-2026
> Up: [Magnetics](index.md)

## Summary

Hydrodynamic flux callers build twice-reconstructed velocity/B states and
characteristic speeds, then `IllinoisGRMHD_A_flux_rhs` computes HLL induction
terms for A. A later declared RHS stage adds gauge gradients and replaces the
`phitilde` RHS using interpolated metric/A data and GRHayL's configured Lorenz
damping input.

Claim evidence:

- Claim: Local callers and helpers establish the A-flux direction order, HLL input assembly, gauge-gradient addition, scratch reuse, and exact Lorenz damping input; external HLL, characteristic, interpolation, and gauge internals are not claimed.
- Role: descriptive behavior
- Deciding authority: registered `IllinoisGRMHD/src/A_flux_rhs.c` and `evaluate_phitilde_and_A_gauge_rhs.c` named functions
- Corroboration: registered four family `evaluate_fluxes_rhs.c` callers and `IllinoisGRMHD/schedule.ccl` RHS ordering
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-run; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-run; options=all four EOS callers; date=07-17-2026`

## Detail

### Caller and direction relationship

Every family `*_evaluate_fluxes_rhs` has the same magnetic call structure.
After each directional hydro reconstruction it uses the shared reconstruction
helper for transverse staggered B and already reconstructed velocities:

- after x and y hydro work, call A-flux direction 2 for `Az_rhs`;
- after y and z hydro work, call A-flux direction 0 for `Ax_rhs`;
- after z plus final transverse reconstruction, call A-flux direction 1 for
  `Ay_rhs`.

This interleaving reuses the family caller's `cmin/cmax` arrays and
`grmhd_primitives_reconstructed_temps`. Family EOS differences affect the
hydrodynamic states/speeds upstream; the local A-flux helper itself is shared.

### HLL induction contribution

`IllinoisGRMHD_A_flux_rhs` maps the requested A direction to the other two
velocity and B components. At every interior point it assembles
`ghl_HLL_vars` from four twice-reconstructed states for each transverse
velocity, two reconstructed states for each transverse densitized B, and
minimum/maximum speeds from the two transverse directions. It writes the
result of `ghl_HLL_flux_with_Btilde` into the selected A RHS.

The source explicitly notes that each speed is sampled at a slightly
different face location than the edge-centered A flux and attributes an
expected negligible effect to prior discussion. This is a source-commented
placement caveat, not locally executed validation. The external HLL formula
and characteristic-speed internals are not inferred.

### Gauge and `phitilde` stage

`IllinoisGRMHD_evaluate_phitilde_and_A_gauge_rhs` is declared after the family
flux group. It intentionally aliases reconstruction scratch arrays:

- `vxr/vyr/vzr` for interpolated shift;
- `vxrr` for interpolated lapse;
- `vxll` for interpolated `alpha*Phi - beta^j*A_j`; and
- `vxl/vyl/vzl` for interpolated `sqrt(gamma) A^i` components.

An interpolation loop extends two points into ghost zones. At each point it
builds a 2-by-2-by-2 cell-centered metric stencil and 3-by-3-by-3 A stencils,
then calls `ghl_interpolate_with_cell_centered_ADM`. The local source comment
says interpolation, not reconstruction, is intentional for this gauge stage.

Over interior points, backward differences of the interpolated
`alpha*Phi - beta^j*A_j` value are added to each `A_i` RHS. Five-point shift
and `phitilde` stencils plus two-point `sqrt(gamma) A^i` data are then passed
to `ghl_calculate_phitilde_rhs`; its return value replaces
`phitilde_rhs[index]`. The damping input is exactly
`ghl_params->Lorenz_damping_factor`. Local code establishes this input and
stencil handoff, not GRHayL's internal gauge formula.

The source-stage family functions zero A and `phitilde` RHS arrays first;
family flux calls assign induction A terms; this later gauge stage adds A
gradients and sets `phitilde_rhs`. Schedule declarations preserve that order.

## Sources

- `IllinoisGRMHD/src/A_flux_rhs.c` — `IllinoisGRMHD_A_flux_rhs`, HLL variable
  assembly, placement caveat, and `ghl_HLL_flux_with_Btilde` call.
- `IllinoisGRMHD/src/evaluate_phitilde_and_A_gauge_rhs.c` —
  `IllinoisGRMHD_evaluate_phitilde_and_A_gauge_rhs`, interpolation, scratch
  aliases, gauge gradients, damping input, and `phitilde` RHS call.
- `IllinoisGRMHD/src/reconstruction_loop.c` —
  `IllinoisGRMHD_reconstruction_loop` used for magnetic transverse states.
- Registered aggregates `illinoisgrmhd-hybrid`,
  `illinoisgrmhd-hybrid-entropy`, `illinoisgrmhd-tabulated`, and
  `illinoisgrmhd-tabulated-entropy` — four `*_evaluate_fluxes_rhs` callers and
  identical direction ordering.
- `IllinoisGRMHD/src/IllinoisGRMHD.h` — reconstruction indices and A-flux
  helper signature.
- `IllinoisGRMHD/schedule.ccl` — family flux blocks and
  `IllinoisGRMHD_evaluate_phitilde_and_A_gauge_rhs` ordering.

## See Also

- Parent: [Magnetics](index.md)
- Depends on: [Reconstruction, Fluxes, and Sources](../evolution/reconstruction-fluxes-and-sources.md)
- Depends on: [Staggered State and Magnetic Reconstruction](staggered-state-and-magnetic-reconstruction.md)
- See also: [Schedule Lifecycle](../architecture/schedule-lifecycle.md)
