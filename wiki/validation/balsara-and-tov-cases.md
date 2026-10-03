# Balsara and TOV Cases

> Status: confirmed · Last reconciled: 10-02-2026
> Up: [Validation](index.md)

## Summary

Six example configurations are visible: five named Balsara 1D initial-data
cases and one magnetized TOV case. Test wrappers exist for all six; active
`test.ccl` rows and oracle directories cover Balsara1, 2, 3, and 5.
`magnetizedTOV` has files/oracles but no case-specific declaration. This is configured/file
coverage only, not executed validation.

Claim evidence:

- Claim: Visible case files explicitly select only Simple and Hybrid EOS,
  never set `evolve_entropy`, and therefore do not establish effective entropy
  family or Tabulated execution.
- Role: generated evidence
- Deciding authority: `IllinoisGRMHD/par/*.par` and
  `IllinoisGRMHD/test/**/*.par` fixture roles and selector assignments
- Corroboration: `IllinoisGRMHD/schedule.ccl` EOS/entropy variant gates;
  `IllinoisGRMHD/test/test.ccl` supplies only configured test declarations
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-run; tool_version=not-run; backend=not-run; precision=not-run; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=static case comparison; date=07-17-2026`

## Detail

### Shared Balsara shape

All five example and test-wrapper Balsara files activate IllinoisGRMHD,
GRHayLib, and GRHayLID; request `HydroTest1D`, GRHayLID A/B initial data,
Cartesian Minkowski ADM data, three ghost/boundary points on every face, and
staggered Carpet operators. Each selects `GRHayLID::initial_data_1D` matching
case name, x shock direction, and discontinuity position `0.000001`.
Underlying initial states live outside this tree and are not inferred.

All explicitly set:

- `GRHayLib::EOS_type = Simple`;
- `IllinoisGRMHD::rescale_magnetics = no`;
- `IllinoisGRMHD::Matter_BC = copy`;
- `IllinoisGRMHD::update_Tmunu = no`.

Balsara1 sets `GRHayLib::Gamma=2`; Balsara2–5 set `5/3`. No Balsara file
explicitly sets `EM_BC`, `Convert_to_HydroBase_every`, or `evolve_entropy`.

### Balsara example/test matrix

| Case | `test.ccl` / oracle state | `par/` example grid and run | `test/` wrapper grid and run | Visible output |
| --- | --- | --- | --- | --- |
| Balsara1 | active, `NPROCS 1`; 8 profile oracles | bounds `[-0.5,0.5]^3`; spacing `(1/1600,1/8,1/8)`; RK4; final time `0.4`; output every 1280 | x `[-0.5,0.5]`, y/z `[-0.0025,0.0025]`; cells `20x8x8`; RK2; 10 iterations; output every 10 | `rho`, `press`, `vx/vy/vz`, centered `Bx/By/Bz` along x |
| Balsara2 | active, `NPROCS 1`; 8 profile oracles | same as Balsara1 | same as Balsara1 | same eight profiles |
| Balsara3 | active, `NPROCS 1`; 8 profile oracles | same as Balsara1 | same thin grid; RK2; 2 iterations; output every 2 | same eight profiles |
| Balsara4 | `TEST` block commented; no oracle directory | same as Balsara1; additionally sets `GRHayLib::max_Lorentz_factor=25` | x `[-0.5,0.5]`, y/z `[-0.005,0.005]`; cells `20x16x16`; RK2; 2 iterations; output every 2; same Lorentz setting | configured same eight profiles, but no checked-in `.asc` evidence |
| Balsara5 | active, `NPROCS 1`; 8 profile oracles | x `[-0.5,0.8]`, y/z `[-0.5,0.5]`; same spacing; RK4; final time `0.55`; output every 1760 | same thin grid as Balsara1; RK2; 10 iterations; output every 10 | same eight profiles |

Checked-in nested Balsara fixtures reflect test-wrapper-scale configurations,
not longer `par/` examples. Their headers identify Cactus 4.15.0 generation
on June 24, 2024. Date/header presence does not prove present reproducibility.

### Magnetized TOV shape

Both TOV files activate IllinoisGRMHD, GRHayLib, Seed_Magnetic_Fields,
TOVSolver, ML_BSSN, and related infrastructure. They configure a
`[-12,12]^3` domain with spacing one, three ghost/boundary points, four Carpet
refinement levels, RK4 with four intermediate steps, TOV density/EOS data,
and `GRHayLib::EOS_type=Hybrid`. They enable IllinoisGRMHD staggered A fields
with `A_b=0.64428596382321`. Example additionally writes
`Afield_type="Pressure_prescription"`; test wrapper does not explicitly set
that external parameter.

| Aspect | `par/magnetizedTOV.par` | `test/magnetizedTOV.par` and fixture |
| --- | --- | --- |
| Iterations | `cctk_itlast=4` | `cctk_itlast=2` |
| TmunuBase | storage and RHS enabled, one timelevel | same |
| Illinois-owned overrides | none explicit | none explicit |
| Scalar output | current HydroBase `rho`; `maximum` reduction every 2 | HydroBase `rho` and centered B; `maximum minimum` every 2 |
| Oracle files | none under `par/` | min/max for `rho`, `Bx_center`, `By_center`, `Bz_center` |

`test.ccl` has no block for `magnetizedTOV`; the case uses the thorn-wide
tolerances and has no test-specific `NPROCS`. No harness run is inferred.

### Selector and coverage limits

Visible case files explicitly select only `Simple` (Balsara) and `Hybrid`
(TOV). None selects `Tabulated`. None explicitly sets `evolve_entropy`.
Because shared parameter's default belongs to external GRHayLib declaration,
effective non-entropy versus entropy family—thus Hybrid versus HybridEntropy
execution—is unknown for these fixtures. No Tabulated, entropy-family,
convergence, AMR-correctness, restart, divergence-control, or external-library
execution is established.

`test.ccl` global comparison tolerance is `1e-11` absolute and relative, but
only a recorded run/comparison could show whether any case currently meets it.

### Modern TOV output storage

The modern example requests HydroBase rho/press/eps and current IllinoisGRMHD
velocities. It does not need deprecated primitive storage supplied only by
`ID_converter_ILGRMHD`. The existing regression parfile is unchanged.

Claim evidence:
- Claim: Modern magnetizedTOV output lists use current populated primitive groups rather than legacy density storage.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/par/magnetizedTOV.par, IO output lists`
- Corroboration: registered `IllinoisGRMHD/schedule.ccl, affected declarations`
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=inspected-not-run; options=local source and declaration inspection; date=10-02-2026`

## Sources

- [`IllinoisGRMHD/par/`](../../IllinoisGRMHD/par/) — six shipped example
  parameter files; stable fixture roles Balsara1–5 and `magnetizedTOV`.
- [`IllinoisGRMHD/test/test.ccl`](../../IllinoisGRMHD/test/test.ccl) — active
  Balsara declarations, tolerances, process counts, and Balsara4 rationale.
- [`IllinoisGRMHD/test/Balsara1.par`](../../IllinoisGRMHD/test/Balsara1.par),
  [`Balsara2.par`](../../IllinoisGRMHD/test/Balsara2.par),
  [`Balsara3.par`](../../IllinoisGRMHD/test/Balsara3.par),
  [`Balsara4.par`](../../IllinoisGRMHD/test/Balsara4.par), and
  [`Balsara5.par`](../../IllinoisGRMHD/test/Balsara5.par) — top-level test wrappers.
- [`IllinoisGRMHD/test/magnetizedTOV.par`](../../IllinoisGRMHD/test/magnetizedTOV.par) —
  top-level TOV test wrapper.
- [`IllinoisGRMHD/test/Balsara1/`](../../IllinoisGRMHD/test/Balsara1/) —
  representative checked-in Balsara per-case fixture and eight profile oracles;
  Balsara2, 3, and 5 have same filename roles.
- [`IllinoisGRMHD/test/magnetizedTOV/`](../../IllinoisGRMHD/test/magnetizedTOV/) —
  checked-in per-case fixture and min/max scalar oracles.

## See Also

- Parent: [Validation](index.md)
- Depends on: [Test Harness and Oracles](test-harness-and-oracles.md)
- Contrasts with: [State and EOS Modes](../evolution/state-and-eos-modes.md)
- See also: [Staggered State and Magnetic Reconstruction](../magnetics/staggered-state-and-magnetic-reconstruction.md)
