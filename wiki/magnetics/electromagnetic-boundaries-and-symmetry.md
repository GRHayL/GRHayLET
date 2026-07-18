# Electromagnetic Boundaries and Symmetry

> Status: confirmed · Last reconciled: 07-17-2026
> Up: [Magnetics](index.md)

## Summary

`EM_BC=frozen` skips A/`phitilde` outer filling; public `copy` selection runs
linear extrapolation for those evolved potentials, after which magnetic
reconstruction supplies copy-style B edges. Public `Symmetry` accepts only
`none`. Dormant equatorial branches exist but do not constitute supported
equatorial operation.

Claim evidence:

- Claim: Public `Symmetry` accepts only `none`; forced `equatorial` enters `IllinoisGRMHD_InitSymBound` parity setup rather than its impossible-string error, while a separate staggered helper explicitly aborts equatorial use. This is not supported equatorial operation.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/param.ccl` `Symmetry`; `IllinoisGRMHD/src/InitSymBound.c` branch; `symmetry_gzs_staggered.c` equatorial guard
- Corroboration: registered `IllinoisGRMHD/src/set_gz_symmetries.c` and `compute_B_and_Bstagger_from_A.c` manual parity call sites
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=inspected-not-run; options=none public,equatorial source-forced; date=07-17-2026`

## Detail

### EM outer boundaries

`IllinoisGRMHD_A_i_outer_boundaries` returns immediately for
`EM_BC=frozen`. Otherwise, including the public `copy` selection, it:

1. returns at iteration zero or on any nonzero refinement level;
2. rejects unequal x/y/z ghost widths;
3. processes every physical ghost layer in x-max, y-max, z-max, x-min, y-min,
   z-min order; and
4. sets `phitilde` and all three A components by linear extrapolation,
   `2*adjacent - next-adjacent`.

z-min filling runs only for `Symmetry=none`. `cctk_bbox` gates every physical
face. B itself is rebuilt later: `IllinoisGRMHD_compute_B_and_Bstagger_from_A`
uses clamped shifted indices for its documented lower/upper copy handling and
then applies centered averaging. Shared frozen pairing and minimum ghost-zone
requirements belong to
[Matter Boundaries and Perturbations](../evolution/matter-boundaries-and-perturbations.md).

### Driver and sync declarations

`IllinoisGRMHD_specify_driver_BCs` registers Driver BC `none` for `Ax`, `Ay`,
`Az`, `phitilde`, centered B, and staggered B because local routines handle
them manually. The function also ends with an unconditional `printf("3\n")`;
this is a code-visible setup output, not validation that registration worked.
The recurring Con2Prim schedule declares a `SYNC` of
conservatives, conditional extras, A, and `phitilde` through the empty C entry
point `IllinoisGRMHD_sync`; reconstruction blocks separately declare `SYNC`
for both B groups. These are CCL/Driver registrations and declarations, not an
observed communication trace.

### Public symmetry boundary and dormant code

`param.ccl` exposes only `Symmetry="none"` and describes equatorial support as
work in progress. If an out-of-schema caller nevertheless forces
`equatorial`, `IllinoisGRMHD_InitSymBound` enters its equatorial branch: it
assigns B-center z parities using `Sym_Bz` and odd z parity for `Stildez` and
`vz`. Its “impossible symmetry” error applies only to strings that are neither
`equatorial` nor `none`; forcing `equatorial` does not enter that error branch.

This does not make equatorial symmetry operational. The separate
`IllinoisGRMHD_set_symmetry_gzs_staggered` helper explicitly calls
`CCTK_VERROR` when `Symmetry=equatorial`, with a source warning that support is
unsafe and the error would need removal. When allowed past that guard, it
derives the number of negative-z ghost points and reflects with
`(2*num_gzs - stagger_z) - k`; only `stagger_z` changes that formula, while
its `stagger_x` and `stagger_y` arguments are source-marked unused.

Parity data are spread across local routines:

- `IllinoisGRMHD_InitSymBound` initially registers even group parity, then the
  forced-equatorial overrides described above.
- `IllinoisGRMHD_compute_B_and_Bstagger_from_A` supplies component parity
  arrays for A, `phitilde`, centered B, and staggered B around curl work.
- `IllinoisGRMHD_set_gz_symmetries`, scheduled after initial recovery, calls
  the staggered helper for `phitilde`, A, and centered/staggered B only inside
  its equatorial branch.

Therefore valid public input supports `none`; partial forced-equatorial code
and an explicit helper abort are a support limitation, not a contradiction.

## Sources

- `IllinoisGRMHD/param.ccl` — `EM_BC`, `Symmetry`, and `Sym_Bz` declarations.
- `IllinoisGRMHD/src/A_i_outer_boundaries.c` —
  `IllinoisGRMHD_A_i_outer_boundaries`, guards and linear extrapolation.
- `IllinoisGRMHD/src/compute_B_and_Bstagger_from_A.c` — copy-style B edge
  handling and magnetic parity arrays.
- `IllinoisGRMHD/src/specify_driver_BCs.c` —
  `IllinoisGRMHD_specify_driver_BCs` EM registrations.
- `IllinoisGRMHD/src/InitSymBound.c` — `IllinoisGRMHD_InitSymBound`, frozen
  pairing, public/dormant symmetry branches, and parity registration.
- `IllinoisGRMHD/src/set_gz_symmetries.c` —
  `IllinoisGRMHD_set_gz_symmetries` post-initial equatorial calls.
- `IllinoisGRMHD/src/symmetry_gzs_staggered.c` — equatorial abort and reflected
  ghost-index formula.
- `IllinoisGRMHD/src/sync.c` — empty `IllinoisGRMHD_sync` entry point.
- `IllinoisGRMHD/schedule.ccl` — Driver, sync, A-boundary, B-reconstruction,
  and post-initial symmetry blocks.

## See Also

- Parent: [Magnetics](index.md)
- Depends on: [Staggered State and Magnetic Reconstruction](staggered-state-and-magnetic-reconstruction.md)
- See also: [Matter Boundaries and Perturbations](../evolution/matter-boundaries-and-perturbations.md)
- See also: [Parameters and Runtime Controls](../integration/parameters-and-runtime-controls.md)
