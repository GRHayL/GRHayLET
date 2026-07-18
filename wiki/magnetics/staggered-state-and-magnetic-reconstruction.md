# Staggered State and Magnetic Reconstruction

> Status: confirmed · Last reconciled: 07-17-2026
> Up: [Magnetics](index.md)

## Summary

IllinoisGRMHD evolves semi-staggered `Ax/Ay/Az` and fully staggered
`phitilde`, derives densitized face-centered B from a discrete curl, then
averages and undensitizes it into vertex-centered B. This reconstruction is
declared during initial conversion and before each recurring Con2Prim.

Claim evidence:

- Claim: Local declarations and curl code establish field placement plus densitized-staggered to undensitized-centered reconstruction; README divergence/AMR language remains attributed intent, not a tested guarantee.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/interface.ccl` magnetic groups and `IllinoisGRMHD/src/compute_B_and_Bstagger_from_A.c` named function
- Corroboration: registered `IllinoisGRMHD/schedule.ccl` initial/recurring blocks and `IllinoisGRMHD/README` design statement
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-run; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-run; options=Symmetry none inspected; date=07-17-2026`

## Detail

### Declared placement and storage

| Field | Logical placement declared by `interface.ccl` | Timelevels / tags |
| --- | --- | --- |
| `Ax` | `(i, j+1/2, k+1/2)` | 3; `Prolongation="STAGGER011"` |
| `Ay` | `(i+1/2, j, k+1/2)` | 3; `Prolongation="STAGGER101"` |
| `Az` | `(i+1/2, j+1/2, k)` | 3; `Prolongation="STAGGER110"` |
| `phitilde = sqrt(gamma) Phi` | `(i+1/2, j+1/2, k+1/2)` | 3; `Prolongation="STAGGER111"` |
| `Bx/By/Bz_stagger` | `(i+1/2,j,k)`, `(i,j+1/2,k)`, `(i,j,k+1/2)` | one interpolated timelevel; prolongation `none` |
| `Bx/By/Bz_center` | vertices | one interpolated timelevel; prolongation `none` |

The interface warns that A requires semi-staggered prolongation/restriction
and `phitilde` requires fully staggered prolongation/restriction. It gives B
staggered fields as densitized and B centered fields as ordinary vertex
components. These declarations do not prove a particular Driver's support.

The README states repository design intent: evolving staggered vector
potential instead of B directly is intended to keep B divergenceless,
including across AMR boundaries. Initial static inspection neither tests nor
upgrades that statement into an AMR/divergence guarantee.

### `IllinoisGRMHD_compute_B_and_Bstagger_from_A`

The routine first calls the manual staggered symmetry helper for each A
component and `phitilde`, with component-specific parity arrays and staggering
arguments. Under the publicly selectable `Symmetry=none`, the helper returns
without filling.

Its first full-grid loop computes a discrete curl:

- `Bx_stagger = d_y Az - d_z Ay`;
- `By_stagger = d_z Ax - d_x Az`;
- `Bz_stagger = d_x Ay - d_y Ax`.

At lower index edges, shifted indices clamp into a copy pattern; the same
index construction supplies the source's stated upper-boundary copy handling.
No metric determinant division occurs here, so staggered B remains
densitized.

The second full-grid loop evaluates `sqrt(abs(det(gamma_ij)))` directly from
the six metric components. Each centered component is the average of the two
neighboring staggered values along its own axis, divided by that local square
root determinant. Thus centering and undensitization occur together. The
routine finishes by invoking the symmetry helper for centered and staggered B
with B-component parity arrays and component staggering.

Magnetic normalization/rescaling during HydroBase exchange is separate; see
[HydroBase, GRHayLib, and Tmunu](../integration/hydrobase-grhaylib-and-tmunu.md).

### Declared invocation points

`schedule.ccl` places this routine:

- in `IllinoisGRMHD_Prim2Con2Prim`, after HydroBase ingress and before
  Prim2Con; and
- in recurring `IllinoisGRMHD_Con2Prim`, after A outer-boundary handling and
  before conservative perturbation/recovery.

Both schedule blocks declare writes to A, `phitilde`, centered B, and
staggered B, then `SYNC` both B groups. This is a declared schedule/call
contract, not proof of runtime execution.

## Sources

- `IllinoisGRMHD/interface.ccl` — groups `Ax`, `Ay`, `Az`, `phitilde`,
  `grmhd_B_stagger`, and `grmhd_B_center`, including placement/tags.
- `IllinoisGRMHD/src/compute_B_and_Bstagger_from_A.c` —
  `IllinoisGRMHD_compute_B_and_Bstagger_from_A`, discrete curl, centering,
  determinant conversion, and symmetry calls.
- `IllinoisGRMHD/src/IllinoisGRMHD.h` — declaration of
  `IllinoisGRMHD_set_symmetry_gzs_staggered` and reconstruction indices.
- `IllinoisGRMHD/src/symmetry_gzs_staggered.c` — reflected-index formula used
  by magnetic reconstruction.
- `IllinoisGRMHD/schedule.ccl` — initial and recurring reconstruction blocks.
- `IllinoisGRMHD/README` — documented vector-potential design claim.

## See Also

- Parent: [Magnetics](index.md)
- Depends on: [Cactus Surface and Build](../architecture/cactus-surface-and-build.md)
- See also: [Induction and Lorenz-Gauge RHS](induction-and-lorenz-gauge-rhs.md)
- See also: [Electromagnetic Boundaries and Symmetry](electromagnetic-boundaries-and-symmetry.md)
