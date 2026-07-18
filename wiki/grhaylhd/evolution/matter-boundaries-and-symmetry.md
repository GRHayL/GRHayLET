# Matter Boundaries and Symmetry

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Evolution](index.md)

## Scope and Non-Scope

This page records selectable matter-boundary declarations, visible six-face
boundary dataflow in all variants, and local symmetry declarations/code.
Boundary correctness, external symmetry APIs, and dormant equatorial support
are not established.

## Summary

`Matter_BC` permits `copy`, `outflow`, or `frozen`. Variant routines return for
frozen, otherwise copy adjacent primitives at six physical faces; outflow also
zeros a normal velocity whose sign points inward. Helpers then limit primitives
and recompute boundary conservatives. `Symmetry` permits only `none`, although
local code retains dormant equatorial branches and spelling mismatches.

## Variant Applicability

| Applicability | Copied primitive extras | Recomputed conservative extras |
| --- | --- | --- |
| Common | Density, pressure, three velocities; `BU` fixed zero | Five core conservatives and `u0` |
| Hybrid/Simple | None | None |
| Hybrid/Simple+Entropy | Entropy | `ent_star` |
| Tabulated | `Y_e`, temperature | `Ye_star` |
| Tabulated+Entropy | Entropy, `Y_e`, temperature | `ent_star`, `Ye_star` |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `EV-BC-01` | Parameter permits copy, outflow, and frozen matter boundaries. | declared | Boundary parameter | `ccl:GRHayLHD/param.ccl#parameter=Matter_BC` |
| `EV-BC-02` | Hybrid boundary visibly handles six faces and recomputes core state. | visible-implementation | Variant boundary function | `c:GRHayLHD/src/Hybrid/outer_boundaries.c#symbol=GRHayLHD_hybrid_outer_boundaries` |
| `EV-BC-03` | HybridEntropy boundary additionally maps entropy. | visible-implementation | Variant boundary function | `c:GRHayLHD/src/HybridEntropy/outer_boundaries.c#symbol=GRHayLHD_hybrid_entropy_outer_boundaries` |
| `EV-BC-04` | Tabulated boundary additionally maps electron fraction and temperature. | visible-implementation | Variant boundary function | `c:GRHayLHD/src/Tabulated/outer_boundaries.c#symbol=GRHayLHD_tabulated_outer_boundaries` |
| `EV-BC-05` | TabulatedEntropy boundary maps both optional state sets. | visible-implementation | Variant boundary function | `c:GRHayLHD/src/TabulatedEntropy/outer_boundaries.c#symbol=GRHayLHD_tabulated_entropy_outer_boundaries` |
| `EV-SYM-01` | Parameter permits only `none` symmetry. | declared | Symmetry parameter | `ccl:GRHayLHD/param.ccl#parameter=Symmetry` |
| `EV-SYM-02` | Initialization visibly contains a dormant equatorial branch and calls an error API below three ghost zones. | visible-implementation | Symmetry initialization | `c:GRHayLHD/src/InitSymBound.c#symbol=GRHayLHD_InitSymBound` |

## Details

### Common gates

All four boundary routines return immediately when `Matter_BC` is `frozen`.
They set `do_outflow` only for `outflow`; therefore `copy` retains adjacent
normal velocity. They also return at iteration zero or when
`GetRefinementLevel(cctkGH) != 0`, require equal ghost-zone counts in all three
directions, and process each ghost layer. Work is gated by matching
`cctk_bbox` face flags. Lower-z work additionally requires local
`Symmetry_none`.

### Six-face copy and inflow signs

| Face | Adjacent source | Outflow clamp on normal component |
| --- | --- | --- |
| x maximum | one x point inward | clamp `vx < 0` to zero |
| x minimum | one x point inward | clamp `vx > 0` to zero |
| y maximum | one y point inward | clamp `vy < 0` to zero |
| y minimum | one y point inward | clamp `vy > 0` to zero |
| z maximum | one z point inward | clamp `vz < 0` to zero |
| z minimum | one z point inward | clamp `vz > 0` to zero; only when symmetry is none |

Tangential velocities and mode scalar primitives copy unchanged from adjacent
point. Each local primitive object has `BU[0..2]` set zero. Mode helper builds
metric and auxiliaries at boundary point, calls primitive-limit/`u0` helper,
computes conservatives, and writes limited primitives plus core and optional
conservatives.

### Declared sync boundary

Schedule CCL places concrete outer-boundary routine after recovery and declares
sync for primitive fields. It declares no conservative sync in these blocks,
despite visible recomputation. This is declaration inventory, not an assertion
about external synchronization behavior.

### Selectable symmetry and dormant code

`param.ccl` offers only `"none"`. `GRHayLHD_InitSymBound` calls `SetCartSymGN`
with positive Cartesian flags for core state and conditional optional groups,
and calls an error API when any ghost-zone count is below three. It contains
an `equatorial` branch that negates z momentum/velocity symmetry and an error
for any value other than `none`; equatorial cannot be selected through local
parameter range. That branch spells momentum `Stilde_z`, while interface
declares `Stildez`. Three recovery equatorial branches also request
`grhd_conservatives` rather than declared `grmhd_conservatives`.

All four boundary file headers describe vector-potential and magnetic-field
stages absent from local declared state/lifecycle. These comments do not expand
scope or prove magnetic behavior.

## Caveats

- Only `Symmetry=none` is locally selectable; see
  [GRH-0001](../contradictions.md#grh-0001).
- Dormant momentum spelling differs; see
  [GRH-0003](../contradictions.md#grh-0003).
- Recovery group spelling differs; see
  [GRH-0002](../contradictions.md#grh-0002).
- Boundary headers retain magnetic stages; see
  [GRH-0005](../contradictions.md#grh-0005).
- External `GetRefinementLevel`, symmetry, ghost-zone, and synchronization
  semantics remain unverified.

## Sources

- [Parameter declarations](../../../GRHayLHD/param.ccl)
- [Schedule declarations](../../../GRHayLHD/schedule.ccl)
- [Symmetry initialization](../../../GRHayLHD/src/InitSymBound.c)
- [Hybrid boundary](../../../GRHayLHD/src/Hybrid/outer_boundaries.c)
- [HybridEntropy boundary](../../../GRHayLHD/src/HybridEntropy/outer_boundaries.c)
- [Tabulated boundary](../../../GRHayLHD/src/Tabulated/outer_boundaries.c)
- [TabulatedEntropy boundary](../../../GRHayLHD/src/TabulatedEntropy/outer_boundaries.c)

## Related Pages

- [Conservative Recovery](conservative-recovery.md)
- [EOS and Entropy Variants](eos-entropy-variants.md)
- [Architecture Variables and Storage](../architecture/variables-and-storage.md)
