# Primitive-Conservative Conversion

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Evolution](index.md)

## Scope and Non-Scope

This page records visible primitive-to-conservative dataflow in all four
variant implementations and its declared initial ordering. GRHayLib helper
semantics and numerical validity remain external.

## Summary

Every variant loops over the local grid, initializes ADM metric and auxiliary
objects, sets all three magnetic components to zero, loads mode primitives,
calls primitive limiting plus `u0` computation, computes conservatives, and
writes limited primitives plus core and optional conservative state.

## Variant Applicability

| Applicability | Extra primitive input/writeback | Extra conservative writeback | Declared sync |
| --- | --- | --- | --- |
| Common | `rho`, `press`, `eps`, native velocity, `u0` | `rho_star`, `tau`, three momentum components | Core conservatives |
| Hybrid/Simple | None | None | Core conservatives |
| Hybrid/Simple+Entropy | `entropy` | `ent_star` | Core plus `ent_star` |
| Tabulated | `Y_e`, `temperature` | `Ye_star` | Core plus `Ye_star` |
| Tabulated+Entropy | `entropy`, `Y_e`, `temperature` | `ent_star`, `Ye_star` | Core plus both optional groups |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `EV-P2C-01` | Hybrid conversion visibly zeros `BU`, limits primitives, computes conservatives, and writes core state. | visible-implementation | Variant function | `c:GRHayLHD/src/Hybrid/prims_to_conservs.c#symbol=GRHayLHD_hybrid_prims_to_conservs` |
| `EV-P2C-02` | HybridEntropy conversion additionally maps entropy. | visible-implementation | Variant function | `c:GRHayLHD/src/HybridEntropy/prims_to_conservs.c#symbol=GRHayLHD_hybrid_entropy_prims_to_conservs` |
| `EV-P2C-03` | Tabulated conversion additionally maps electron fraction and temperature. | visible-implementation | Variant function | `c:GRHayLHD/src/Tabulated/prims_to_conservs.c#symbol=GRHayLHD_tabulated_prims_to_conservs` |
| `EV-P2C-04` | TabulatedEntropy conversion maps both optional state families. | visible-implementation | Variant function | `c:GRHayLHD/src/TabulatedEntropy/prims_to_conservs.c#symbol=GRHayLHD_tabulated_entropy_prims_to_conservs` |
| `EV-P2C-05` | Initial schedule orders HydroBase velocity conversion before Prim2Con. | declared | Initial converter schedule | `ccl:GRHayLHD/schedule.ccl#schedule=convert_HydroBase_to_GRHayLHD` |
| `EV-P2C-06` | Hybrid Prim2Con schedule declares `HydroBase::eps` in `READS`. | declared | Hybrid schedule block | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_hybrid_prims_to_conservs` |
| `EV-P2C-07` | HybridEntropy Prim2Con schedule declares `HydroBase::eps` in `READS`. | declared | HybridEntropy schedule block | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_hybrid_entropy_prims_to_conservs` |
| `EV-P2C-08` | Tabulated Prim2Con schedule declares `HydroBase::eps` in `READS`. | declared | Tabulated schedule block | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_tabulated_prims_to_conservs` |
| `EV-P2C-09` | TabulatedEntropy Prim2Con schedule declares `HydroBase::eps` in `READS`. | declared | TabulatedEntropy schedule block | `ccl:GRHayLHD/schedule.ccl#schedule=GRHayLHD_tabulated_entropy_prims_to_conservs` |
| `EV-P2C-10` | Hybrid body does not load `eps[index]` into `prims` and later writes `prims.eps`. | visible-implementation | Hybrid body | `c:GRHayLHD/src/Hybrid/prims_to_conservs.c#symbol=GRHayLHD_hybrid_prims_to_conservs` |
| `EV-P2C-11` | HybridEntropy body does not load `eps[index]` into `prims` and later writes `prims.eps`. | visible-implementation | HybridEntropy body | `c:GRHayLHD/src/HybridEntropy/prims_to_conservs.c#symbol=GRHayLHD_hybrid_entropy_prims_to_conservs` |
| `EV-P2C-12` | Tabulated body does not load `eps[index]` into `prims` and later writes `prims.eps`. | visible-implementation | Tabulated body | `c:GRHayLHD/src/Tabulated/prims_to_conservs.c#symbol=GRHayLHD_tabulated_prims_to_conservs` |
| `EV-P2C-13` | TabulatedEntropy body does not load `eps[index]` into `prims` and later writes `prims.eps`. | visible-implementation | TabulatedEntropy body | `c:GRHayLHD/src/TabulatedEntropy/prims_to_conservs.c#symbol=GRHayLHD_tabulated_entropy_prims_to_conservs` |

## Details

### Initial ordering

Local schedule declares `convert_HydroBase_to_GRHayLHD` first in
`GRHayLHD_Prim2Con2Prim`. Optional primitive perturbation follows converter and
precedes Prim2Con; concrete Prim2Con group follows converter. Velocity formula
and external HydroBase convention are owned by Integration; ordering here is a
CCL declaration only.

### Common visible algorithm

Each routine constructs `ghl_metric_quantities` from lapse, shift, and six
spatial-metric components, then calls `ghl_compute_ADM_auxiliaries`. It loads a
`ghl_primitive_quantities` object and explicitly assigns
`prims.BU[0] = prims.BU[1] = prims.BU[2] = 0.0`. It calls
`ghl_enforce_primitive_limits_and_compute_u0`, passes any returned error to
`ghl_abort_if_error`, then calls `ghl_compute_conservs`.

Writeback stores possibly changed `rho`, `press`, `eps`, native velocities,
and `u0`. Core conservative writeback maps `cons.rho`, `cons.tau`, and
`cons.SD[0..2]` to five evolved fields. Variant extras follow canonical mode
matrix. CCL declares sync of core conservatives plus each active optional
conservative group.

### Independent variant differences

- Hybrid loads base primitives only.
- HybridEntropy also loads and writes `entropy`, then writes
  `cons.entropy` to `ent_star`.
- Tabulated loads and writes `Y_e` and `temperature`, then writes `cons.Y_e`
  to `Ye_star`.
- TabulatedEntropy performs both optional mappings.

No routine loads magnetic state from a gridfunction; zero `BU` is local visible
invariant for these conversion paths.

## Caveats

- Helper names and call order are visible; their math, limit policy, error
  contract, and conservation properties are delegated to GRHayLib.
- CCL `SYNC` proves declared synchronization intent, not completed exchange.
- All four Prim2Con schedule blocks declare `HydroBase::eps` in `READS`, but
  none of the four visible variant bodies loads `eps[index]` into `prims`;
  each later writes `prims.eps` after the GRHayLib limit helper. This
  declaration/body read-set mismatch is tracked under
  [GRH-0013](../contradictions.md#grh-0013); runtime consequences remain
  unverified.

## Sources

- [Schedule declarations](../../../GRHayLHD/schedule.ccl)
- [Hybrid conversion](../../../GRHayLHD/src/Hybrid/prims_to_conservs.c)
- [HybridEntropy conversion](../../../GRHayLHD/src/HybridEntropy/prims_to_conservs.c)
- [Tabulated conversion](../../../GRHayLHD/src/Tabulated/prims_to_conservs.c)
- [TabulatedEntropy conversion](../../../GRHayLHD/src/TabulatedEntropy/prims_to_conservs.c)

## Related Pages

- [EOS and Entropy Variants](eos-entropy-variants.md)
- [Declared Schedule Lifecycle](../architecture/schedule-lifecycle.md)
- [Architecture Variables and Storage](../architecture/variables-and-storage.md)
