# Primitive-Conservative Conversion

> Status: confirmed · Last reconciled: 07-17-2026
> Up: [Evolution](index.md)

## Summary

Each family has a full-grid Prim2Con kernel: initialize metric data, assemble
family primitives, enforce primitive limits while computing `u0`, compute
conservatives, then write any corrected primitives and family conservatives.
After Con2Prim, every family performs a separate full-grid conservative
recomputation; failure recovery and boundary-fill repair have different
owners.

Claim evidence:

- Claim: All four local families use the stated Prim2Con stages and a separate post-Con2Prim conservative recomputation loop; external GRHayL internals are not claimed.
- Role: descriptive behavior
- Deciding authority: registered aggregates `illinoisgrmhd-hybrid`, `illinoisgrmhd-hybrid-entropy`, `illinoisgrmhd-tabulated`, and `illinoisgrmhd-tabulated-entropy`; four `*_prims_to_conservs` functions and post-recovery loops
- Corroboration: registered `IllinoisGRMHD/schedule.ccl`, `IllinoisGRMHD_prims_to_conservs` and `IllinoisGRMHD_conservs_to_prims` groups
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-run; options=Hybrid,HybridEntropy,Tabulated,TabulatedEntropy; date=07-17-2026`

## Detail

### Prim2Con kernel

The four stable entry points are
`IllinoisGRMHD_hybrid_prims_to_conservs`,
`IllinoisGRMHD_hybrid_entropy_prims_to_conservs`,
`IllinoisGRMHD_tabulated_prims_to_conservs`, and
`IllinoisGRMHD_tabulated_entropy_prims_to_conservs`. Each loops over all local
grid points and:

1. initializes `ghl_metric_quantities` from lapse, shift, and spatial metric,
   then computes local ADM auxiliaries;
2. assembles a `ghl_primitive_quantities` record from HydroBase primitives,
   IllinoisGRMHD velocity, and centered B, adding entropy and/or `Y_e` plus
   temperature according to the selected family;
3. calls `ghl_enforce_primitive_limits_and_compute_u0`, aborting through the
   local error boundary if it reports an error;
4. calls `ghl_compute_conservs`;
5. writes possibly modified rho, pressure, internal energy, velocity, and
   `u0`, plus applicable entropy, `Y_e`, and temperature; and
6. writes `rho_star`, `tau`, `Stildex/y/z`, plus `ent_star` and/or `Ye_star`.

The initial declared schedule places this group after A-to-B reconstruction,
then immediately schedules Con2Prim. This page does not own HydroBase ingress
or magnetic rescaling.

### Ordinary successful Con2Prim handoff

On the ordinary positive-density path, family Con2Prim code undensitizes the
conservatives and invokes its configured multi-method GRHayL call. After a
successful result it applies the same primitive-limit/`u0` call used by
Prim2Con and writes recovered primitives. Solver ordering and fallback policy
belong to [Con2Prim Recovery and Diagnostics](con2prim-recovery-and-diagnostics.md);
external solver behavior is not inferred.

### Common post-recovery recomputation

All four Con2Prim entry points then start a second OpenMP loop. They rebuild
metric and primitive records from the just-written grid functions, call
`ghl_compute_conservs`, and overwrite the hydrodynamic conservatives plus
family extras. They retain original values only long enough to accumulate
absolute-change numerators and normalization sums for verbose diagnostics.

The source comment states why this is a separate loop: neighbor averaging in
the recovery loop must remain deterministic and free of a race with
conservative rewrites. This is a source-stated implementation purpose, not a
general thread-safety guarantee.

### Ownership boundaries

- Recovery retries, atmosphere reset, repair encoding, and counters belong to
  [Con2Prim Recovery and Diagnostics](con2prim-recovery-and-diagnostics.md).
- Primitive limiting and conservative recomputation inside the four
  `hydro_outer_boundaries.c` functions belong to
  [Matter Boundaries and Perturbations](matter-boundaries-and-perturbations.md).
- `ghl_*` calls above describe only inputs, outputs, and error handling visible
  at IllinoisGRMHD call sites.

## Sources

- `IllinoisGRMHD/src/Hybrid/prims_to_conservs.c` —
  `IllinoisGRMHD_hybrid_prims_to_conservs`.
- `IllinoisGRMHD/src/HybridEntropy/prims_to_conservs.c` —
  `IllinoisGRMHD_hybrid_entropy_prims_to_conservs`.
- `IllinoisGRMHD/src/Tabulated/prims_to_conservs.c` —
  `IllinoisGRMHD_tabulated_prims_to_conservs`.
- `IllinoisGRMHD/src/TabulatedEntropy/prims_to_conservs.c` —
  `IllinoisGRMHD_tabulated_entropy_prims_to_conservs`.
- `IllinoisGRMHD/src/Hybrid/conservs_to_prims.c`,
  `IllinoisGRMHD/src/HybridEntropy/conservs_to_prims.c`,
  `IllinoisGRMHD/src/Tabulated/conservs_to_prims.c`, and
  `IllinoisGRMHD/src/TabulatedEntropy/conservs_to_prims.c` — post-recovery
  primitive writes and separate `ghl_compute_conservs` loops.
- `IllinoisGRMHD/schedule.ccl` — groups `IllinoisGRMHD_Prim2Con2Prim`,
  `IllinoisGRMHD_prims_to_conservs`, and `IllinoisGRMHD_conservs_to_prims`.

## See Also

- Parent: [Evolution](index.md)
- Depends on: [State and EOS Modes](state-and-eos-modes.md)
- See also: [Con2Prim Recovery and Diagnostics](con2prim-recovery-and-diagnostics.md)
- See also: [HydroBase, GRHayLib, and Tmunu](../integration/hydrobase-grhaylib-and-tmunu.md)
