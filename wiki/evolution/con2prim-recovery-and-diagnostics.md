# Con2Prim Recovery and Diagnostics

> Status: confirmed · Last reconciled: 10-02-2026
> Up: [Evolution](index.md)

## Summary

Recovery uses configured primary/backup inversion, convex blends with a fixed
neighbor mean, then a supported terminal fallback. Hybrid EOS permits Font1D;
Simple and tabulated EOS reset to atmosphere after exhausted supported retries.
Terminal resets add 100 to the local diagnostic before its final publication.

Claim evidence:
- Claim: The four local recovery bodies preserve constant neighbor means and directly copy the full-neighbor candidate at weight four, supply current-state seeds when automatic guessing is disabled, and retain the terminal-reset marker; this does not establish external solver success.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/src/*/conservs_to_prims.c, recovery ladders and failure_checker writes`
- Corroboration: registered `IllinoisGRMHD/schedule.ccl, affected declarations`
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=inspected-not-run; options=local source and declaration inspection; date=10-02-2026`

## Detail

### Recovery ladder shared by all families

Every primitive/conservative carrier starts with deterministic zero placeholders.
When `calc_primitive_guess` is disabled in GRHayL's parameter object, callers load
current rho, pressure, epsilon, velocities, and active family fields, enforce
primitive limits, and reconstruct `u0`. They do not rely on the noncheckpointed
`u0` grid value after restart. B is always read from the centered grid fields;
tabulated temperature is also seeded from the grid. Each neighboring retry
restores the seed before entering the external multi-method solver.

1. Positive `rho_star` is undensitized and passed to the configured multi-method
   solver; Hybrid/Simple first apply conservative limits.
2. A NaN product over active outputs marks the returned state singular.
3. Nonpositive density selects constant atmosphere and adds the ones marker.
4. On error, the bounded 3-by-3-by-3 neighborhood excludes the central point.
   For fixed neighbor count N, sum S, center C, and attempt w=1..4, each
   conservative blend is `(w/4)*(S/N) + (1-w/4)*C`. For w=4 the neighbor
   mean `S/N` is assigned directly, so a nonfinite center does not enter the
   candidate. Finite results still require finite neighbor components; signed
   zero can differ from the old arithmetic blend. Empty neighborhoods skip
   averaging and proceed to the final policy. Active entropy/Ye participate.
5. Only `ghl_eos_hybrid` permits the local Font1D fallback. Simple and tabulated
   EOS have no local emergency Font attempt. Remaining failure resets atmosphere.
6. Final primitive enforcement and a separate conservative recomputation publish
   the repaired state. The latter loop preserves deterministic neighbor reads.

The compiled TabulatedEntropy body remains present, but the startup guard rejects
its runtime selection. Hybrid entropy supports only one cold segment with equal
cold/thermal Gamma; the proxy is `hybrid_entropy`, and physical HydroBase entropy
is marked unavailable. See [State and EOS Modes](state-and-eos-modes.md).

### Counters and repair encoding

Each call counts every point satisfying `sqrt_detgamma > psi6threshold` in the
selected population, independently of recovery success. Only terminal failures
increment its failure numerator. The label states this threshold criterion;
it does not claim independent apparent-horizon detection.

| Place | Meaning | Final contribution |
| --- | --- | --- |
| ones | nonpositive-density atmosphere reset | local +1 |
| tens | final primitive velocity limiting | local +10 |
| hundreds | exhausted supported recovery, atmosphere reset | local +100 |
| thousands | first backup flag | `1000 * diagnostics.backup[0]` |
| ten-thousands | tau reset | `10000 * diagnostics.tau_fix` |
| hundred-thousands | momentum reset | `100000 * diagnostics.Stilde_fix` |

`failure_checker` is overwritten at every RK substep. Verbose reductions also
report all three backup flags, averaging/Font attempts, inversion iterations,
and differences between original and recomputed conservatives. These local
counters do not establish behavior inside the external solvers.

## Sources

- `IllinoisGRMHD/src/Hybrid/conservs_to_prims.c` —
  `IllinoisGRMHD_hybrid_conservs_to_prims`, including Font1D and final write.
- `IllinoisGRMHD/src/HybridEntropy/conservs_to_prims.c` —
  `IllinoisGRMHD_hybrid_entropy_conservs_to_prims`, including entropy retry and
  diagnostics.
- `IllinoisGRMHD/src/Tabulated/conservs_to_prims.c` —
  `IllinoisGRMHD_tabulated_conservs_to_prims`, direct terminal atmosphere path.
- `IllinoisGRMHD/src/TabulatedEntropy/conservs_to_prims.c` —
  `IllinoisGRMHD_tabulated_entropy_conservs_to_prims`, combined extras and
  direct terminal atmosphere path.
- `IllinoisGRMHD/interface.ccl` — group `failure_checker` and RK-substep
  overwrite warning.
- `IllinoisGRMHD/param.ccl` — parameter `verbose`.
- `IllinoisGRMHD/schedule.ccl` — group `IllinoisGRMHD_conservs_to_prims`.

## See Also

- Parent: [Evolution](index.md)
- Depends on: [Primitive-Conservative Conversion](primitive-conservative-conversion.md)
- See also: [Matter Boundaries and Perturbations](matter-boundaries-and-perturbations.md)
