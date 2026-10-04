# Conservative Recovery

> Page status: reviewed · Last reviewed: 10-04-2026
> Up: [Evolution](index.md)

## Scope and Non-Scope

This page traces visible conservative-to-primitive control flow and final
writes in all four variants. GRHayLib solver math, success criteria, atmosphere
values, and physical validity remain external.

## Summary

All variants route non-positive conserved density to atmosphere, otherwise
undensitize and call `ghl_con2prim_multi_method`, screen recovered fields for
NaN, retry failures with bounded-neighborhood weighted conservative inputs,
and use atmosphere after terminal failure. Hybrid families alone visibly call
an explicit Font1D fallback restricted to Hybrid EOS. After recovery, all variants limit primitives,
write them, then use a separate loop to recompute and overwrite conservatives
while accumulating change diagnostics.

## Variant Applicability

| Applicability | Recovered extras | Visible pre-solver conservative limits | Explicit post-retry fallback | Recomputed extras |
| --- | --- | --- | --- | --- |
| Common | Base thermodynamics and velocity | Mode-dependent | Atmosphere after terminal error | Five core conservatives |
| Hybrid/Simple | None | Yes | Hybrid-only `ghl_hybrid_Font1D`, then atmosphere; Simple goes directly to atmosphere | None |
| Hybrid/Simple+Entropy | Entropy | Yes | Hybrid-only `ghl_hybrid_Font1D`, then atmosphere; Simple goes directly to atmosphere | `ent_star` |
| Tabulated | `Y_e`, temperature | No local `ghl_apply_conservative_limits` call | Atmosphere; no explicit Font1D call | `Ye_star` |
| Tabulated+Entropy | Entropy, `Y_e`, temperature | No local `ghl_apply_conservative_limits` call | Atmosphere; no explicit Font1D call | `ent_star`, `Ye_star` |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `EV-C2P-01` | Hybrid recovery visibly applies conservative limits, retries weighted neighborhoods, calls Font1D only for Hybrid EOS, writes primitives, and recomputes core conservatives. | visible-implementation | Full variant function | `c:GRHayLHD/src/Hybrid/conservs_to_prims.c#symbol=GRHayLHD_hybrid_conservs_to_prims` |
| `EV-C2P-02` | HybridEntropy follows Hybrid fallback structure and includes entropy in recovery and recomputation. | visible-implementation | Full variant function | `c:GRHayLHD/src/HybridEntropy/conservs_to_prims.c#symbol=GRHayLHD_hybrid_entropy_conservs_to_prims` |
| `EV-C2P-03` | Tabulated recovery includes electron fraction and temperature, but no explicit Font1D call. | visible-implementation | Full variant function | `c:GRHayLHD/src/Tabulated/conservs_to_prims.c#symbol=GRHayLHD_tabulated_conservs_to_prims` |
| `EV-C2P-04` | TabulatedEntropy includes both optional state sets and no explicit Font1D call. | visible-implementation | Full variant function | `c:GRHayLHD/src/TabulatedEntropy/conservs_to_prims.c#symbol=GRHayLHD_tabulated_entropy_conservs_to_prims` |
| `EV-C2P-05` | Interface declares `failure_checker` as per-substep-overwritten diagnostic. | declared | Diagnostic group | `ccl:GRHayLHD/interface.ccl#group=failure_checker` |

## Details

### Primary path and atmosphere path

Each function rejects `ghl_params->calc_prim_guess == false` before recovery.
It initializes complete primitive/conservative objects to zero, then loads
active mode conservatives, diagnostics, metric, and auxiliaries. Active
guesses are delegated to GRHayL; the guard does not implement previous-state guesses. Condition `cons.rho > 0.0`
enters solver path; its complement sets constant atmosphere, adds one to local
failure code, increments density-fix count, and marks success. The visible atmosphere gate and decoder both use non-positive density.

Positive-density Hybrid paths call `ghl_apply_conservative_limits` before
undensitization. Tabulated paths visibly skip that helper. All call
`ghl_undensitize_conservatives` and `ghl_con2prim_multi_method`. A product of
required recovered primitive fields is tested with `isnan`; optional entropy,
electron fraction, and temperature join product in relevant modes. NaN forces
local singular-error code.

### Bounded-neighborhood retries

On error, code bounds each coordinate to local grid and scans clipped
`i-1..i+1`, `j-1..j+1`, `k-1..k+1` neighborhood, excluding center. It sums
every active conservative field. A loop with `avg_weight` values 1 through 4
constructs `w*(neighbor_sum/n_avg) + (1-w)*center` for each active field,
with fixed neighbor count and `w = avg_weight/4`. With no neighbors it skips
these retries. It undensitizes, retries multi-method recovery, and repeats the NaN screen.
This is visible algorithmic structure; no numerical-quality or race-free claim
is inferred from comments.

For Hybrid EOS only, Hybrid and HybridEntropy reload/apply conservative
limits, undensitize, and call `ghl_hybrid_Font1D`. Simple skips this cold
fallback. Remaining errors or a failed NaN screen lead to atmosphere. Tabulated and TabulatedEntropy proceed directly from exhausted
weighted retries to atmosphere. All terminal paths update aggregate failures and the failure count above
`psi6threshold`. The regional denominator counts every grid point above that
threshold, including successful recoveries, and is labelled `AbovePsi6Threshold`.

### Post-recovery writes and recomputation

All paths call `ghl_enforce_primitive_limits_and_compute_u0` and abort through
helper on returned error. They write `rho`, `press`, `eps`, `u0`, and native
velocity. Entropy and tabulated families additionally write their mode
primitives.

A second OpenMP loop reconstructs metric and primitives, again zeros `BU`,
calls `ghl_compute_conservs`, and overwrites five core conservative fields plus
active `ent_star` and/or `Ye_star`. Before overwrite it snapshots original
conservatives; after computation it accumulates absolute differences and
denominators. Under `verbose == yes`, each variant prints aggregate backup,
limit, failure, iteration, and conservative-difference summaries; exact fields
follow mode.

### `failure_checker` legend versus write order

Source comments assign 1 to a nonpositive-density atmosphere reset, 10 to speed
limiting, 100 to exhaustion of all allowed recovery attempts, 1000 to backup
use, 10000 to a tau fix, and 100000 to a momentum fix. Every terminal branch
adds 100 to `local_failure_checker`, and the single final point assignment
preserves that contribution. This static reading does not establish a current
Cactus runtime diagnostic result.

### Dormant symmetry group name

All four variants visibly request `GRHayLHD::grmhd_conservatives` in their
equatorial blocks, matching the interface declaration. Only `Symmetry=none`
is locally selectable; the group name does not establish supported equatorial
symmetry behavior.

## Caveats

- Helper calls establish visible order, not solver semantics or recovery
  correctness.
- Absence of explicit Font1D in tabulated files does not exclude fallback
  inside external multi-method implementation.
- Comments calling second loop deterministic do not prove thread safety.

## Sources

- [Interface diagnostic declaration](../../../GRHayLHD/interface.ccl)
- [Hybrid recovery](../../../GRHayLHD/src/Hybrid/conservs_to_prims.c)
- [HybridEntropy recovery](../../../GRHayLHD/src/HybridEntropy/conservs_to_prims.c)
- [Tabulated recovery](../../../GRHayLHD/src/Tabulated/conservs_to_prims.c)
- [TabulatedEntropy recovery](../../../GRHayLHD/src/TabulatedEntropy/conservs_to_prims.c)

## Related Pages

- [EOS and Entropy Variants](eos-entropy-variants.md)
- [Primitive-Conservative Conversion](primitive-conservative-conversion.md)
- [Architecture Variables and Storage](../architecture/variables-and-storage.md)
