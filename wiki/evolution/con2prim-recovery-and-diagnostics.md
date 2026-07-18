# Con2Prim Recovery and Diagnostics

> Status: contested · Last reconciled: 07-17-2026
> Up: [Evolution](index.md)

## Summary

Con2Prim uses primary inversion, weighted neighbor retries, then a
family-specific terminal fallback. Hybrid families try Font1D before resetting
to atmosphere; tabulated families reset directly. Active issue
[`CONTR-0002`](../contradictions.md#contr-0002):
the terminal `+100` update is overwritten by the final assignment, so current
code does not retain the documented hundreds marker.

Claim status: contested; contradiction: CONTR-0002.
Backlink: [`CONTR-0002` register entry](../contradictions.md#contr-0002).

Claim evidence:

- Claim: All four local variants overwrite the earlier terminal `failure_checker += 100` with a final assignment that omits 100; no runtime value was observed.
- Role: descriptive behavior
- Deciding authority: registered exact sources `illinoisgrmhd-hybrid-con2prim`, `illinoisgrmhd-hybrid-entropy-con2prim`, `illinoisgrmhd-tabulated-con2prim`, and `illinoisgrmhd-tabulated-entropy-con2prim`; terminal branches and final `failure_checker[index]` assignments
- Corroboration: same functions' decoder comments state the intended hundreds marker; Tabulated functions also lack a local Font1D call
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=inspected-not-run; options=Hybrid,HybridEntropy,Tabulated,TabulatedEntropy; date=07-17-2026`

## Detail

### Recovery ladder shared by all families

Each of the four `*_conservs_to_prims` functions initializes per-point GRHayL
diagnostics and metric data, loads centered B and family conservatives, and
uses the library's default guess; source comments explicitly say a prior-
timelevel guess is not implemented.

1. For positive `rho_star`, conservatives are divided by
   `sqrt_detgamma` and passed to `ghl_con2prim_multi_method`. Hybrid families
   first call `ghl_apply_conservative_limits`; tabulated families do not make
   that local call on the primary path.
2. A NaN product across every expected output field converts an otherwise
   returned result into `ghl_error_c2p_singular`. Entropy and tabulated fields
   are included only in families that carry them.
3. For nonpositive `rho_star`, the point is set to constant atmosphere,
   `local_failure_checker` gains 1, the rho-reset counter increments, and the
   retry ladder is skipped.
4. On inversion error, code sums the available neighbors in the bounded
   3-by-3-by-3 neighborhood, excluding the point itself. For original neighbor
   count `N`, neighbor sum `S`, center value `C`, and attempt `w=1..4`, each
   retried conservative is `((w/4)S + (1-w/4)C)/D_w`, with
   `D_w=(N+1,N+2,N+3,N+3)`. Thus current fourth retry is `S/(N+3)`, not a
   normalized full-neighbor average. Entropy and/or `Y_e` participate where
   evolved.
5. If retries still fail, Hybrid and HybridEntropy reload/limit the original
   conservatives, undensitize, and call `ghl_hybrid_Font1D`. Tabulated and
   TabulatedEntropy have no local Font1D call and proceed directly to a
   constant-atmosphere reset.
6. Terminal fallback increments failure counters and, when
   `sqrt_detgamma > ghl_params->psi6threshold`, the two locally named horizon
   counters. It then follows the common final primitive limiting and write
   path.

After the ladder, every variant calls
`ghl_enforce_primitive_limits_and_compute_u0`; a speed-limited result adds 10
to the local repair value. Entropy families write entropy, tabulated families
write `Y_e` and temperature, and the combined family writes all three.
All variants then run the separate conservative recomputation described in
[Primitive-Conservative Conversion](primitive-conservative-conversion.md).

### Counters and verbose output

When `verbose` equals `yes`, reductions report grid size, three GRHayL backup
counters, speed/rho fixes, averaged points, terminal failures, locally named
horizon counts, average inversion iterations, and relative/summed
conservative changes. Hybrid output additionally reports Font1D attempt count;
entropy/tabulated output adds entropy and/or `Y_e` change diagnostics. These
are current-call reductions. `interface.ccl` warns that `failure_checker` is
overwritten at every RK substep.

No local code establishes what each external multi-method attempt does, which
methods a configuration supports, or the semantics of GRHayL internals beyond
returned diagnostics consumed here.

### Actual `failure_checker` encoding

The decoder comment in every variant states:

| Place | Commented meaning | Executed final assignment |
| --- | --- | --- |
| ones | atmosphere reset when `rho_star < 0` | code tests `cons.rho <= 0`, then `local_failure_checker += 1` |
| tens | velocity limiting | `local_failure_checker += 10` |
| hundreds | C2P and Font failure | **not retained** |
| thousands | backups used | `1000 * diagnostics.backup[0]` |
| ten-thousands | tau reset | `10000 * diagnostics.tau_fix` |
| hundred-thousands | momentum reset | `100000 * diagnostics.Stilde_fix` |

On terminal failure, all four functions execute
`failure_checker[index] += 100`. Later in the same point iteration they
unconditionally execute:

```text
failure_checker[index] = local_failure_checker
                       + 1000*diagnostics.backup[0]
                       + 10000*diagnostics.tau_fix
                       + 100000*diagnostics.Stilde_fix;
```

Terminal failure never adds 100 to `local_failure_checker`; the earlier update
is discarded. Tabulated variants also use direct atmosphere fallback, so the
decoder's “Font Fix” wording does not describe their local path. Until
`CONTR-0002` is resolved by source reconciliation plus a targeted forced-
fallback run, this page remains contested and must not promise a hundreds
digit.

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
- See also: [`CONTR-0002`](../contradictions.md#contr-0002)
- See also: [Matter Boundaries and Perturbations](matter-boundaries-and-perturbations.md)
