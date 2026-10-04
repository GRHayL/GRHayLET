# Beta Equilibrium

> Page status: reviewed · Last reviewed: 10-02-2026
> Up: [Initial Data](index.md)

## Scope and Non-Scope

This page owns the local declarations and visible dataflow described below.
Only GRHayLID sources are domain evidence. Framework execution, library
semantics, production-table provenance, and numerical validation remain
out of scope.

## Summary

The standalone routine now searches Ye at each actual density and requested
temperature using the visible residual mu_e-mu_n+mu_p. It recognizes
node/tolerance-qualified roots, rejects missing roots and failed lookups,
and checks the final residual. Atmosphere points reset density as well as
thermodynamics; README explicitly exempts them from equilibrium.

## Mode Applicability

| Applicability | Local surface |
| --- | --- |
| BetaEquilibrium | Standalone post-hydro operation with Ye/T storage and effective-bounds checks. |
| Common | An external hydro producer must initialize density in HydroBase_Initial. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `BETA-EQUILIBRIUM-01` | Beta scheduling declares Ye/T storage, rho reads, and five writes including rho. | declared | Named local source | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_BetaEquilibrium` |
| `BETA-EQUILIBRIUM-02` | The sentinel check precedes temperature-bound use, and the atmosphere branch resets rho. | visible-implementation | Named local source | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_BetaEquilibrium` |
| `BETA-EQUILIBRIUM-03` | The residual helper visibly computes mu_e-mu_n+mu_p from base-potential API outputs. | visible-implementation | Named local source | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_beta_residual` |
| `BETA-EQUILIBRIUM-04` | The local root helper scans effective Ye intervals, accepts tolerance-qualified nodes, and bisects sign changes. | visible-implementation | Named local source | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_beta_root` |
| `BETA-EQUILIBRIUM-05` | Chemical-potential interpolation and table conventions remain delegated. | out-of-scope | Named local source | `c:GRHayLID/src/BetaEquilibrium.c#call=ghl_tabulated_compute_P_eps_muhat_mue_mup_mun_from_T?function=GRHayLID_beta_residual` |
| `BETA-EQUILIBRIUM-06` | Residual tolerance defaults to 1e-8 MeV and is required positive and finite locally. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=beq_residual_tolerance` |

## Details

### Preconditions and ordering

ParamCheck requires Tabulated EOS and a non-none hydro producer. Beta requests
Ye/T storage independently of their selectors and runs after HydroBase_Initial.
The body checks storage and the -1 temperature sentinel before using temperature.
It requires positive finite beq_temperature within effective T bounds and
positive finite beq_residual_tolerance. Atmosphere metadata is checked against
effective bounds before iteration.

### Pointwise state construction

A nonfinite input density is diagnosed. For rho<=1.01*rho_atm, the body
sets density=rho_atm, Ye=Y_e_atm, T=T_atm and evaluates pressure/energy at
that complete tuple. README declares this atmosphere branch exempt from
beta equilibrium. Other points keep their density, use beq_temperature,
and search Ye within effective bounds.

The root helper starts at Y_e_min, scans table Ye nodes inside effective bounds
and the upper effective endpoint, and accepts residual magnitude at or below
tolerance. A sign change triggers at most 80 bisection iterations; no bracket,
nonfinite output, failed API status, or unmet residual tolerance is diagnosed.
The first root in increasing Ye is selected. The final residual is checked
again before publication. Neither the cached density-root builder nor its
minimum-Ye fallback is called. The residual arithmetic uses base potentials
rather than munu; this is visible dataflow, not proof of table provenance.

The density/Ye/T limits agree with post-beta entropy, which rejects invalid
states and does not clamp these fields. Beta does not write entropy or velocity.
GID-0007/0008 record resolved ordering and description mismatches.

## Caveats

Storage and schedule declarations do not prove allocation or execution.
Local guards and calls do not establish external error or interpolation
semantics. No checked-in GRHayLID test/parfile/oracle validates this path.

## Sources

- [schedule.ccl](../../../GRHayLID/schedule.ccl)
- [BetaEquilibrium.c](../../../GRHayLID/src/BetaEquilibrium.c)
- [param.ccl](../../../GRHayLID/param.ccl)

## Related Pages

- [Entropy Computation](entropy-computation.md)
- [Schedule Lifecycle](../architecture/schedule-lifecycle.md)
- [Parameters](../integration/parameters-and-configurations.md)
- [GID-0007](../contradictions.md#gid-0007)
- [GID-0008](../contradictions.md#gid-0008)
