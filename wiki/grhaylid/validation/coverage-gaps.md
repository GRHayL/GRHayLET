# Coverage Gaps

> Page status: reviewed · Last reviewed: 10-02-2026
> Up: [Validation](index.md)

## Scope and Non-Scope

This page owns the local declarations and visible dataflow described below.
Only GRHayLID sources are domain evidence. Framework execution, library
semantics, production-table provenance, and numerical validation remain
out of scope.

## Summary

The local source tree has no checked-in Cactus tests, parfiles, or numerical
oracles. The 1D_tests names identify initial-data generators. Added guards and
solver code do not establish framework execution or production physical validity.
The resolved source mismatches retain separate runtime evidence gaps.

## Mode Applicability

| Applicability | Local surface |
| --- | --- |
| Common | Local declarations/code support static claims only. |
| HydroTest1D | EOS closure, selector and direction coverage need configured execution. |
| HydroTest1D+Magnetic | Placement/curl/normalization need consumer validation. |
| BetaEquilibrium | Physical residual and table provenance require external evidence. |
| Entropy/Hybrid | Native-proxy consumer contract needs integration coverage. |
| Entropy/Tabulated | External input ownership and primitive replacement need scheduled execution. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `COVERAGE-GAPS-01` | Manifest lists implementation units and no test target. | coverage-gap | Named local source | `build:GRHayLID/src/make.code.defn#field=SRCS` |
| `COVERAGE-GAPS-02` | Local selector lists generators rather than a Cactus regression declaration. | coverage-gap | Named local source | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` |
| `COVERAGE-GAPS-03` | README explicitly states that no configured Cactus runtime or production-table validation is included. | declared | Named local source | `doc:GRHayLID/README#section=1. Purpose` |

## Details

Required future evidence includes a configured double-precision Cactus build,
full selector/storage matrix execution, initialized external Ye/T ownership,
metric-consistent Lorentz-factor production for actual consumers, magnetic
curl/boundary/normalization checks, and production-table chemical-potential
residuals at final rho/Ye/T. Atmosphere is deliberately exempt from equilibrium.

Useful local checks cover live sound-wave dispatch, Hybrid energy call selection,
staggering coordinate expressions, finite EOS returns, and parameter guard
coverage. Disposable probes are operational evidence, not new checked-in
GRHayLID numerical oracles. This KB does not import a sibling test suite or
external EOS implementation as branch-owned domain evidence.

Historical GID-0004/0005/0006/0010 are resolved by changed local evidence.
Their resolution does not claim a complete numerical regression pass.

### Ranked evidence proposals

These priorities order future evidence acquisition; they do not assert a
current pass or assign new defect severities. Each needs a configured isolated
Cactus environment and recorded inputs, commands, platform, results, and oracles.

| Priority | Proposed evidence | Decisive checks |
| --- | --- | --- |
| P0 | Double-precision Cactus build and complete selector/storage schedule matrix | Generated argument declarations, EOS API linkage, allocated destinations, rejected incompatible selections, initialized external input ownership, beta/entropy order. |
| P0 | Production-table beta and atmosphere cases | Independent base-potential residual at final off-node rho/T, node/no-root handling, effective limits, atmosphere EOS closure and the declared exception; trustworthy potential conventions. |
| P1 | Eight one-dimensional selectors for Simple/Hybrid with x/y/z coverage | Expected states, longitudinal velocity, unequal thermal/cold exponents and integration constants, admissible cold pressure. |
| P1 | Paired staggered/base magnetic cases and import checks | Component placement, metric/basis restriction, direct versus reconstructed fields, boundary curls, gauge behavior and normalization agreement. |
| P1 | Gas/sphere valid and invalid inputs | Both EOS regions, consistent outputs or diagnostic rejection, radius boundary and finite subluminal interior velocity. |
| P1 | Standalone and combined entropy cases | Native-proxy consumer agreement, vacuum/nonfinite rejection, physical Tabulated entropy, deliberate primitive replacement and preserved post-beta state. |
| P2 | Consumer-specific Lorentz-factor and legacy integration cases | Ordered metric-consistent W producer for an actual reader, replacement initialization and supported legacy settings. |

## Caveats

Storage and schedule declarations do not prove allocation or execution.
Local guards and calls do not establish external error or interpolation
semantics. No checked-in GRHayLID test/parfile/oracle validates this path.

## Sources

- [make.code.defn](../../../GRHayLID/src/make.code.defn)
- [param.ccl](../../../GRHayLID/param.ccl)
- [README](../../../GRHayLID/README)

## Related Pages

- [Hydro Tests](../initial-data/one-d-tests-hydro.md)
- [Magnetic Tests](../initial-data/one-d-tests-magnetic.md)
- [Beta](../initial-data/beta-equilibrium.md)
- [GID-0004](../contradictions.md#gid-0004)
- [GID-0005](../contradictions.md#gid-0005)
- [GID-0006](../contradictions.md#gid-0006)
- [GID-0010](../contradictions.md#gid-0010)
