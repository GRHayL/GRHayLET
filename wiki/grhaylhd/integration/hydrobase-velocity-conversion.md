# HydroBase Velocity Conversion

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Integration](index.md)

## Scope and Non-Scope

This page records locally implemented conversion between HydroBase's Valencia
velocity and GRHayLHD's native `v^i = u^i/u^0`, including declared schedule
contexts and visible Lorentz-factor writeback. HydroBase convention semantics,
successful scheduling, and numerical validity remain external.

## Summary

Forward conversion computes native velocity as
`v_native^i = alpha U_HydroBase^i - beta^i`. Reverse conversion computes
`U_HydroBase^i = (v_native^i + beta^i)/alpha` and writes `w_lorentz` from a
local metric contraction. Initial and analysis reverse calls have a nonzero
cadence guard; leakage scheduling does not, while converter performs modulo by
that cadence.

## Variant Applicability

| Applicability | Conversion behavior |
| --- | --- |
| Common | Same native/Valencia conversion, metric inputs, cadence check, and `w_lorentz` expression in every EOS/entropy mode. |
| Hybrid/Simple | No mode-specific conversion branch. |
| Hybrid/Simple+Entropy | No mode-specific conversion branch. |
| Tabulated | No mode-specific conversion branch. |
| Tabulated+Entropy | No mode-specific conversion branch. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `INT-VEL-01` | Forward converter visibly computes `alpha U^i - beta^i` for all three components. | visible-implementation | Forward converter | `c:GRHayLHD/src/convert_HydroBase_to_GRHayLHD.c#symbol=convert_HydroBase_to_GRHayLHD` |
| `INT-VEL-02` | Reverse converter visibly computes `(v^i + beta^i)/alpha`, a metric-contraction Lorentz factor, and warning paths. | visible-implementation | Reverse converter | `c:GRHayLHD/src/convert_GRHayLHD_to_HydroBase.c#symbol=convert_GRHayLHD_to_HydroBase` |
| `INT-VEL-03` | Initial-data schedule declares HydroBase-to-native conversion before Prim2Con. | declared | Initial conversion schedule | `ccl:GRHayLHD/schedule.ccl#schedule=convert_HydroBase_to_GRHayLHD` |
| `INT-VEL-04` | Analysis schedule declares reverse conversion only when cadence is nonzero. | declared | Analysis schedule context | `ccl:GRHayLHD/schedule.ccl#schedule=convert_GRHayLHD_to_HydroBase?context=CCTK_ANALYSIS` |
| `INT-VEL-05` | Leakage schedule declares reverse conversion without the positive cadence guard used by other reverse contexts. | declared | Leakage schedule context | `ccl:GRHayLHD/schedule.ccl#schedule=convert_GRHayLHD_to_HydroBase?context=GRHayLHD_RHS` |
| `INT-VEL-06` | Leakage declaration omits metric reads and `w_lorentz` write visible in the converter body. | unresolved | Declaration/body mismatch | `c:GRHayLHD/src/convert_GRHayLHD_to_HydroBase.c#symbol=convert_GRHayLHD_to_HydroBase` |

## Details

### Ownership boundary and formulas

GRHayLHD inherits HydroBase and directly reads or writes HydroBase-owned
`vel`, `w_lorentz`, `rho`, `press`, `eps`, and mode-dependent thermodynamic
fields. It owns distinct `vx`, `vy`, and `vz` because local interface text
defines them as `u^i/u^0` rather than HydroBase's Valencia representation.

For lapse `alpha`, shift `beta^i`, HydroBase velocity `U^i`, and native
velocity `v^i`, visible assignments are:

```text
v^i = alpha U^i - beta^i
U^i = (v^i + beta^i) / alpha
```

Reverse conversion defines `q^i = v^i + beta^i` and computes

```text
w_lorentz = 1 / sqrt(1 - gamma_ij q^i q^j / alpha^2).
```

It emits an informational warning when contraction exceeds `1.0`, computes
the square root anyway, and emits another message when `W/alpha` is NaN. Local
code does not visibly clamp this path.

### Declared schedule contexts

- `HydroBase_Prim2ConInitial` declares forward conversion before Prim2Con.
- Same initial group optionally declares reverse conversion after Con2Prim
  when `Convert_to_HydroBase_every` is nonzero.
- `CCTK_ANALYSIS` declares reverse conversion under same nonzero guard, after
  `ML_BSSN_evolCalcGroup` and before named diagnostics.
- NRPyLeakageET-active RHS branch declares reverse conversion after flux RHS,
  outside that positive cadence guard.

Reverse function returns unless `cctk_iteration % Convert_to_HydroBase_every`
is zero. Its visible body writes both HydroBase velocity and `w_lorentz`;
leakage schedule's declared `WRITES` set names only velocity. The body also
reads all six spatial-metric fields, while leakage's declared `READS` names
only lapse, shift, and native velocity. No runtime effect is inferred from
these static declaration/body differences; see
[GRH-0011](../contradictions.md#grh-0011).

## Caveats

- Local formulas establish assignments, not equivalence to external
  HydroBase semantics or safe behavior for all metric/velocity values.
- Leakage guard versus modulo precondition is open
  [GRH-0006](../contradictions.md#grh-0006); no runtime crash is claimed.
- Leakage declaration/body field-set mismatch is open
  [GRH-0011](../contradictions.md#grh-0011).
- Schedule declarations establish intent only; they do not prove calls ran.

## Sources

- [Interface declarations](../../../GRHayLHD/interface.ccl)
- [Schedule declarations](../../../GRHayLHD/schedule.ccl)
- [HydroBase-to-native converter](../../../GRHayLHD/src/convert_HydroBase_to_GRHayLHD.c)
- [Native-to-HydroBase converter](../../../GRHayLHD/src/convert_GRHayLHD_to_HydroBase.c)

## Related Pages

- [Declared Schedule Lifecycle](../architecture/schedule-lifecycle.md)
- [Primitive-Conservative Conversion](../evolution/primitive-conservative-conversion.md)
- [GRHayLib Contract](grhaylib-contract.md)
