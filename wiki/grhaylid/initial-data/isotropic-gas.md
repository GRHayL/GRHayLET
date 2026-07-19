# IsotropicGas Initial Data

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Initial Data](index.md)

## Scope and Non-Scope

This page records IsotropicGas precondition text, local parameter sentinels,
the visible tabulated-EOS call, and pointwise HydroBase assignments. It does
not establish error-macro termination, EOS-table behavior, returned values,
successful scheduling, or physical and numerical validity.

## Summary

`GRHayLID_IsotropicGas` visibly checks for Tabulated `EOS_type` and
`GRHayLID` selections for both `initial_Y_e` and `initial_temperature`. It
then applies `CHECK_PARAMETER` to three local controls, makes one
`ghl_tabulated_compute_P_eps_from_T` call, and writes uniform density,
electron fraction, temperature, pressure, internal energy, and zero velocity
throughout its point loop.

## Mode Applicability

| Applicability | Visible initial-data behavior |
| --- | --- |
| IsotropicGas | One uniform thermodynamic state, one delegated tabulated-EOS call, and zero velocity. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `ID-ISO-01` | Function visibly performs three precondition comparisons, three sentinel checks, one tabulated-EOS call, and uniform pointwise writes. | visible-implementation | Function body | `c:GRHayLID/src/IsotropicGas.c#symbol=GRHayLID_IsotropicGas` |
| `ID-ISO-02` | `IsotropicGas_rho` is a nonnegative real with forbidden `-1` sentinel and default `-1`. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_rho` |
| `ID-ISO-03` | `IsotropicGas_Y_e` is a nonnegative real with forbidden `-1` sentinel and default `-1`. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_Y_e` |
| `ID-ISO-04` | `IsotropicGas_temperature` is a nonnegative real with forbidden `-1` sentinel and default `-1`. | declared | Parameter declaration | `ccl:GRHayLID/param.ccl#parameter=IsotropicGas_temperature` |
| `ID-ISO-05` | ThornGuide prints the EOS selection as lowercase `"tabulated"` and the electron-fraction and temperature selections as `"GRHayLID"`. | declared | IsotropicGas subsection | `doc:GRHayLID/doc/documentation.tex#section=IsotropicGas` |

## Details

### Preconditions and sentinels

Visible function order compares:

1. `EOS_type` against `Tabulated`, calling `CCTK_ERROR` when unequal.
2. `initial_Y_e` against `GRHayLID`, calling `CCTK_ERROR` when unequal.
3. `initial_temperature` against `GRHayLID`, calling `CCTK_ERROR` when
   unequal.

It then expands `CHECK_PARAMETER` for `IsotropicGas_rho`,
`IsotropicGas_Y_e`, and `IsotropicGas_temperature`. Each parameter
declaration admits nonnegative values, marks `-1` forbidden, and defaults to
`-1`. The function visibly invokes the macro; resulting error behavior remains
external.

ThornGuide's IsotropicGas table prints `EOS_type` as `"tabulated"`; the body
compares it with `"Tabulated"`. The HydroBase selections are spelled
identically. Equivalence of the EOS spellings is external; see
[GID-0014](../contradictions.md#gid-0014).

### EOS call and uniform writes

Before the loop, code calls `ghl_tabulated_compute_P_eps_from_T` once with
`ghl_eos`, the three local parameter values, and addresses for local pressure
and internal-energy outputs. No local return code is captured. The loop then
assigns every point:

- `Y_e`, `rho`, and `temperature` from the three parameters;
- `press` and `eps` from the two output variables;
- all three `vel` components to zero.

Call semantics, table interpolation, bounds handling, units beyond parameter
descriptions, and output validity belong to external GRHayLib behavior.

## Caveats

- Three comparisons and calls are visible; whether they terminate execution is
  not established locally.
- One external EOS call is visible, but its returned pressure and internal
  energy are not validated by checked-in local evidence.
- Schedule description calls this three-dimensional setup a "1D test"; see
  [GID-0011](../contradictions.md#gid-0011).
- ThornGuide and the function use different EOS literal capitalization; see
  [GID-0014](../contradictions.md#gid-0014).
- Uniform writes do not prove a physically isotropic or valid state.

## Sources

- [IsotropicGas implementation](../../../GRHayLID/src/IsotropicGas.c)
- [Parameter declarations](../../../GRHayLID/param.ccl)
- [ThornGuide source](../../../GRHayLID/doc/documentation.tex)
- [Schedule declarations](../../../GRHayLID/schedule.ccl)

## Related Pages

- [ConstantDensitySphere Initial Data](constant-density-sphere.md)
- [GRHayLib Contract](../integration/grhaylib-contract.md)
- [GID-0011](../contradictions.md#gid-0011)
- [GID-0014](../contradictions.md#gid-0014)
