# HydroBase, GRHayLib, and Tmunu Boundary

> Status: confirmed · Last reconciled: 07-17-2026
> Up: [Integration](index.md)

## Summary

IllinoisGRMHD inherits ADMBase, HydroBase, TmunuBase, and GRHayLib. Local code
converts velocity and magnetic representations at HydroBase boundaries, calls
GRHayL through declared headers, and optionally adds locally assembled
stress-energy components to TmunuBase. This page describes only that visible
boundary; it does not infer GRHayL, HydroBase, TmunuBase, or NRPyLeakageET
internals.

Claim evidence:

- Claim: Local ingress/egress converts velocity and magnetic representations;
  optional local Tmunu routine adds returned components to inherited storage,
  without establishing external-library internals or observed execution.
- Role: public/scientific contract
- Deciding authority: `IllinoisGRMHD/src/convert_HydroBase_to_IllinoisGRMHD.c::convert_HydroBase_to_IllinoisGRMHD`,
  `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c::convert_IllinoisGRMHD_to_HydroBase`,
  and `IllinoisGRMHD/src/compute_Tmunu.c::IllinoisGRMHD_compute_Tmunu`
- Corroboration: `IllinoisGRMHD/schedule.ccl` conversion and `AddToTmunu` schedule blocks
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-run; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=static source inspection; date=07-17-2026`

## Detail

### Declared dependency boundary

- `configuration.ccl` requires HDF5 and GRHayL.
- `interface.ccl` inherits `ADMBase`, `Tmunubase`, `HydroBase`, and
  `GRHayLib`, and consumes `GRHayLib.h`.
- Calls named `ghl_*` below are opaque external calls at this scope. Local
  argument preparation and result writes are visible; library algorithms and
  guarantees are not.

### HydroBase ingress

`convert_HydroBase_to_IllinoisGRMHD` is scheduled in
`IllinoisGRMHD_Prim2Con2Prim`, inside `HydroBase_Prim2ConInitial`. It reads
HydroBase `vel`, `Avec`, and `Aphi`, plus ADM lapse and shift, then writes the
IllinoisGRMHD velocity and vector/scalar potentials.

IllinoisGRMHD stores `v^i = u^i/u^0`; its interface explicitly says this is
not HydroBase's Valencia velocity. For each component, ingress computes

```text
v_Illinois^i = lapse * vel_HydroBase^i - shift^i
```

With `rescale_magnetics=yes`, ingress multiplies each HydroBase `Avec`
component by `(4*pi)^(-1/2)`; with `no`, factor is one. It copies
`HydroBase::Aphi` directly to `phitilde`. No determinant is computed in this
conversion routine: `interface.ccl` supplies semantic declaration
`phitilde = sqrt(gamma) Phi`, while ingress assumes `Aphi` already represents
quantity to copy. Centered and staggered B are built later from A by
`IllinoisGRMHD_compute_B_and_Bstagger_from_A`.

HydroBase density, pressure, internal energy, entropy, electron fraction, and
temperature are used directly by scheduled evolution variants; this ingress
routine does not duplicate them.

### HydroBase egress and cadence

`convert_IllinoisGRMHD_to_HydroBase` computes

```text
vel_HydroBase^i = (v_Illinois^i + shift^i) / lapse
```

It computes `w_lorentz` from ADM spatial metric and converted velocity. It
also writes centered B to `HydroBase::Bvec`, multiplying by `(4*pi)^(1/2)`
when `rescale_magnetics=yes`, or by one otherwise.

Locally declared call sites are:

- initial conversion after `IllinoisGRMHD_conservs_to_prims`, present when
  local `Convert_to_HydroBase_every` is nonzero;
- `CCTK_ANALYSIS`, with declared ordering before named diagnostics, also
  present when local cadence is nonzero;
- after flux RHS evaluation when thorn `NRPyLeakageET` is active;
- two equivalent initial/analysis sites inside retained
  `ID_converter_ILGRMHD` compatibility gate.

If thorn `Convert_to_HydroBase` is active, the routine reads that thorn's
cadence, returns when it is zero, and otherwise runs only on divisible
iterations. When that thorn is inactive, it evaluates
`cctk_iteration % IllinoisGRMHD::Convert_to_HydroBase_every` without first
guarding zero.

Leakage schedule declaration names only `HydroBase::vel` as written, although
the routine also assigns `w_lorentz` and `Bvec`. This schedule site is not
gated by IllinoisGRMHD cadence. Therefore, with NRPyLeakageET active,
`Convert_to_HydroBase` inactive, and local cadence at its default zero, local
code evaluates integer remainder by zero: undefined C behavior that may trap.
A positive cadence or code-level zero guard is required. No NRPyLeakageET
internals or observed run are inferred.

Claim evidence:

- Claim: The leakage-gated call can evaluate integer remainder by zero when `Convert_to_HydroBase` is inactive and local cadence is zero; this is a local undefined-behavior path, not an observed run result.
- Role: descriptive behavior
- Deciding authority: `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c::convert_IllinoisGRMHD_to_HydroBase`, unguarded local-cadence remainder
- Corroboration: `IllinoisGRMHD/param.ccl::Convert_to_HydroBase_every` default zero and `IllinoisGRMHD/schedule.ccl::NRPyLeakageET` call site
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=inspected-not-run; options=NRPyLeakageET active, Convert_to_HydroBase inactive, local cadence zero; date=07-17-2026`

ThornGuide recommends matching conversion cadence to diagnostics and says
more frequent copying slows a simulation. This is attributed design advice,
not a performance measurement made by this KB.

### Tmunu handoff

When `update_Tmunu=yes`, schedule places `IllinoisGRMHD_compute_Tmunu` in
`AddToTmunu`. `IllinoisGRMHD_RegisterVars` also registers TmunuBase scalar,
vector, and tensor groups as constrained MoL groups under same condition.

For every local grid point, routine:

1. builds GRHayL metric and auxiliary objects from ADM lapse, shift, and
   spatial metric;
2. builds primitive object from HydroBase `rho`, `press`, and `eps`, plus
   IllinoisGRMHD `vx/vy/vz`, `u0`, and centered B;
3. calls `ghl_compute_TDNmunu`;
4. adds, using `+=`, ten returned components to `eTtt`, `eTtx`, `eTty`,
   `eTtz`, `eTxx`, `eTxy`, `eTxz`, `eTyy`, `eTyz`, and `eTzz`.

This establishes additive local writes, not stress-energy initialization,
external tensor conventions, or runtime execution.

## Sources

- [`IllinoisGRMHD/configuration.ccl`](../../IllinoisGRMHD/configuration.ccl) —
  `requires HDF5 GRHayL`.
- [`IllinoisGRMHD/interface.ccl`](../../IllinoisGRMHD/interface.ccl) —
  `inherits`, `grmhd_velocities`, `phitilde`, and external include declarations.
- [`IllinoisGRMHD/src/convert_HydroBase_to_IllinoisGRMHD.c`](../../IllinoisGRMHD/src/convert_HydroBase_to_IllinoisGRMHD.c) —
  `convert_HydroBase_to_IllinoisGRMHD`.
- [`IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c`](../../IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c) —
  `convert_IllinoisGRMHD_to_HydroBase`.
- [`IllinoisGRMHD/src/compute_Tmunu.c`](../../IllinoisGRMHD/src/compute_Tmunu.c) —
  `IllinoisGRMHD_compute_Tmunu`.
- [`IllinoisGRMHD/src/MoL_registration.c`](../../IllinoisGRMHD/src/MoL_registration.c) —
  `IllinoisGRMHD_RegisterVars` Tmunu condition.
- [`IllinoisGRMHD/schedule.ccl`](../../IllinoisGRMHD/schedule.ccl) —
  `IllinoisGRMHD_Prim2Con2Prim`, `IllinoisGRMHD_RHS`, `CCTK_ANALYSIS`,
  `AddToTmunu`, and compatibility schedule blocks.
- [`IllinoisGRMHD/param.ccl`](../../IllinoisGRMHD/param.ccl) — cadence,
  rescaling, and Tmunu controls.
- [`IllinoisGRMHD/doc/documentation.tex`](../../IllinoisGRMHD/doc/documentation.tex) —
  `Parameters` cadence guidance.

## See Also

- Parent: [Integration](index.md)
- Depends on: [Parameters and Runtime Controls](parameters-and-runtime-controls.md)
- Implements: [Schedule Lifecycle](../architecture/schedule-lifecycle.md)
- See also: [Staggered State and Magnetic Reconstruction](../magnetics/staggered-state-and-magnetic-reconstruction.md)
