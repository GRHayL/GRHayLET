# HydroBase, GRHayLib, and Tmunu Boundary

> Status: confirmed · Last reconciled: 10-02-2026
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

HydroBase density, pressure, internal energy, tabulated entropy, electron fraction, and
temperature are used directly by scheduled evolution variants; this ingress
routine does not duplicate them.

### HydroBase egress and cadence

`convert_IllinoisGRMHD_to_HydroBase` computes

```text
vel_HydroBase^i = (v_Illinois^i + shift^i) / lapse
```

It computes `w_lorentz` from ADM metric and converted velocity and copies
canonical normalized centered B to `HydroBase::Bvec`. The import compatibility
switch does not multiply exported B by `sqrt(4*pi)`.

The mandatory converter performs no cadence arithmetic. Leakage RHS and legacy
initialization call it directly. `IllinoisGRMHD_convert_HydroBase_diagnostics`
handles analysis cadence, including the retained old-thorn parameter lookup,
and guards zero before modulo. Modern initial conversion remains conditional
on requesting export. Modern guarded occurrences are excluded when the legacy
initializer owns the corresponding call, avoiding duplicate scheduling.

Every converter occurrence declares metric, lapse, shift, velocity, and B reads
and velocity, Lorentz-factor, and Bvec writes. Tmunu's additive destination arrays
are also declared as reads at its own occurrence. These declarations do not
prove communication/validity behavior in any driver.

Claim evidence:
- Claim: Mandatory and diagnostic conversion are separate, canonical B export is independent of legacy import, and schedule access includes actual fields.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c, both entry points`
- Corroboration: registered `IllinoisGRMHD/schedule.ccl, affected declarations`
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=inspected-not-run; options=local source and declaration inspection; date=10-02-2026`


IllinoisGRMHD installs no startup check that rejects smallbPoynET and cannot
identify a consumer's revision. A consumer reads the exported `Bvec` as
the canonical normalized field and must not divide it by `sqrt(4*pi)` again.
The ThornGuide documents that the updated smallbPoynET diagnostic follows this
convention, that earlier smallbPoynET revisions do not, and that a diagnostic
reads the exported fields only at its own cadence, so the governing
`Convert_to_HydroBase_every` should be a positive divisor of that cadence.

Locally, the default `Convert_to_HydroBase_every=0` disables the diagnostic
export, so IllinoisGRMHD does not refresh `Bvec` for a diagnostic consumer; the
mandatory leakage and compatibility conversions write only on their own
schedules. A positive value exports on divisible iterations. When the legacy
`Convert_to_HydroBase` thorn is active, the scheduled diagnostic wrapper applies
that thorn's `Convert_to_HydroBase::Convert_to_HydroBase_every` instead, with the
same zero guard and modulo test; unless `ID_converter_ILGRMHD` is active, the
wrapper is scheduled only if `IllinoisGRMHD::Convert_to_HydroBase_every` is also
nonzero. Consumer behavior and coupled schedule execution are not established by
this tree.

Claim evidence:
- Claim: The local exporter writes `Bvec` with a unit factor independently of import normalization and export cadence. The default `Convert_to_HydroBase_every=0` disables the diagnostic export, so IllinoisGRMHD does not refresh `Bvec` for a diagnostic consumer; a positive value exports on divisible iterations; and, when the legacy `Convert_to_HydroBase` thorn is active, the scheduled diagnostic wrapper uses that thorn's `Convert_to_HydroBase::Convert_to_HydroBase_every` instead, being scheduled without `ID_converter_ILGRMHD` only if the IllinoisGRMHD parameter is also nonzero. No local startup check rejects smallbPoynET or identifies a consumer's revision. Documented intent, not locally verified consumer behavior: the updated smallbPoynET consumes this convention, and a diagnostic reads exported fields only at its own cadence, so the governing `Convert_to_HydroBase_every` should be a positive divisor of that cadence.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c`, `convert_IllinoisGRMHD_to_HydroBase` (`mag_factor` and `Bvec` writes) and `IllinoisGRMHD_convert_HydroBase_diagnostics` (zero guard, modulo test, and legacy parameter lookup) for the local facts; registered `IllinoisGRMHD/doc/documentation.tex`, `Updating Old Parfiles` magnetic migration paragraph and the export-cadence paragraph that follows it, for the documented consumer convention and cadence guidance
- Corroboration: registered `IllinoisGRMHD/param.ccl`, `Convert_to_HydroBase_every` (range `0:*`, default `0`); registered `IllinoisGRMHD/schedule.ccl`, `CCTK_WRAGH` and `CCTK_ANALYSIS` occurrences, which declare no smallbPoynET check and gate the modern `CCTK_ANALYSIS` diagnostic export on nonzero `Convert_to_HydroBase_every`; the documented consumer convention and cadence guidance have no local corroboration because the consumer is external
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=not-applicable; options=local source and documentation inspection; date=10-04-2026`

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

### Hybrid entropy boundary

Hybrid/Simple recovery uses `IllinoisGRMHD::hybrid_entropy`, not physical
specific entropy. Its publication sites set HydroBase entropy to NaN to mark
that diagnostic unavailable; tabulated entropy evolution is rejected at startup.

Claim evidence:
- Claim: Hybrid/Simple proxy values are kept out of the physical HydroBase entropy field.
- Role: public/scientific contract
- Deciding authority: registered `IllinoisGRMHD/src/HybridEntropy/prims_to_conservs.c, primitive publication`
- Corroboration: registered `IllinoisGRMHD/schedule.ccl, affected declarations`
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-run; precision=not-applicable; GPU=not-applicable; restart=not-run; distributed=not-run; error_path=inspected-not-run; options=local source and declaration inspection; date=10-02-2026`

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
  `Parameters` cadence guidance and `Updating Old Parfiles` magnetic migration
  and export-cadence paragraphs.

## See Also

- Parent: [Integration](index.md)
- Depends on: [Parameters and Runtime Controls](parameters-and-runtime-controls.md)
- Implements: [Schedule Lifecycle](../architecture/schedule-lifecycle.md)
- See also: [Staggered State and Magnetic Reconstruction](../magnetics/staggered-state-and-magnetic-reconstruction.md)
