# Test Inventory and Oracles

> Page status: reviewed · Last reviewed: 07-17-2026
> Up: [Validation](index.md)

## Scope and Non-Scope

This page inventories local test declarations, five parfile roles, and eleven
checked-in ASCII observations. It records declared tolerances and file content
only. No case was run or regenerated; current test success, numerical validity,
convergence, artifact lineage, and external test-harness semantics remain
unproven.

## Summary

`test/test.ccl` declares Balsara0 and TOV, each with `RELTOL 1e-12`.
Checked-in evidence includes one shipped example, two authored test inputs,
two checked-in companion configurations whose headers assert Cactus
generation, five Balsara0 x-direction ASCII files, and six TOV
density/pressure ASCII files across x/y/z. Authored inputs select Simple and
Hybrid EOS respectively but do not explicitly select shared entropy evolution.

## Variant Applicability

| Applicability | Local test evidence |
| --- | --- |
| Common | Two case declarations, five parfiles, and eleven ASCII files exist independently of observed execution. |
| Hybrid/Simple | Balsara0 authored input sets Simple; TOV authored input sets Hybrid; entropy state remains locally unresolved for both. |
| Hybrid/Simple+Entropy | No authored test input explicitly selects entropy. |
| Tabulated | No authored test input explicitly selects Tabulated EOS. |
| Tabulated+Entropy | No authored test input explicitly selects this combination. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `VAL-INV-01` | Test CCL declares Balsara0 with relative tolerance `1e-12`. | declared | Balsara0 test block | `test:GRHayLHD/test/test.ccl#case=Balsara0` |
| `VAL-INV-02` | Test CCL declares TOV with relative tolerance `1e-12`. | declared | TOV test block | `test:GRHayLHD/test/test.ccl#case=TOV` |
| `VAL-INV-03` | Authored Balsara0 input sets Simple EOS, disables Tmunu update, and requests five x-direction fields. | checked-in-observation | Authored Balsara0 input | `par:GRHayLHD/test/Balsara0.par#parameter=GRHayLib::EOS_type` |
| `VAL-INV-04` | Authored TOV input sets Hybrid EOS, names external Font1D as first backup routine, and requests HydroBase density and pressure. | checked-in-observation | Authored TOV input | `par:GRHayLHD/test/TOV.par#parameter=GRHayLib::con2prim_backup_routines[0]` |
| `VAL-INV-05` | Balsara0 authored input contains no explicit `GRHayLib::evolve_entropy` assignment. | checked-in-observation | Full authored-input inspection | `par:GRHayLHD/test/Balsara0.par#file` |
| `VAL-INV-06` | Five Balsara0 files contain x-direction numeric observations at iterations 0 and 10. | checked-in-observation | Representative Balsara0 file | `oracle:GRHayLHD/test/Balsara0/rho.x.asc#file` |
| `VAL-INV-07` | Six TOV files contain density/pressure numeric observations in x/y/z at iterations 0 and 2. | checked-in-observation | Representative TOV file | `oracle:GRHayLHD/test/TOV/hydrobase-rho.x.asc#file` |
| `VAL-INV-08` | Balsara0 companion contains checked-in Cactus-generation and original-path header assertions. | checked-in-observation | Whole companion file | `par:GRHayLHD/test/Balsara0/Balsara0.par#file` |
| `VAL-INV-09` | TOV authored input contains no explicit `GRHayLib::evolve_entropy` assignment. | checked-in-observation | Full authored-input inspection | `par:GRHayLHD/test/TOV.par#file` |
| `VAL-INV-10` | Shared entropy selection remains locally unresolved because authored-input omissions do not establish the external shared parameter value. | unresolved | Shared selection relationship | `ccl:GRHayLHD/param.ccl#parameter=evolve_entropy` |
| `VAL-INV-11` | TOV companion contains checked-in Cactus-generation and original-path header assertions. | checked-in-observation | Whole companion file | `par:GRHayLHD/test/TOV/TOV.par#file` |

## Details

### Declared cases

| Case | Test declaration | Authored input observations | Requested output | Checked-in oracle family |
| --- | --- | --- | --- | --- |
| Balsara0 | `TEST Balsara0`; `RELTOL 1e-12` | Iteration termination at 10; Simple EOS; `update_Tmunu=no`; copy matter boundary; no explicit entropy selection | `rho`, `press`, `vx`, `vy`, `vz`; x direction; every 10 iterations | Five x-direction files with iteration 0 and 10 blocks |
| TOV | `TEST TOV`; `RELTOL 1e-12` | Iteration limit 2; Hybrid EOS; external `con2prim_backup_routines[0]` set to `Font1D`; no explicit entropy selection; local Tmunu/boundary controls omitted | HydroBase `rho` and `press`; every 2 iterations | Six density/pressure files across x/y/z with iteration 0 and 2 rows |

Case declarations and tolerances are CCL intent. Adjacent authored parfiles and
oracle directories share case names, but local static evidence does not prove
a current harness invocation or result.

### Five parfile roles

| File | Role | Distinguishing evidence |
| --- | --- | --- |
| `GRHayLHD/par/Balsara0.par` | Shipped example | Time-terminated Balsara-labeled configuration; not a separate `test.ccl` case. |
| `GRHayLHD/test/Balsara0.par` | Authored test input | Iteration-10 input and five requested x-direction output fields. |
| `GRHayLHD/test/TOV.par` | Authored test input | Iteration-2 input and two requested HydroBase output fields. |
| `GRHayLHD/test/Balsara0/Balsara0.par` | Checked-in companion | Header asserts Cactus 4.15.0 generation and names an original Balsara0 input path. |
| `GRHayLHD/test/TOV/TOV.par` | Checked-in companion | Header asserts Cactus 4.15.0 generation and names an original TOV input path. |

Both authored test inputs set `IOUtil::parfile_write = "generate"`. This and
companion comments are checked-in observations, not proof of current authored-
input identity or oracle production.

### Eleven checked-in oracle files

| File | Visible dataset/direction | Header and numeric-row observations |
| --- | --- | --- |
| `Balsara0/press.x.asc` | `press`, x | CarpetIOASCII header names source parameter path and nine columns; 399 rows each at iterations 0/time 0 and 10/time 0.0125. |
| `Balsara0/rho.x.asc` | `rho`, x | Same header structure and block sizes/times. |
| `Balsara0/vx.x.asc` | `vx`, x | Same header structure and block sizes/times. |
| `Balsara0/vy.x.asc` | `vy`, x | Same header structure and block sizes/times. |
| `Balsara0/vz.x.asc` | `vz`, x | Same header structure and block sizes/times. |
| `TOV/hydrobase-press.x.asc` | HydroBase-labeled pressure, x | Only generic CarpetIOASCII header; numeric rows have nine fields: 168 at iteration 0/time 0 and 98 at iteration 2/time 0.125. |
| `TOV/hydrobase-press.y.asc` | HydroBase-labeled pressure, y | Same minimal header and row counts/times. |
| `TOV/hydrobase-press.z.asc` | HydroBase-labeled pressure, z | Same minimal header and row counts/times. |
| `TOV/hydrobase-rho.x.asc` | HydroBase-labeled density, x | Same minimal header and row counts/times. |
| `TOV/hydrobase-rho.y.asc` | HydroBase-labeled density, y | Same minimal header and row counts/times. |
| `TOV/hydrobase-rho.z.asc` | HydroBase-labeled density, z | Same minimal header and row counts/times. |

Balsara headers explicitly define fields as iteration, grid indices, time,
coordinates, and data. TOV files provide no comparable variable, parameter-
filename, or column-format header beyond generic CarpetIOASCII identification;
their filenames and nine-field numeric rows are recorded without assigning
undocumented column semantics.

### Provenance boundary

Balsara oracle headers name `repos/GRHayLET/GRHayLHD/test/Balsara0.par`, while
companion header names a different absolute original path. TOV companion names
an absolute original path, while TOV oracle headers omit parameter-path detail.
These assertions do not establish current-input identity or companion-to-
oracle production; see [GRH-0009](../contradictions.md#grh-0009) and
[GRH-0010](../contradictions.md#grh-0010).

## Caveats

- `RELTOL` is a declared tolerance; local files do not establish how an
  external harness applies it.
- Checked-in numeric content is an observation, not proof of correctness or
  current test status.
- TOV's Font1D assignment is configuration evidence only; selection and
  algorithm behavior are delegated to GRHayLib.
- No tabulated or explicitly entropy-selected authored input is visible.
- File naming and headers do not prove artifact chronology or lineage.

## Sources

- [Test declarations](../../../GRHayLHD/test/test.ccl)
- [Shipped example](../../../GRHayLHD/par/Balsara0.par)
- [Authored Balsara0 input](../../../GRHayLHD/test/Balsara0.par)
- [Authored TOV input](../../../GRHayLHD/test/TOV.par)
- [Balsara0 companion and oracles](../../../GRHayLHD/test/Balsara0)
- [TOV companion and oracles](../../../GRHayLHD/test/TOV)

## Related Pages

- [Parameters and Configurations](../integration/parameters-and-configurations.md)
- [EOS and Entropy Variants](../evolution/eos-entropy-variants.md)
- [Perturbations and Diagnostics](../evolution/perturbations-and-diagnostics.md)
