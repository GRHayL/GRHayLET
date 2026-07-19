# Coverage Gaps

> Page status: reviewed · Last reviewed: 07-19-2026
> Up: [Validation](index.md)

## Scope and Non-Scope

This page records regression dimensions not owned by checked-in GRHayLID test
declarations, parfiles, companions, or numeric oracles, then ranks proposed
future evidence. Checked-in absence means “not visibly covered here,” never
broken, unsupported, or untested elsewhere. No proposal was implemented or
executed.

## Summary

The complete 14-file source tree has no `test.ccl`, test/example parfile,
companion configuration, or numeric oracle. Consequently no checked-in
artifact observes the seven live one-dimensional hydro selections across
three directions, magnetic selection and staggering, either three-dimensional
setup, beta equilibrium, or entropy dispatch. The dead sound-wave dispatch and
the unconsumed staggering control are especially narrow targets for future
declaration and observation evidence.

## Mode Applicability

| Applicability | Visible checked-in gap |
| --- | --- |
| Common | No local test declaration, parfile, companion, oracle, build result, or runtime trace. |
| HydroTest1D | No case covers seven live arms across x/y/z; sound-wave selection has no dispatch evidence. |
| HydroTest1D+Magnetic | No case distinguishes magnetic on/off, either keyword selection, both output arrays, or staggering intent. |
| IsotropicGas | No checked-in setup or observation covers preconditions, sentinels, uniform state, or EOS-call boundary. |
| ConstantDensitySphere | No checked-in setup or observation covers interior/exterior branches, velocity controls, or EOS-call boundary. |
| BetaEquilibrium | No checked-in setup or observation covers guard, ordering, atmosphere threshold, writes, or missing entropy write. |
| Entropy/Hybrid | No checked-in setup or observation covers Hybrid schedule selection or entropy output. |
| Entropy/Tabulated | No checked-in setup or observation covers bounds-call dataflow, six declared writes, or non-Hybrid/Tabulated no-arm selection. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `VAL-GAP-01` | No checked-in configuration or regression declaration selects any of the eight declared `initial_data_1D` keyword values. | coverage-gap | Parameter exposes selection axis; complete registry has no par/test source class | `ccl:GRHayLID/param.ccl#parameter=initial_data_1D` |
| `VAL-GAP-02` | No checked-in artifact observes whether declared staggering intent changes magnetic data. | coverage-gap | Declared but visibly unconsumed control | `ccl:GRHayLID/param.ccl#parameter=stagger_A_fields` |
| `VAL-GAP-03` | No checked-in artifact observes the beta-equilibrium schedule, ordering, threshold branch, or declared writes. | coverage-gap | Beta-equilibrium schedule boundary | `ccl:GRHayLID/schedule.ccl#schedule=GRHayLID_BetaEquilibrium` |
| `VAL-GAP-04` | No checked-in artifact resolves the commented sound-wave dispatch arm against the later loop-side branch. | coverage-gap | Complete hydro function control flow | `c:GRHayLID/src/1D_tests_hydro_data.c#symbol=GRHayLID_1D_tests_hydro_data` |
| `VAL-GAP-05` | The build manifest lists six implementation units; the registered tree contains no test implementation or oracle artifact. | coverage-gap | Complete `SRCS` field plus registry audit | `build:GRHayLID/src/make.code.defn#field=SRCS` |
| `VAL-GAP-06` | Complete registry expansion contains no checked-in configuration supplying dispatch evidence for `initial_entropy="GRHayLID"`. | coverage-gap | Declared selection axis plus complete-registry absence | `ccl:GRHayLID/param.ccl#parameter=initial_entropy` |

## Details

### Checked-in evidence inventory

Registry expansion contains one README, one ThornGuide, four CCL files, one
build manifest, one header, and six C units. There is no path under
`GRHayLID/**` for `test.ccl`, a `.par` file, a test/example directory, a
companion configuration, or a numeric output/oracle. This absence can change
when local files are added; source registration and reconciliation must occur
before any coverage status changes.

### Untested axes visible in local declarations

| Area | Visible selection/dataflow axis lacking checked-in observation |
| --- | --- |
| One-dimensional hydro | Seven live dispatch arms—Balsara1-5, equilibrium, shock tube—times x/y/z; discontinuity position; EOS rejection wording; energy formula |
| Sound wave | Declared keyword, commented pre-loop arm, final error arm ordered before loop, sine expression, and acknowledged pressure/kinetic-energy incompleteness |
| One-dimensional magnetic | Magnetic schedule boolean; either/or Avec/Bvec selection; unconditional writes of both; five Balsara field-state pairs; rotations |
| Staggering | Default-enabled `stagger_A_fields`, no visible consumer, direct coordinate-array assignments with no visible offset arithmetic, and contrary source comment |
| IsotropicGas | Three compatibility checks, three sentinels, one external EOS call, uniform writes, zero velocity |
| ConstantDensitySphere | Compatibility checks, seven sentinels, unchecked velocity parameters, interior/exterior EOS calls, radius branch, duplicate zero writes |
| Beta equilibrium | Tabulated guard, bounds-before-sentinel order, setup/error calls, atmosphere threshold, four writes, and no entropy write |
| Entropy | Hybrid/Tabulated schedule arms, external calls and mutation boundary, README any-EOS wording, and no declared arm/error for other EOS values |

The list identifies local axes only. It makes no statement about external test
suites, downstream configurations, scientific adequacy, or defect presence.

### Ranked proposals

“Expected evidence” names a missing checked-in artifact, not a predicted pass.

| Rank / risk | Smallest proposed addition | Expected checked-in evidence | Main dependency |
| --- | --- | --- | --- |
| P0 / dispatch contradiction | Minimal sound-wave parfile and declared regression case | Selection record plus bounded density/pressure/velocity observations distinguishing pre-loop dispatch from loop formula | Any source correction needed to define intended control flow; authorized runtime environment |
| P0 / declaration-consumer drift | Paired magnetic cases with `stagger_A_fields=yes/no` and explicit magnetic keyword selections | Avec/Bvec observations at coordinates where a half-cell offset would be distinguishable | Defined staggering convention and source consumer; HydroBase output support |
| P0 / lifecycle ambiguity | Entropy cases for Hybrid, Tabulated, and one externally admitted other EOS value | Schedule/runtime record and entropy-related fields showing which local arm is selected | External EOS domain and initialized GRHayLib; possible explicit fallback policy |
| P0 / write-set mismatch | Beta-equilibrium case near and above atmosphere threshold | Electron fraction, temperature, pressure, energy, and entropy observations plus bounded diagnostic output | Tabulated EOS data and external beta-equilibrium helper semantics |
| P1 / 1D matrix | One declared case per seven live arm, with x/y/z rotation coverage distributed across cases | Hydrodynamic fields and magnetic fields where applicable, with exact tolerances | Agreed reference data and output variables |
| P1 / 3D setups | Minimal IsotropicGas and ConstantDensitySphere parfiles and cases | Uniform gas fields; sphere interior/exterior fields and nonzero interior velocity probe | Tabulated EOS data and small grid configuration |
| P2 / build boundary | Supported build configuration with declared HDF5 dependency | Reproducible build-system record scoped to revision and environment | Explicit build authority and installed external dependencies |

### Issue-linked targets

[GID-0004](../contradictions.md#gid-0004) needs staggering-consumer evidence.
[GID-0005](../contradictions.md#gid-0005) and
[GID-0006](../contradictions.md#gid-0006) need sound-wave dispatch and state
evidence. [GID-0010](../contradictions.md#gid-0010) needs explicit dispatch or
error policy for a selected entropy keyword outside the two named EOS arms.
None is resolved by a proposed test alone; resolution requires admissible
changed local evidence and full backlink reconciliation.

## Caveats

- Missing local artifacts do not prove missing external tests or invalid code.
- A future parfile proves checked-in selection only; a test declaration proves
  declaration, and an oracle proves only stored content in its fixture context.
- Runtime reproduction needs explicit authority, isolated configuration,
  recorded revision/command/platform/input/result, and bounded resources.
- No Cactus build, test, external dependency, or GRHayLID routine was run in
  this static review.

## Sources

- [Parameter declarations](../../../GRHayLID/param.ccl)
- [Schedule declarations](../../../GRHayLID/schedule.ccl)
- [One-dimensional hydro implementation](../../../GRHayLID/src/1D_tests_hydro_data.c)
- [Build manifest](../../../GRHayLID/src/make.code.defn)

## Related Pages

- [One-D Hydro Tests](../initial-data/one-d-tests-hydro.md)
- [One-D Magnetic Tests](../initial-data/one-d-tests-magnetic.md)
- [Beta Equilibrium](../initial-data/beta-equilibrium.md)
- [Entropy Computation](../initial-data/entropy-computation.md)
- [Parameters and Configurations](../integration/parameters-and-configurations.md)
- [GID-0004: unconsumed staggering control](../contradictions.md#gid-0004)
- [GID-0005: sound-wave dispatch mismatch](../contradictions.md#gid-0005)
- [GID-0006: sound-wave state hazard](../contradictions.md#gid-0006)
- [GID-0010: entropy no-arm ambiguity](../contradictions.md#gid-0010)
