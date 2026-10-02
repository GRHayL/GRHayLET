## GRHayLID issue dispositions

All 20 ranked groups have local repairs or explicit compatibility contracts as
of 10-02-2026. The historical audit is retained below; its original line
numbers and GRHayL-main paths refer to the audit environment, not this checkout.
The referenced GRHayLID_fps.md companion is absent here. External-library
findings were rechecked against the public GRHayL sources before replacing the
local beta call path; no external dependency was modified.

| Issue | Local disposition | Implementation / contract |
| --- | --- | --- |
| 1 | Repaired | Gas and both sphere regions validate finite effective-bounds inputs, check EOS status/output, and reject failures before grid writes. |
| 2 | Repaired | ParamCheck and magnetic scheduling require both GRHayLID selectors; body checks both destination allocations. Mixed ownership is rejected. |
| 3 | Repaired with input contract | Beta explicitly requests Ye/T storage; standalone routines check actual storage. Entropy requires producers initialized in HydroBase_Initial; missing selectors are diagnosed without beta. |
| 4 | Unsupported precision rejected | Header requires CCTK_REAL_PRECISION_8, bounding the GRHayLib parameter-array interface too; EOS output calls also stage doubles. |
| 5 | Repaired | Hybrid uses the cold-plus-thermal helper; below-cold-curve benchmark states are rejected. Simple retains its gamma-law energy. |
| 6 | Repaired | Live sound-wave dispatch initializes all pre-loop state, uses the selected velocity axis, and requires finite amplitude below unity. README identifies it as a velocity perturbation. |
| 7 | Repaired | Simple joins the native entropy schedule arm, subject to the explicit proxy contract in Issue 14. |
| 8 | Repaired | ParamCheck rejects magnetic and Ye/T selectors without a matching local producer and validates the hydro/EOS matrix. |
| 9 | Local call path replaced | Beta solves at the actual requested temperature using interpolated base potentials, then verifies the final residual. The external cached-root builder remains unchanged. |
| 10 | Local call path replaced | Local search accepts node/tolerance-qualified roots and rejects no-root/nonconvergent cases; it never silently reports minimum Ye as success. |
| 11 | Repaired | Atmosphere resets rho with Ye/T and recomputes pressure/epsilon at that complete tuple; rho is declared WRITES. Atmosphere is explicitly exempt from equilibrium. |
| 12 | Repaired | Beta uses effective bounds. Post-beta entropy preserves rho/Ye/T and rejects incompatible states instead of clamping them. |
| 13 | Repaired | Sphere checks finite combined Cartesian speed strictly below one; native tests check the identity spatial metric. |
| 14 | Explicit compatibility contract | Simple/Hybrid proxy output requires allow_native_entropy_proxy=yes (default no). Generic physical-entropy consumers must use a compatible physical producer. Native helper semantics are preserved. |
| 15 | Restricted geometry enforced | Native test bodies require an ADMBase identity spatial metric; README requires Cartesian coordinates/components. Curved-background construction is unsupported. |
| 16 | Local call path replaced | Beta uses mu_e-mu_n+mu_p rather than shifted munu. Correct base-potential conventions/provenance remain a table-owner obligation. |
| 17 | Documented normalization contract | README and ThornGuide require IllinoisGRMHD::rescale_magnetics=no for the normalized Balsara data; legacy converter defaults are preserved. |
| 18 | Repaired | Native entropy rejects zero/negative/nonfinite density and negative/nonfinite pressure, and checks finite output before publication. |
| 19 | Repaired | Avec receives component-specific transverse half-cell offsets when requested; Bvec stays at base coordinates. |
| 20 | Repaired | Guide, README, selector help and schedule descriptions match actual controls, outputs, capabilities, entropy representation and primitive side effects. |

Changed implementation and declaration owners are
[GRHayLID/src](GRHayLID/src), [schedule.ccl](GRHayLID/schedule.ccl),
[param.ccl](GRHayLID/param.ccl), [README](GRHayLID/README), and
[ThornGuide](GRHayLID/doc/documentation.tex). The
[GRHayLID KB](wiki/grhaylid/index.md), canonical source-map dependencies, and
historical GID issues are reconciled with the changed local sources.

Verification uses disposable actual-source probes outside the repository,
compiled with GCC 13.3, C11, strict warnings, AddressSanitizer and
UndefinedBehaviorSanitizer. Cactus argument/storage/error facilities are stubs;
selected GRHayL EOS interpolation, bounds and Hybrid routines are actual source.
Synthetic tables independently specify temperature/density-dependent roots,
exact node roots, absent roots, and a shifted munu dataset. The probes cover
EOS rejection, both sphere regions, Hybrid unequal exponents/integration
constants, all sound-wave directions, selector/storage failures, entropy
input/representation guards, geometry and component staggering, final beta
residuals and post-beta entropy. The baseline gas probe fails the expected
invalid-input rejection; the repaired probe suite passes. Alternate-real
syntax probes require an early unsupported-precision diagnostic.

The default beta residual tolerance is an absolute numerical criterion of
1e-8 MeV, configurable by beq_residual_tolerance; it is not an established
production physical-error tolerance. No configured Cactus build/scheduler run,
production EOS-table comparison, or current Einstein Toolkit regression pass
is claimed. These remain validation gaps, not newly discovered code defects.
G01 Lorentz-factor ownership and G02 legacy runtime compatibility remain
consumer/configuration obligations; G03 retains the absent local ET regression
suite. G04's density-node-root interpolation mechanism is bypassed locally by
the final-state solve, with production residual validation still required.

## Historical audit record

The following 20 groups retain the original mechanism, trigger and proposed
verification for traceability. Their wording describes the pre-repair audit
state. High/medium/low correspond to P1/P2/P3; prevalence was not measured.

<details>
<summary><strong>Issue 1 (high):</strong> Gas and sphere initializers consume uninitialized outputs after tabulated EOS errors.</summary>

<br>

Original finding: `F01`.

- **Description**: Gas and sphere initializers declare pressure and internal-energy locals without initialization, discard the tabulated EOS return status, and copy those locals into HydroBase. The current EOS implementation can return an error before writing either output.

- **Production trigger and likelihood**: Zero density or temperature, or another out-of-table primitive triple, passes the local sentinel checks and allowed parameter ranges. This is a reachable malformed-input path; valid in-table inputs do not establish the defect.

- **Impact**: Indeterminate pressure and energy are read and published instead of the input error being handled. Later native primitive repair can limit persistence, but does not make this initializer error path valid.

- **Locations**: [GRHayLID/src/IsotropicGas.c](GRHayLID/src/IsotropicGas.c):29–36,47–51; [GRHayLID/src/ConstantDensitySphere.c](GRHayLID/src/ConstantDensitySphere.c):33–49,65–83; [GRHayLID/param.ccl](GRHayLID/param.ccl):83–99,110–126,143–159; [NRPyEOS_P_and_eps_from_rho_Ye_T.c](GRHayL-main/GRHayL/EOS/Tabulated/interpolators/NRPyEOS_P_and_eps_from_rho_Ye_T.c):19–29.

- **Repository evidence and checks**: Actual-source probes observed untouched EOS outputs after failure. A historical MemorySanitizer gas probe detected uninitialized-value use; its exact compiler command was not retained. Source tracing independently establishes the mechanism. Compiler-filled diagnostic patterns in another probe do not predict production values. Tabulated entropy's ordinary finite-input clamp protects its interpolation domain; its ignored status alone was not established as a separate ordinary failure.

- **Confidence**: High for the unchecked-error mechanism and stated input boundary.

- **Desired postcondition and fix direction**: Check every EOS return before entering the grid loop. Reject the requested state with its region and primitive triple, or deliberately clamp and recompute a consistent state. Initializing the output locals alone does not resolve the failed calculation.

- **Verification**: Exercise valid and invalid gas inputs and both sphere regions with actual EOS calls under a compatible memory checker. Require valid initialized outputs on success and the intended rejection or consistent recovery on error.

</details>

<details>
<summary><strong>Issue 2 (high, conditional storage):</strong> Selecting either magnetic destination admits a routine that writes both.</summary>

<br>

Original finding: `F02`.

- **Description**: The 1D magnetic routine accepts either GRHayLID Avec or Bvec selection, but unconditionally writes both groups. Its schedule also declares both writes without establishing both allocations or reconciling competing owners.

- **Production trigger and likelihood**: Select one GRHayLID destination and leave the other at none with no alternative storage provider. IllinoisGRMHD allocates both groups and masks this missing-storage case; shipped Balsara configurations select both. Selecting another initializer for the second group instead creates competing writes without an explicit local order.

- **Impact**: The admitted routine accesses an unavailable destination, or the driver rejects the declaration before execution. With another producer, one initializer can overwrite the other's output. A universal IllinoisGRMHD crash is not claimed.

- **Locations**: [GRHayLID/src/1D_tests_magnetic_data.c](GRHayLID/src/1D_tests_magnetic_data.c):8–9,81–89,101–128; [GRHayLID/schedule.ccl](GRHayLID/schedule.ccl):12–18; [IllinoisGRMHD/schedule.ccl](IllinoisGRMHD/schedule.ccl):3–4; official [HydroBase schedule](https://bitbucket.org/einsteintoolkit/einsteinbase/raw/master/HydroBase/schedule.ccl):16–23.

- **Repository evidence and checks**: Source and HydroBase declarations establish independent destination allocation. Actual-body probes faulted with either missing buffer represented as NULL. These stubs do not establish every driver's buffer representation or schedule-check behavior.

- **Confidence**: High for the selection/write mismatch; runtime consequences remain conditional on storage and driver behavior.

- **Desired postcondition and fix direction**: Write and declare only owned destinations, or require both selectors and diagnose mismatches before initialization. Establish ordering if another producer intentionally shares a destination.

- **Verification**: Check both single-selector cases without an alternative allocator, both-selector initialization, and mixed-producer selection in a real ET schedule. Require a clear precondition failure or correctly limited writes without unavailable accesses or unordered competing initialization.

</details>

<details>
<summary><strong>Issue 3 (high, conditional storage):</strong> Standalone beta and Tabulated entropy omit optional-field storage preconditions.</summary>

<br>

Original finding: `F03`.

- **Description**: Beta equilibrium is scheduled from its boolean option and writes Ye and temperature without checking their storage. Tabulated entropy similarly reads and mutates those optional groups without establishing allocation or initialized input ownership.

- **Production trigger and likelihood**: An external hydro initializer is combined with either standalone operation while a required Ye/T selector remains none and no other storage provider exists. HydroBase defaults these optional selectors to none. External-ID use alone is not evidence of missing storage.

- **Impact**: Missing read/write buffers can be accessed. Even allocated fields are unsafe entropy inputs if no earlier producer initialized them.

- **Locations**: [GRHayLID/schedule.ccl](GRHayLID/schedule.ccl):41–48,59–66; [GRHayLID/src/BetaEquilibrium.c](GRHayLID/src/BetaEquilibrium.c):10–21,31–34,50–53; [GRHayLID/src/ComputeEntropy.c](GRHayLID/src/ComputeEntropy.c):28–31; official [HydroBase storage](https://bitbucket.org/einsteintoolkit/einsteinbase/raw/master/HydroBase/schedule.ccl):8–10,28–31 and [selector defaults](https://bitbucket.org/einsteintoolkit/einsteinbase/raw/master/HydroBase/param.ccl).

- **Repository evidence and checks**: Source tracing establishes the standalone schedule and accesses without local allocation checks. A beta actual-body probe with absent buffers faulted. No complete-driver missing-storage reproduction was performed.

- **Confidence**: High for the missing local prerequisites, with the stated storage/provider condition.

- **Desired postcondition and fix direction**: Require or establish active Ye/T storage independently of the hydro family. Document and enforce the initialized inputs required by entropy; storage allocation alone is insufficient.

- **Verification**: Exercise external-ID beta and entropy with each required group absent, then with allocated and initialized inputs. Require clear diagnostics or established storage, and verify that entropy cannot consume fields before their producer.

</details>

<details>
<summary><strong>Issue 4 (high, alternate precision):</strong> CCTK_REAL output pointers are passed to double-pointer EOS APIs.</summary>

<br>

Original finding: `F14`.

- **Description**: Gas/sphere EOS outputs use CCTK_REAL locals, and Tabulated entropy passes CCTK_REAL grid-element addresses directly. The supplied GRHayL API requires double pointers; there is no representation bridge for a non-double CCTK_REAL.

- **Production trigger and likelihood**: An attempted Cactus configuration with a non-double real representation. Cactus lists alternate real sizes subject to architecture support. Normal eight-byte double precision is unaffected; no supported alternate-precision ET run was established.

- **Impact**: Incompatible pointer diagnostics or incorrectly represented writes. If a four-byte call is compiled through the diagnostics, a double output can overwrite adjacent storage; a distinct sixteen-byte representation is also incompatible.

- **Locations**: [GRHayLID/src/IsotropicGas.c](GRHayLID/src/IsotropicGas.c):29–36; [GRHayLID/src/ConstantDensitySphere.c](GRHayLID/src/ConstantDensitySphere.c):33–49; [GRHayLID/src/ComputeEntropy.c](GRHayLID/src/ComputeEntropy.c):28–31; [ghl_eos_functions.h](GRHayL-main/GRHayL/include/ghl_eos_functions.h):91–106,338–342; [GRHayLib initialization](GRHayL-main/implementations/GRHayLib/src/initialize_and_shutdown.c):224–228.

- **Repository evidence and checks**: Actual signatures and [Cactus precision options](https://www.cactuscode.org/documentation/usersguide/UsersGuidech6.html) establish the conditional mismatch. A historical four-byte syntax probe rejected the translation units; its exact command was not retained. No actual alternate-precision Cactus build was run. GRHayLib's parameter-array boundary also requires adaptation.

- **Confidence**: High for incompatibility when the representations differ; supported-build prevalence is unestablished.

- **Desired postcondition and fix direction**: Declare and enforce double precision, or stage through double temporaries and explicit value conversions across all affected boundaries. Pointer casts do not repair representation differences.

- **Verification**: Confirm the normal double build, then test any deliberately supported alternate representation through complete GRHayLib/GRHayLID initialization. Require either an early supported-precision diagnostic or correct bounded output writes and conversions.

</details>

<details>
<summary><strong>Issue 5 (medium):</strong> General Hybrid initial data uses the ideal-fluid internal-energy relation.</summary>

<br>

Original finding: `F04`.

- **Description**: The 1D initializer accepts Hybrid but computes epsilon as P/[rho*(Gamma_piece−1)]. The canonical Hybrid relation includes cold energy, piece integration constants, and the thermal Gamma: epsilon_cold+(P−P_cold)/[rho*(Gamma_th−1)].

- **Production trigger and likelihood**: An accepted Hybrid EOS whose thermal Gamma differs from the cold piece's Gamma or whose active piece has a nonzero energy integration constant. Matching single-piece gamma-law settings and the shipped Simple examples can avoid the mismatch.

- **Impact**: The published HydroBase tuple is inconsistent with its selected EOS. Native limiters can recompute epsilon later, so persistence into every evolution consumer is not established.

- **Locations**: [GRHayLID/src/1D_tests_hydro_data.c](GRHayLID/src/1D_tests_hydro_data.c):125–128; [NRPyEOS_hybrid_compute_epsilon.c](GRHayL-main/GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_epsilon.c):20–26; [NRPyEOS_compute_P_cold_and_eps_cold.c](GRHayL-main/GRHayL/EOS/Hybrid/NRPyEOS_compute_P_cold_and_eps_cold.c):22–31.

- **Repository evidence and checks**: For rho=P=1, Gamma_piece=2, Gamma_th=1.5 and K=.1, actual-source calculations give initializer epsilon=1 versus canonical epsilon=1.9. The accepted parameter family and source formulas establish the discrepancy without invalid inputs.

- **Confidence**: High for general accepted Hybrid settings at the initializer boundary.

- **Desired postcondition and fix direction**: Use the existing Hybrid energy helper, retaining the simple relation for Simple, or restrict the benchmarks explicitly to a matching gamma-law EOS. Define handling below the admissible cold pressure curve.

- **Verification**: Compare initialized energy to an independently evaluated cold-plus-thermal relation for unequal exponents and a multi-piece integration constant, then retain matching Simple/gamma-law coverage before downstream repair.

</details>

<details>
<summary><strong>Issue 6 (medium):</strong> The advertised sound-wave selector terminates before reaching its implementation.</summary>

<br>

Original finding: `F05`.

- **Description**: The parameter permits sound wave, but its pre-loop dispatch arm is commented out. Selection falls into CCTK_VERROR before the implemented wave loop. A y/z velocity-direction problem is latent in that unreachable loop.

- **Production trigger and likelihood**: Request sound wave with otherwise supported Simple/Hybrid data. The dispatch failure follows directly from this accepted selector; the orientation problem is exposed only if dispatch is restored.

- **Impact**: Initialization aborts for an advertised setup. Restoring dispatch alone would still put the longitudinal perturbation in vel[x] for y/z waves and can expose uninitialized pre-loop left/right state.

- **Locations**: [GRHayLID/param.ccl](GRHayLID/param.ccl):51–61; [GRHayLID/src/1D_tests_hydro_data.c](GRHayLID/src/1D_tests_hydro_data.c):55–59,67–69,99–111.

- **Repository evidence and checks**: Full dispatch tracing and the [Cactus terminating-error contract](https://einsteintoolkit.org/referencemanual/ReferenceManual.html) establish the abort. An actual-body stub probe observed that route. No executed y/z sound-wave failure is claimed for the current fatal dispatch.

- **Confidence**: High for the live dispatch failure; the orientation finding is explicitly latent.

- **Desired postcondition and fix direction**: Restore complete initialization and direction-aware wave state, with defined admissible amplitude, or remove the unsupported selector. Avoid uninitialized rotations; kinetic energy is not required in thermal pressure.

- **Verification**: Exercise x/y/z sound waves through the live dispatch, check finite EOS-consistent states and perturbation direction, and verify admissible amplitude bounds. If removed, require explicit rejection rather than an advertised broken option.

</details>

<details>
<summary><strong>Issue 7 (medium):</strong> Simple EOS entropy selection schedules no requested producer.</summary>

<br>

Original finding: `F06`.

- **Description**: GRHayLID's entropy schedule branches only for Hybrid and Tabulated, although GRHayLib supports Simple and initializes compatible native entropy machinery for it.

- **Production trigger and likelihood**: EOS_type=Simple and HydroBase::initial_entropy=GRHayLID with no independent entropy producer. The selection is declared and the destination can be allocated without a writer; later native recomputation can mask the missing initialization.

- **Impact**: Requested initial entropy remains uninitialized without an unsupported-selection diagnostic. The any-EOS capability description exceeds the schedule implementation.

- **Locations**: [GRHayLID/schedule.ccl](GRHayLID/schedule.ccl):51–68; [GRHayLib/param.ccl](GRHayL-main/implementations/GRHayLib/param.ccl):73–79; [Simple EOS initialization](GRHayL-main/GRHayL/GRHayL_Core/initialize_eos.c):123–155; [GRHayLID/README](GRHayLID/README):26–30.

- **Repository evidence and checks**: Schedule, Simple metadata and actual EOS dispatch were compared directly. Source establishes the omitted producer; no inevitable native-converter failure is inferred.

- **Confidence**: High for the selected-producer omission.

- **Desired postcondition and fix direction**: Dispatch Simple to the compatible routine, or reject the selection explicitly and narrow the advertised capability.

- **Verification**: Select Simple entropy without another writer and inspect the field before downstream conversion. Require the documented initialized representation or an explicit unsupported-selection error.

</details>

<details>
<summary><strong>Issue 8 (medium, configuration):</strong> Other GRHayLID field selectors can name an inactive producer.</summary>

<br>

Original finding: `F07`.

- **Description**: Ye, temperature, Avec and Bvec selectors advertise GRHayLID independently of the schedule branches that can produce them. No local parameter check validates the complete family/option/selector matrix.

- **Production trigger and likelihood**: Select GRHayLID magnetics with gas/sphere, or disable magnetic initialization while retaining its magnetic selectors; alternatively select its Ye/T with HydroTest1D and beta off. A separately ordered, initialized external owner can avoid the omission.

- **Impact**: HydroBase can allocate requested fields while no selected GRHayLID writer initializes them. This is a missing-producer configuration defect, separate from unavailable-storage accesses in Issues 2 and 3.

- **Locations**: [GRHayLID/param.ccl](GRHayLID/param.ccl):14–37; [GRHayLID/schedule.ccl](GRHayLID/schedule.ccl):3–48; [GRHayLID/src/1D_tests_hydro_data.c](GRHayLID/src/1D_tests_hydro_data.c).

- **Repository evidence and checks**: The selector-to-producer matrix was traced against all local setup branches. Only the enabled HydroTest1D branch supplies magnetics; the 1D hydro body does not write Ye/T. Storage declarations do not prove initialized data.

- **Confidence**: High for the locally admitted inactive-producer combinations.

- **Desired postcondition and fix direction**: Diagnose selector combinations without an available initialized owner, or provide the selected producer with explicit ordering. Document external ownership where deliberately allowed.

- **Verification**: Exercise the named incompatible combinations and valid counterparts. Require every selected field to have a proven producer before use or a clear configuration error.

</details>

<details>
<summary><strong>Issue 9 (medium, upstream):</strong> Beta composition and thermodynamics use different temperatures.</summary>

<br>

Original finding: `F08`.

- **Description**: GRHayLID requests beta equilibrium at beq_temperature, writes that temperature, and queries thermodynamics at it. The supplied beta builder instead truncates the log-temperature index and solves using one lower table slice without temperature interpolation.

- **Production trigger and likelihood**: Requested T lies between nodes and equilibrium composition depends on T. Exact-node or temperature-independent cases can avoid this mechanism; production-table prevalence was not measured.

- **Impact**: The initialized composition need not be in equilibrium at the written temperature. This is not a consistently snapped state.

- **Locations**: [GRHayLID/src/BetaEquilibrium.c](GRHayLID/src/BetaEquilibrium.c):19–21,36–53; [NRPyEOS_tabulated_get_index.c](GRHayL-main/GRHayL/EOS/Tabulated/NRPyEOS_tabulated_get_index.c):9–16; [beta builder](GRHayL-main/GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c):235–257.

- **Repository evidence and checks**: An actual-source synthetic table with T={1,4} and equilibrium Ye={.2,.4} returned Ye=.2 at T=2, while the table-consistent root is .3 and the returned residual was −.1. This diagnoses the mechanism, not a production physical error magnitude. The caller's argument order is correct.

- **Confidence**: High for the upstream slice mismatch under the stated temperature dependence.

- **Desired postcondition and fix direction**: Interpolate chemical potential at requested T before solving. A restricted workaround may require verified T nodes or consistently snap and document the written T; it does not resolve Issue 10 or gap G04.

- **Verification**: Use a temperature-dependent independent root oracle at node and off-node T, checking the residual at the final written state against a specified tolerance. Separately retain exact-root and off-density-node checks.

</details>

<details>
<summary><strong>Issue 10 (medium, upstream):</strong> Exact roots and absent brackets can return minimum Ye as successful equilibrium.</summary>

<br>

Original finding: `F09`.

- **Description**: The root helper accepts only strictly negative neighboring chemical-potential products. It misses isolated exact node roots; when no strict bracket is found, it returns minimum Ye while the public builder reports success.

- **Production trigger and likelihood**: An exact root at an interior Ye node, or a slice without an in-range root. A root already at minimum Ye or another valid strict bracket can avoid a wrong result. No production-table frequency was established.

- **Impact**: Non-equilibrium composition can be reported as solved equilibrium, and GRHayLID's return checks cannot detect the falsely successful result.

- **Locations**: [find_Ye_st_munu_is_zero and public beta builder](GRHayL-main/GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c):7–31,253–261; [GRHayLID/src/BetaEquilibrium.c](GRHayLID/src/BetaEquilibrium.c):19–21.

- **Repository evidence and checks**: Actual-source probes returned minimum Ye with success and nonzero residuals for both interior exact-root and no-root cases. The failure is possible at exact T nodes, independently of Issue 9. An intentional minimum-Ye clipping policy does not prove physical equilibrium.

- **Confidence**: High for the root-exclusion and indistinguishable-fallback mechanisms.

- **Desired postcondition and fix direction**: Recognize exact or tolerance-qualified node roots and return an explicit no-bracket outcome. If clipping is intended, expose a distinct policy/status rather than calling it solved equilibrium.

- **Verification**: Check minimum/interior/maximum node roots, ordinary strict brackets, and same-sign slices against an independent oracle. Require correct roots or distinguishable no-equilibrium/fallback status.

</details>

<details>
<summary><strong>Issue 11 (medium):</strong> Beta's atmosphere branch preserves density inconsistent with its thermodynamics.</summary>

<br>

Original finding: `F10`.

- **Description**: For rho≤1.01*rho_atm, beta equilibrium retains input density while copying pressure, epsilon, Ye and temperature computed for rho_atm. The schedule declares density read-only.

- **Production trigger and likelihood**: Retained rho differs from rho_atm for a density-dependent EOS. Positive in-table densities can establish the discrepancy; vacuum or below-table density can also remain. Later repair depends on the consumer.

- **Impact**: The published atmosphere tuple can violate EOS closure at its stored density.

- **Locations**: [GRHayLID/src/BetaEquilibrium.c](GRHayLID/src/BetaEquilibrium.c):28–34; [GRHayLID/schedule.ccl](GRHayLID/schedule.ccl):45–47; [atmosphere EOS construction](GRHayL-main/GRHayL/GRHayL_Core/initialize_eos.c):446–457.

- **Repository evidence and checks**: An actual-body synthetic P=rho*T probe retained rho=.15 but wrote P=.2 at T=1 for rho_atm=.2. Source traces the copied values to the atmosphere EOS evaluation; downstream repair does not disprove the mixed tuple.

- **Confidence**: High for the density/thermodynamics mismatch.

- **Desired postcondition and fix direction**: Reset the complete atmosphere state, including density and its WRITES declaration, or evaluate thermodynamics at the retained valid density. Document any deliberate atmosphere exception to equilibrium.

- **Verification**: Exercise below, at and just above rho_atm within the branch, checking EOS closure at the stored density before conversion and confirming accurate read/write declarations.

</details>

<details>
<summary><strong>Issue 12 (medium, combined selections):</strong> Entropy clipping can undo imposed beta equilibrium.</summary>

<br>

Original finding: `F11`.

- **Description**: Beta admits temperatures within raw table bounds, then Tabulated entropy runs afterward and clamps rho/Ye/T to configured effective bounds. It replaces pressure and epsilon without solving equilibrium again.

- **Production trigger and likelihood**: Enable both operations when a correct equilibrium state lies outside narrower configured bounds. These effective bounds are supported; the mechanism can occur at an exact table-temperature node independently of Issues 9 and 10.

- **Impact**: The final initialized state can cease to satisfy beta equilibrium despite the preceding successful equilibrium operation.

- **Locations**: [GRHayLID/src/BetaEquilibrium.c](GRHayLID/src/BetaEquilibrium.c):13–15; [GRHayLID/schedule.ccl](GRHayLID/schedule.ccl):60–66; [GRHayLID/src/ComputeEntropy.c](GRHayLID/src/ComputeEntropy.c):28–31; [NRPyEOS_enforce_table_bounds.c](GRHayL-main/GRHayL/EOS/Tabulated/NRPyEOS_enforce_table_bounds.c):3–16.

- **Repository evidence and checks**: An actual-source probe began at correct T=4, Ye=.4. Configured T_max=2 changed T to 2 while retaining Ye=.4, yielding chemical-potential residual .1. Source tracing establishes valid narrower configured bounds and the post-beta mutation order.

- **Confidence**: High for the clipping/equilibrium interaction under simultaneous selections.

- **Desired postcondition and fix direction**: Agree on effective limits before constructing equilibrium. Reject incompatible requests or resolve equilibrium at the final state, and document entropy's primitive mutations.

- **Verification**: Combine beta and entropy with narrower temperature, density and Ye bounds using independent residual checks at the final stored state. Require the specified equilibrium/atmosphere policy or an explicit incompatibility diagnostic.

</details>

<details>
<summary><strong>Issue 13 (medium, input boundary):</strong> Sphere velocity parameters admit non-timelike initialized data.</summary>

<br>

Original finding: `F12`.

- **Description**: Sphere velocity components have unrestricted parameter ranges and are copied without a physical-speed check. Eulerian velocity requires gamma_ij*v^i*v^j<1.

- **Production trigger and likelihood**: On the documented Cartesian Minkowski background, permitted inputs such as vx=1.2, vy=vz=0 violate that condition. This is a malformed-input boundary, not a failure of ordinary subluminal data.

- **Impact**: Published velocity has no finite real Lorentz factor. Native speed limiting can repair it later; propagation into every consumer is not established.

- **Locations**: [GRHayLID/param.ccl](GRHayLID/param.ccl):128–141; [GRHayLID/src/ConstantDensitySphere.c](GRHayLID/src/ConstantDensitySphere.c):7–9,74–83.

- **Repository evidence and checks**: The actual initializer accepted and wrote vx=1.2, giving v²=1.44 in Minkowski space. The [GRHydro Lorentz-factor definition](https://einsteintoolkit.org/thornguide/EinsteinEvolve/GRHydro/documentation.html) supplies the admissibility relation. No universal downstream NaN result was demonstrated.

- **Confidence**: High for the local admissibility omission on the stated background.

- **Desired postcondition and fix direction**: Reject nonfinite or luminal/superluminal combinations for the supported background, or validate against the actual metric if extending geometry support.

- **Verification**: Check ordinary subluminal inputs, combined-component speeds at and above unity, and nonfinite inputs. Require finite admissible data or a clear rejection before publication.

</details>

<details>
<summary><strong>Issue 14 (medium, interoperability):</strong> Hybrid entropy publishes a native proxy under HydroBase's physical-specific-entropy interface.</summary>

<br>

Original finding: `F15`.

- **Description**: GRHayLID stores the native Hybrid quantity P/rho^(Gamma−1) directly in HydroBase::entropy. HydroBase defines that field as physical specific entropy in k_b/baryon.

- **Production trigger and likelihood**: A generic consumer interprets the initialized field according to HydroBase's physical entropy contract. Native GRHayL evolution/recovery intentionally expects the proxy and is outside this failure allegation.

- **Impact**: Matching EOS settings alone does not align entropy meanings. Along P=K*rho^Gamma, the stored proxy is K*rho rather than constant physical specific entropy on an isentrope.

- **Locations**: [GRHayLID/src/ComputeEntropy.c](GRHayLID/src/ComputeEntropy.c):12; [NRPyEOS_compute_entropy_function.c](GRHayL-main/GRHayL/EOS/Hybrid/NRPyEOS_compute_entropy_function.c):22–30; official [HydroBase entropy definition](https://einsteintoolkit.org/thornguide/EinsteinBase/HydroBase/documentation.html).

- **Repository evidence and checks**: The actual helper formula, native consumer representation and primary interface were compared directly. Isentropic scaling establishes the semantic difference. No failing generic runtime consumer or native entropy-recovery failure was reproduced.

- **Confidence**: High for the representation mismatch under the named consumer interpretation.

- **Desired postcondition and fix direction**: Document and constrain interoperability, or specify a physical conversion and separately owned native-proxy field. Preserve native contracts through an explicit migration; replacing the canonical helper directly is insufficient.

- **Verification**: Check the agreed representation and units on an isentropic family and trace both physical-entropy and native-proxy consumers. Require each to receive its documented quantity without breaking native recovery.

</details>

<details>
<summary><strong>Issue 15 (medium, geometry compatibility):</strong> Magnetic initial data depends on flat Cartesian background and component conventions.</summary>

<br>

Original finding: `F16`.

- **Description**: GRHayLID directly writes B and builds Cartesian linear A without reading the metric. IllinoisGRMHD reconstructs densitized curl fields from A and divides by sqrt(det gamma) to obtain physical B.

- **Production trigger and likelihood**: Compose these data with sqrt(det gamma)≠1 or another component basis without the required adaptation. Standard Cartesian Minkowski Balsara configurations avoid this condition; curved-background benchmark support was not established.

- **Impact**: Direct B and curl-derived physical B can disagree. For sqrt(det gamma)=2 and curl A=1, reconstructed physical B=.5 while direct B is 1.

- **Locations**: [GRHayLID/src/1D_tests_magnetic_data.c](GRHayLID/src/1D_tests_magnetic_data.c):81–128; [IllinoisGRMHD/src/compute_B_and_Bstagger_from_A.c](IllinoisGRMHD/src/compute_B_and_Bstagger_from_A.c):53–79,102–138.

- **Repository evidence and checks**: The initializer and actual metric-densitized reconstruction were traced, including the determinant factor. This establishes a compatibility precondition, not a demonstrated failure of supplied flat-space tests or general discrete boundary/stencil equality.

- **Confidence**: High for the stated geometry/component precondition and determinant mismatch.

- **Desired postcondition and fix direction**: State and enforce the intended background and component basis, identifying the external ADMBase owner, or implement a specified metric-aware construction and transformation.

- **Verification**: Verify the documented flat-background convention and the determinant-factor example independently. If broader geometry is supported, compare direct and reconstructed physical components with the correct metric/basis transformations and boundary conventions.

</details>

<details>
<summary><strong>Issue 16 (medium, affected upstream tables):</strong> Original shifted LS-table munu can produce the wrong physical equilibrium root.</summary>

<br>

Original finding: `F17`.

- **Description**: The StellarCollapse distribution warns that affected original LS tables contain shifted munu in the LS density region above 1e8 g/cm³. Current GRHayL ingestion copies that dataset unchanged, and beta equilibrium solves its zero rather than rebuilding the physical potential combination.

- **Production trigger and likelihood**: Use an affected original table in the affected density region. Repaired LS tables and other EOS tables are excluded; no affected production table was downloaded or exercised during the audit.

- **Impact**: The helper can solve a shifted zero instead of physical mu_e−mu_n+mu_p=0. The physical consequence is an inference from the documented table defect and visible ingestion path, without a measured production error magnitude.

- **Locations**: [NRPyEOS_stellarcollapse.c](GRHayL-main/GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse.c):42–46; [NRPyEOS_stellarcollapse_to_ghl.c](GRHayL-main/GRHayL/EOS/Tabulated/stellarcollapse/NRPyEOS_stellarcollapse_to_ghl.c):72–87; [upstream beta helper](GRHayL-main/GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c).

- **Repository evidence and checks**: The primary [StellarCollapse advisory](https://stellarcollapse.org/equationofstate.html) identifies the shifted field and corrected base-potential combination. The actual dataset mapping/copy and beta-key use support the conditional inference. No production-table runtime validation is claimed.

- **Confidence**: High for the restricted provenance/ingestion mechanism; production magnitude remains unmeasured.

- **Desired postcondition and fix direction**: Correct the table with established provenance or consistently derive the intended chemical-potential combination. Assign the remedy to table provenance/upstream EOS integration and check physical residuals.

- **Verification**: With a known affected table and a repaired counterpart, evaluate physical chemical-potential residuals in and outside the affected density region. Require the specified physical equilibrium tolerance and preserve unaffected-table behavior.

</details>

<details>
<summary><strong>Issue 17 (medium, downstream normalization):</strong> Default IllinoisGRMHD import rescales normalized Balsara fields a second time.</summary>

<br>

Original finding: `F19`.

- **Description**: GRHayLID's Balsara magnetic coefficients already use the normalized comparison convention. IllinoisGRMHD defaults rescale_magnetics to yes and multiplies imported A by another 1/sqrt(4*pi); reconstructed B is linear in A.

- **Production trigger and likelihood**: Compose a new Balsara parfile using the default Illinois magnetic conversion. Every supplied Balsara parfile explicitly disables it, so the shipped examples avoid the mismatch.

- **Impact**: For fixed fluid and metric state, the imported/reconstructed magnetic amplitude is smaller by sqrt(4*pi), and quadratic magnetic stress terms by 4*pi, relative to the intended benchmark. No measured later-time evolution error is claimed.

- **Locations**: [GRHayLID/src/1D_tests_magnetic_data.c](GRHayLID/src/1D_tests_magnetic_data.c):13–34; [IllinoisGRMHD/param.ccl](IllinoisGRMHD/param.ccl):18–20; [convert_HydroBase_to_IllinoisGRMHD.c](IllinoisGRMHD/src/convert_HydroBase_to_IllinoisGRMHD.c):7,26–28; [Balsara1.par](IllinoisGRMHD/test/Balsara1.par):76; [GRHayL derivation](GRHayL-main/docs/raw/derivation.md):79–94.

- **Repository evidence and checks**: Primary [HydroBase magnetic definitions](https://einsteintoolkit.org/thornguide/EinsteinBase/HydroBase/documentation.html), [Giacomazzo/Rezzolla equations 2.1–2.3 and Table 4](https://arxiv.org/pdf/gr-qc/0507102), coefficients, importer and all five parfiles were compared. This supports the conditional normalization finding; the supplied coefficients and examples are not alleged wrong.

- **Confidence**: High for the additional import factor under the named default setting.

- **Desired postcondition and fix direction**: Document rescale_magnetics=no for these normalized data, or enforce an explicit normalization agreement. Preserve legacy initial-data conventions when considering converter defaults.

- **Verification**: Compare intended initialized amplitudes through import and reconstruction with rescaling enabled and disabled for fixed fluid/metric data. Require the documented benchmark setting and retain compatible legacy conversions.

</details>

<details>
<summary><strong>Issue 18 (medium, external zero-density input):</strong> Standalone Hybrid entropy can publish NaN on an initialized vacuum state.</summary>

<br>

Original finding: `F20`.

- **Description**: The standalone Hybrid entropy body calls P/rho^(Gamma−1) without a local positive-density or finite-input guard. External-ID compatibility is advertised without stating that domain prerequisite.

- **Production trigger and likelihood**: An external producer deliberately initializes rho=P=0 with Gamma=2 before entropy initialization. This is a zero-density input boundary, not a positive-density scientific failure. HydroBase's zero selector alone does not prove that pressure was initialized to zero.

- **Impact**: The helper evaluates 0/0 and publishes NaN; positive P at zero rho instead gives infinity. Allocated storage and initialized input ownership do not ensure finite entropy. Native atmosphere limiting/recomputation can mask the result later.

- **Locations**: [GRHayLID/schedule.ccl](GRHayLID/schedule.ccl):51–58; [GRHayLID/src/ComputeEntropy.c](GRHayLID/src/ComputeEntropy.c):7–12; [NRPyEOS_compute_entropy_function.c](GRHayL-main/GRHayL/EOS/Hybrid/NRPyEOS_compute_entropy_function.c):27–30; [GRHayLID/README](GRHayLID/README):26–30.

- **Repository evidence and checks**: An actual entropy-body/helper probe wrote NaN for initialized rho=P=0 and Gamma=2. The formula and absent local check independently establish the mechanism. This does not imply finite physical specific entropy exists in vacuum or that every downstream consumer preserves the NaN.

- **Confidence**: High for the local zero-density/domain-guard mechanism.

- **Desired postcondition and fix direction**: Enforce positive finite density and admissible finite pressure, or define a deliberate consistent atmosphere/vacuum policy before entropy evaluation. Reject unsupported external states instead of silently publishing nonfinite output.

- **Verification**: Check positive valid inputs, explicitly initialized zero-density cases and nonfinite inputs before downstream repair. Require the documented rejection or consistent atmosphere policy without accidental NaN/infinity.

</details>

<details>
<summary><strong>Issue 19 (low, interface):</strong> The advertised magnetic staggering switch has no implementation effect.</summary>

<br>

Original finding: `F13`.

- **Description**: stagger_A_fields promises selectable +1/2 placement, but no implementation reads it. Magnetic initialization uses the unshifted base coordinates for both settings.

- **Production trigger and likelihood**: Toggle the declared switch in an otherwise valid 1D magnetic setup. The lack of parameter use applies to either setting.

- **Impact**: The promised placement cannot be selected. For these special axis-aligned potentials, omitted offsets can be curl-free gauge changes; necessarily wrong or divergent reconstructed B is not established.

- **Locations**: [GRHayLID/param.ccl](GRHayLID/param.ccl):46–48; [GRHayLID/src/1D_tests_magnetic_data.c](GRHayLID/src/1D_tests_magnetic_data.c):90–92.

- **Repository evidence and checks**: A complete source search found no parameter read. An actual-body toggle probe produced identical A values. Continuum curl and cyclic rotations are correct under the intended flat Cartesian conventions.

- **Confidence**: High for the ineffective documented option; no automatic magnetic-divergence defect is inferred.

- **Desired postcondition and fix direction**: Implement the documented component-specific placement, or remove/redefine the unused option. Preserve and document the downstream gauge/placement contract.

- **Verification**: Compare output coordinates and A values against the chosen component placement for both settings, then independently check the intended curl and gauge behavior before alleging a B correction.

</details>

<details>
<summary><strong>Issue 20 (low, documentation):</strong> Guides and interface help misstate controls, outputs, capabilities and entropy side effects.</summary>

<br>

Original finding: `F18`.

- **Description**: The guide names nonexistent compute_entropy instead of HydroBase::initial_entropy=GRHayLID, claims only A rather than both A/B output, and advertises an absent cylindrical explosion. Parameter help uses test_1D_initial_data instead of initial_data_1D. Schedule descriptions call gas/sphere 1D and say beta sets entropy. Entropy documentation omits Tabulated clamping and P/epsilon replacement.

- **Production trigger and likelihood**: A user follows the affected guide, generated help or schedule descriptions when configuring or interpreting initialization. No independent additional numerical-state failure is attributed to this documentation group.

- **Impact**: Misleading configuration instructions, capability expectations and ownership/side-effect descriptions. Related behavioral defects remain separate above rather than counted again.

- **Locations**: [GRHayLID/doc/documentation.tex](GRHayLID/doc/documentation.tex):25–32,44–48,59–64; [GRHayLID/param.ccl](GRHayLID/param.ccl):9,51–61; [GRHayLID/schedule.ccl](GRHayLID/schedule.ccl):28,38,48; [GRHayLID/README](GRHayLID/README):26–30; [GRHayLID/src/ComputeEntropy.c](GRHayLID/src/ComputeEntropy.c):28–31.

- **Repository evidence and checks**: Each description was compared with declared parameters, schedule and implementations. Line 9 is incorrect help; line 51 correctly declares the actual selector and is not a second typo. The false-positive review retains all these facts but lowers the standalone group from P2/P3 to P3.

- **Confidence**: High for the source-supported documentation discrepancies.

- **Desired postcondition and fix direction**: Use actual controls and setup names, describe A/B ownership accurately, remove unsupported capability claims, and document beta/entropy separation and primitive mutations. Keep claims consistent with any behavioral fixes.

- **Verification**: Compare the revised guide, generated parameter help, schedule descriptions and README against the implemented selector/producer matrix and all modified fields. Confirm documentation does not duplicate separately tracked behavior findings.

</details>

### Remaining ownership and validation gaps

These four items remain evidence gaps rather than additional ranked runtime issues.

- **G01 — HydroBase W ownership and ordering remain consumer-specific.** GRHayLID writes velocity but omits w_lorentz; native GRHayL consumers calculate their own u0. Ordinary initial exporters are guarded by Convert_to_HydroBase_every, whose defaults are zero. IllinoisGRMHD's legacy branch separately declares an initial W exporter after recovery regardless of that schedule guard at [schedule.ccl](IllinoisGRMHD/schedule.ccl):590–599, with the W assignment at [converter](IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c):91–94. This corrects the former blanket guard claim. A declared later exporter does not prove initialized W before an earlier consumer or actual publication in every configuration. Identify and order a metric-consistent producer for a consumer that actually reads HydroBase W; no unavoidable native failure or failing generic early consumer was demonstrated.

- **G02 — Legacy replacement ownership is visible; runtime compatibility remains unvalidated.** When ID_converter_ILGRMHD is active, [GRHayLib](GRHayL-main/implementations/GRHayLib/schedule.ccl):4–16 suppresses its normal initializer but retains termination. With IllinoisGRMHD active, [its legacy schedule](IllinoisGRMHD/schedule.ccl):565–573 selects IllinoisGRMHD_backward_compatible_initialize at CCTK_WRAGH. The [body](IllinoisGRMHD/src/backward_compatible_initialize.c):16–17,26–31,42–50 allocates the shared globals and initializes Hybrid parameters/functions; it rejects neos>1. This disproves missing replacement ownership for the visible path and removes the requirement that the absent converter itself allocate the globals. Actual legacy execution, admissible settings and initialization success remain unvalidated; other configurations still need their own owner. No runtime null dereference or Tabulated replacement was established.

- **G03 — Available validation does not establish the advertised mode matrix.** GRHayLID has no local automated test/parfile/oracle suite; 1D_tests routines generate data. Coupled [Illinois Balsara parfiles](IllinoisGRMHD/test/Balsara1.par) cover narrower Simple/x/both-magnetic/matching-normalization cases. [test.ccl](IllinoisGRMHD/test/test.ccl):23–28 disables Balsara4 with a PPM-sensitivity explanation; Balsara1/2/3/5 are declared active. Balsara4 explicitly sets max_Lorentz_factor=25, so unavoidable clipping by the generic default is excluded. These declarations prove neither a current pass nor an intrinsically incorrect Balsara4 setup. Actual ET family/selector/storage execution, standalone operations, precision/table compatibility, W ownership and the per-issue checks remain needed.

- **G04 — Density interpolation of beta roots lacks an established physical residual tolerance.** The [beta helper](GRHayL-main/GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c):173–181,253–257 interpolates density-node roots, while the ordinary [table interpolator](GRHayL-main/GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h):106 onward interpolates chemical potential itself. Those operations need not share a zero. At exact T=1, synthetic slices munu_0=Ye−.2 and munu_1=2*(Ye−.4), densities {.1,10} and Ye nodes {.1,.3,.5} give strict brackets without exact node roots. An actual-source probe returned Ye=.3 at rho=1, with table-interpolated residual −.05; the exact interpolant root is 1/3. Issues 9 and 10 are absent in this example. The [library test](GRHayL-main/Unit_Tests/unit_test_tabulated_eos.c):626–634 checks root averaging, not a physical residual tolerance. Establish that tolerance at final off-node rho/T; solve the interpolant zero there if exact table equilibrium is required. No production magnitude or scientifically unacceptable error was determined.
