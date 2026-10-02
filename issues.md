## Issues Found by Agentic LLM Review

## Remediation status

The source corrections below address all 19 defect findings and the dormant
maintenance cleanup (G03). [The owning thorn manual](NRPyLeakageET/README.md)
addresses G01. Original findings and proposed framework verification remain
below as historical context; they do not describe a current measured failure.

| Findings | Implemented correction |
| --- | --- |
| F01 | All three rate wrappers check status before output consumption; every error, including sanitized-output status, aborts with cell/state diagnostics. |
| F02, F03 | POLR uses y/z neighbors for y/z metrics and minus/plus neighbor argument order. |
| F04, F05, F06, F07 | Correct energy lapse factor, normalized metric-norm Lorentz cap, checks before invalid square roots and cell writes, allocated vector indexing. |
| F08 | Both hydro hosts refresh required HydroBase state with leakage active independently of diagnostic cadence; all converter access declarations and analysis ordering are reconciled. |
| F15, F16 | Luminosity cadence and optical-depth method are startup-only; cadence is guarded before modulo. |
| F17 | Luminosity output constructs a dynamically sized pathname. |
| F09 | Checked constrained-group registration moves to recovery-capable `MoL_Register`. |
| F10, F19 | Physical writes exclude SymBase faces; Driver/Boundary scalar symmetry selection and initialization/evolution application are added; post-boundary sync names both groups. |
| F11, F12 | POLR rejects truncated finest-level windows and reversed resolved bounds. |
| F13, F14 | Ownership uses Carpet's current local component; metadata includes helper closure, an early inaccessible-past-timelevel rejection, explicit hierarchy/timelevel validity publication and final ghost refill, plus dynamic RHS validity/modification calls. |
| F18 | Only compatibility key 0 is accepted; unsupported constants key 1 is rejected. |
| G01, G03 | Owning usage manual, corrected internal declaration/comment, dimension-aware dormant debug indexing, and structured overflow-safe iteration exit. |

G02 remains an integration-evidence gap: disposable compiled wrapper probes
provide local regression evidence, but this checkout has no configured Einstein
Toolkit executable, real EOS table, or target test/parfile suite. Full evolution,
MPI/AMR, symmetry-engine, PreSync enforcement, steering-provider, and checkpoint
results are still required. No current full-framework validation is claimed.
No maintained test suite, oracle regeneration, or CI run is supplied by this repair.

This meta-ticket collects 19 NRPyLeakageET defect findings: ten high-severity, eight medium-severity, and one low-severity finding. Three additional documentation, validation, and dormant maintenance gaps are listed separately. Three independent review seats examined the thorn, the user-designated current GRHayL implementation, direct hydro integrations, relevant framework source, and documentation; three fresh seats subsequently challenged the report for false positives. No complete defect finding was disproved. The F15/F16 steering claims below incorporate the required narrower scope from the prior false-positive analysis (not present in this checkout).

Issues are ordered by severity, with original finding IDs retained for traceability. Severity describes impact under the stated trigger, not occurrence frequency. Invalid-state, configuration-dependent, and unexecuted framework consequences are distinguished from ordinary production failures. Verification bullets describe proposed follow-up checks, not completed executions or authorization to add maintained tests.

<details>
<summary><strong>Issue 1 (high):</strong> GRHayL errors are ignored before untouched outputs are consumed. [F01]</summary>

<br>

- **Description**: The opacity, luminosity, and combined source wrappers discard the returned `ghl_error_codes_t` and consume uninitialized automatic outputs. Current kernels can return on EOS, blocking, or Fermi-integral errors before assigning those outputs.

- **Production trigger and likelihood**: An admitted cell whose rate evaluation fails. Density and psi6 thresholds do not establish valid EOS bounds/composition or successful evaluation; failure frequency was not measured.

- **Impact**: Indeterminate reads can contaminate opacities, luminosities, and evolution RHSs with arbitrary values. Ignoring the separate sanitized-output error also conceals a scientific failure.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_compute_neutrino_opacities.c:29,78,83–88`; `NRPyLeakageET/src/NRPyLeakageET_compute_neutrino_luminosities.c:63–72`; `NRPyLeakageET/src/NRPyLeakageET_compute_neutrino_opacities_and_add_source_terms_to_MHD_rhss.c:179–223`; `GRHayL-main/GRHayL/include/ghl_nrpyleakage.h:67–149`.

- **Repository evidence and checks**: Header/implementation tracing establishes the early-return contract. The original HDF5-enabled probe linked the three actual kernels with a controlled failing EOS dispatch: each returned status 10 and preserved sentinel outputs. This was error propagation, not a real-table run. Fresh false-positive review confirmed the source contract without rerunning that receipt. Upstream unchanged-output assertions were inspected, not executed.

- **Confidence**: High for caller misuse and the admitted-failure mechanism; production prevalence is unmeasured.

- **Desired postcondition and fix direction**: Check every status before reading outputs and apply an explicit abort or documented recovery policy, including sanitized-output errors. Initializing locals alone does not satisfy the contract. Preserve the upstream intentional early-return behavior.

- **Verification**: Exercise successful and early-failing calls for all three wrappers. Require no output consumption after failure, useful cell/state diagnostics, and explicit handling of sanitized-output errors.

</details>

<details>
<summary><strong>Issue 2 (high):</strong> POLR y/z metric stencils sample x neighbors. [F02]</summary>

<br>

- **Description**: `stencil_gyy` and `stencil_gzz` use `im1_j_k`/`ip1_j_k`, although the current kernel requires the corresponding y/z neighbors to compute proper face lengths.

- **Production trigger and likelihood**: A positive nonuniform metric with an affected y/z edge supplying the least-cost path. Uniform metrics or other winning edges can hide the defect; modest variation is sufficient.

- **Impact**: Incorrect optical-depth costs bias leakage suppression and derived rates.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_optical_depths_PathOfLeastResistance.c:94–96`; `GRHayL-main/GRHayL/include/ghl_nrpyleakage.h:151–177`; `GRHayL-main/GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c:33–49`.

- **Repository evidence and checks**: Fresh review compiled the unchanged wrapper with the actual current kernel. With unit opacity/spacing, true y-face metric 1, x-sampled gyy 100, and competing neighbor depths 100, the wrapper returned 7.1063352017759476 instead of 1; compile/run exited 0. An independent kernel caller with wrong-axis metric 1.1 returned 1.0246950765959599 instead of 1. These are local face-cost probes, not ET evolution.

- **Confidence**: High for the stencil mismatch and isolated numerical consequence.

- **Desired postcondition and fix direction**: Sample gyy at `i_jm1_k`/`i_jp1_k` and gzz at `i_j_km1`/`i_j_kp1`, following the current API.

- **Verification**: Use independent y-only and z-only nonuniform-metric face-cost oracles with other paths explicitly nonwinning. Include modest variations; a constant-metric case alone cannot expose the defect.

</details>

<details>
<summary><strong>Issue 3 (high):</strong> POLR reverses the current minus/plus neighbor argument order. [F03]</summary>

<br>

- **Description**: The wrapper passes plus then minus for all six opacity/depth neighbor pairs; the current API requires minus then plus. Each neighbor's opacity and depth stay paired, but the pair is attached to the opposite metric face.

- **Production trigger and likelihood**: Unequal opposite proper face lengths with an affected edge winning the minimum. Equal lengths conceal the reversal.

- **Impact**: Wrong face associations produce incorrect optical depths and subsequent leakage suppression.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_optical_depths_PathOfLeastResistance.c:124–131`; `GRHayL-main/GRHayL/include/ghl_nrpyleakage.h:164–175`; `GRHayL-main/GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_optical_depths_PathOfLeastResistance.c:18–29`.

- **Repository evidence and checks**: Fresh unchanged-wrapper/current-kernel execution with x-minus/center/x-plus metrics 1/1/9, minus-neighbor depth 0, and all other depths 100 returned 2.2360679774997898 instead of 1; compile/run exited 0. Constant y/z metrics isolate this from F02. A separate mild-variation kernel caller returned reversed cost 1.014889156509222 versus ordered cost 1.004987562112089.

- **Confidence**: High for the signature mismatch and independently isolated face error.

- **Desired postcondition and fix direction**: Pass opacity and optical-depth neighbors in the exact current minus/plus signature order. The library follows its declared contract.

- **Verification**: Isolate winning minus and plus edges on each axis with asymmetric lengths and independent face-cost expectations. Retain an equal-length case as a masking control.

</details>

<details>
<summary><strong>Issue 4 (high):</strong> The conservative energy source is missing a lapse factor. [F04]</summary>

<br>

- **Description**: For the advertised host energy `tau_tilde = alpha^2 sqrt(gamma) T^{00} - rho_star`, isotropic cooling contributes `alpha^2 sqrt(gamma) Q u^0`. The wrapper instead adds `alpha sqrt(gamma) Q u^0`.

- **Production trigger and likelihood**: Nonunit lapse with nonzero cooling in the advertised IllinoisGRMHD/GRHayLHD conservative-energy equation. Unit lapse hides the error; custom string-selected RHS conventions require separate review.

- **Impact**: Cooling is over-applied by `1/alpha` when alpha is below one.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_compute_neutrino_opacities_and_add_source_terms_to_MHD_rhss.c:201–208`; `GRHayL-main/GRHayL/Con2Prim/compute_conservs.c:82–96`; `GRHayL-main/docs/raw/derivation.md:848–875`; `GRHayL-main/GRHayL/Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h:87–103`.

- **Repository evidence and checks**: Independent conservative-source derivation gives normal projection QW and coordinate-time source `alpha sqrt(gamma) QW`. Q has no lapse/velocity argument to compensate. The original unchanged-wrapper probe, using controlled successful rates at rest with alpha=0.5 and Q=-1, added -1 instead of -0.5. Fresh review independently confirmed the derivation and algebra; no new real-EOS execution is claimed.

- **Confidence**: High for the advertised conservative-variable contract.

- **Desired postcondition and fix direction**: Add the missing lapse to the energy increment, for example `tau_rhsL = alpL * sqrtmgQ * u4U[0]`. Preserve the distinct Ye and momentum projections; neither needs this additional factor. Substituting covariant u0 is not a remedy.

- **Verification**: Compare against an independently projected covariant source at several nonunit lapses, shifts, and normalized velocities. Check Ye/momentum separately and document supported custom RHS conventions.

</details>

<details>
<summary><strong>Issue 5 (high):</strong> The leakage Lorentz cap constructs an unnormalized four-velocity. [F05]</summary>

<br>

- **Description**: The cap branch scales Eulerian velocity by `W_max/W` and independently sets W to W_max. This does not produce squared speed `1-W_max^-2` and therefore fails four-velocity normalization.

- **Production trigger and likelihood**: Finite subluminal input above the leakage cap. Aligned upstream and leakage caps can hide this branch; it does not require superluminal input.

- **Impact**: Incorrect finite momentum/source contractions despite apparently enforcing the selected cap.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_compute_neutrino_opacities_and_add_source_terms_to_MHD_rhss.c:114–137`; `GRHayL-main/GRHayL/GRHayL_Core/limit_v_and_compute_u0.c:21–45,76–92`.

- **Repository evidence and checks**: The original flat-space wrapper probe with W20/cap10 constructed u0=10 and ux≈4.993746, giving norm -75.0625 rather than -1; normalized capped ux is sqrt(99)≈9.949874. Fresh independent algebra reproduced this, plus failures for W5/cap2 and cap1. The current GRHayL helper uses norm-based scaling and returns a checked error status.

- **Confidence**: High for finite above-cap input and the normalization defect.

- **Desired postcondition and fix direction**: Reuse the current checked `ghl_limit_v_and_compute_u0` with the intended cap and transport-velocity convention, or implement equivalent norm-based scaling and recomputation. Preserve error handling and do not silently substitute a different host cap.

- **Verification**: Check normalization and the requested cap with finite above-cap states, cap1, nonzero shift, and nontrivial positive spatial metrics. Confirm below-cap states remain correct.

</details>

<details>
<summary><strong>Issue 6 (high, invalid-state failure path):</strong> Newly generated nonfinite sources escape the NaN guard. [F06]</summary>

<br>

- **Description**: The wrapper takes the Lorentz-factor square root before validating metric speed. Its subsequent `isnan` product includes rates/opacities and old RHS values, but excludes W, four-velocity, and newly formed increments; it also fails to reject infinity.

- **Production trigger and likelihood**: An invalid velocity/geometry combination, such as finite Eulerian speed 1.1, reaches the wrapper. Healthy upstream recovery can prevent this trigger; an ordinary successful recovery is not alleged to supply it. The product can also produce misleading diagnostics through overflow and multiplication by zero.

- **Impact**: Nonfinite energy/momentum increments can be committed without the intended error response.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_compute_neutrino_opacities_and_add_source_terms_to_MHD_rhss.c:119–129,185–223,247`.

- **Repository evidence and checks**: The original unchanged-wrapper probe supplied speed1.1, valid geometry, finite successful controlled rates, and zero initial RHSs. It wrote NaN energy/momentum and returned without its error callback. Fresh review confirmed the control flow. This is distinct from F05's wrong finite result on subluminal input.

- **Confidence**: High for the stated invalid-state mechanism; occurrence after healthy host recovery is not established.

- **Desired postcondition and fix direction**: Validate geometry/speed before the square root, use checked velocity construction, and check individual actual outputs/increments for finiteness before committing them. If fast-math is supported, use classification that remains correct under that configuration.

- **Verification**: Supply superluminal speed, nonfinite inputs, and finite-rate/nonfinite-increment combinations; require the chosen failure policy before writes. Include a finite case that exposes misleading product-based diagnostics. No fast-math run is claimed.

</details>

<details>
<summary><strong>Issue 7 (high, padded allocation):</strong> HydroBase velocity component offsets use the local rather than allocated volume. [F07]</summary>

<br>

- **Description**: Manual velx/vely/velz bases use `prod(cctk_lsh)`, while Cactus vector components are separated by `prod(cctk_ash)`. Correct scalar point indexing cannot repair an incorrect component base.

- **Production trigger and likelihood**: Allocation padding makes allocated and local volumes differ. Equal volumes hide the defect; no particular Carpet padding frequency was measured.

- **Impact**: Wrong velocity components change the Lorentz factor and energy/momentum leakage contractions.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_compute_neutrino_opacities_and_add_source_terms_to_MHD_rhss.c:5–7,99–101`; the direct hydro exporters' `CCTK_VECTGFINDEX3D` accesses.

- **Repository evidence and checks**: The [Cactus indexing header](https://bitbucket.org/cactuscode/cactus/raw/master/src/include/cctk_core.h) defines allocated vector stride. An original framework-consistent shim with lsh=(3,3,3), ash=(4,3,3), and velocity (0.1,0.2,0.3) produced Sz≈-0.209657 instead of -0.323498. Fresh review confirmed the allocation contract; the shim was not a Carpet launch.

- **Confidence**: High for a supported padded layout.

- **Desired postcondition and fix direction**: Use `vel[CCTK_VECTGFINDEX3D(cctkGH,i,j,k,component)]` or a framework-derived allocated component base.

- **Verification**: Exercise unequal allocated/local shapes with distinct component sentinels and compare leakage contractions to an independent component-wise oracle. Retain the unpadded case.

</details>

<details>
<summary><strong>Issue 8 (high, hydro integration):</strong> The HydroBase refresh required by leakage is gated by diagnostic cadence. [F08]</summary>

<br>

- **Description**: IllinoisGRMHD and GRHayLHD schedule their converters during hydro RHS evaluation when leakage is active, but the converter bodies still gate execution by an analysis cadence whose default is zero.

- **Production trigger and likelihood**: Leakage with default-zero or larger-than-one diagnostic cadence. GRHayLHD reaches modulo zero; IllinoisGRMHD does likewise without the legacy converter thorn. Its legacy zero-cadence branch instead returns without refreshing required data. Cadence1 mitigates these paths.

- **Impact**: Crash or stale HydroBase velocity/W during leakage source evaluation. The defect belongs to direct hydro integrations, not the GRHayL leakage kernel.

- **Locations**: `IllinoisGRMHD/schedule.ccl:197–205`; `IllinoisGRMHD/src/convert_IllinoisGRMHD_to_HydroBase.c:8–16`; `IllinoisGRMHD/param.ccl:9–12`; `GRHayLHD/schedule.ccl:137–145`; `GRHayLHD/src/convert_GRHayLHD_to_HydroBase.c:7–8`; `GRHayLHD/param.ccl:9–12`.

- **Repository evidence and checks**: Direct sibling schedule/body tracing establishes the unconditional diagnostic gate on the leakage-required RHS hook. Its metadata also omits some full converter accesses, including metric reads and W writes. Neither the original audit nor fresh review executed a coupled ET run.

- **Confidence**: High for the declared hook and control-flow mismatch; coupled stage freshness remains unexecuted.

- **Desired postcondition and fix direction**: Refresh the velocity/W needed by leakage at every required RHS stage independently of optional diagnostics. Safely guard optional output and declare the converter's actual accesses. Check analysis W freshness separately for luminosities.

- **Verification**: Trace stages with diagnostic cadences 0, 1, and greater than one, including the Illinois legacy branch. Require current velocity/W before every leakage use and no modulo-zero path.

</details>

<details>
<summary><strong>Issue 9 (high, conditional runtime steering):</strong> An installed luminosity routine can divide by zero or ignore disabling changes. [F15]</summary>

<br>

- **Description**: `compute_luminosities_every` accepts every integer and is STEERABLE ALWAYS, with negatives advertised as disabling output. Once installed, the routine performs modulo without a nonpositive guard; the thorn supplies no runtime allocation/rebinding mechanism.

- **Production trigger and likelihood**: A positive startup cadence is changed to zero/negative through an update that leaves installed schedule/storage unchanged. Positive-to-zero reaches modulo zero on the next call; -1 still outputs every iteration. Starting disabled and changing only parameter data to positive cannot install the omitted path. No steering-provider execution or frequency estimate is available. The disabled startup default itself does not immediately divide by zero.

- **Impact**: Crash, ineffective disabling, or ineffective enabling under that unchanged-schedule/storage condition.

- **Locations**: `NRPyLeakageET/param.ccl:124–127`; `NRPyLeakageET/schedule.ccl:83–95`; `NRPyLeakageET/src/NRPyLeakageET_compute_neutrino_luminosities_global_sum_and_output_to_file.cc:14`.

- **Repository evidence and checks**: The original exact-expression UBSan probe reported division by zero and exited1; this was a defect demonstration, not a steering run. Fresh review traced the [binding generator](https://bitbucket.org/cactuscode/cactus/raw/master/lib/sbin/CreateScheduleBindings.pl), [initializer](https://bitbucket.org/cactuscode/cactus/raw/master/src/main/InitialiseCactus.c), and [parameter setter](https://bitbucket.org/cactuscode/cactus/raw/master/src/main/Parameters.c): native assignment has no intrinsic schedule/storage reconstruction, but registered callbacks are possible. The [Users Guide](https://einsteintoolkit.org/usersguide/UsersGuide.html) contains both interface-driven dynamic-scheduling prose and startup-only external-condition prose; it cannot prove universal provider behavior.

- **Confidence**: High for the unsafe expression and stated unchanged-schedule/storage path; arbitrary rebinding providers remain unverified.

- **Desired postcondition and fix direction**: Guard cadence<=0 before modulo and provide real runtime scheduling/storage support, or make the parameter startup-only. A guard alone cannot enable functionality omitted at startup.

- **Verification**: Name and inspect/run the actual steering provider. Check schedule/storage and output before and after positive-to-zero, positive-to-negative, and disabled-to-positive changes. Distinguish new-process recovery from live steering.

</details>

<details>
<summary><strong>Issue 10 (high, long output path):</strong> Luminosity filename construction can overflow a stack buffer. [F17]</summary>

<br>

- **Description**: Unrestricted output-directory and filename strings are joined with `sprintf` into `char filename[512]`. A joined pathname of at least 512 characters, excluding the terminator, exceeds the buffer.

- **Production trigger and likelihood**: Sufficiently long configured directory/name combinations. Valid nested paths can trigger it; invalid names also corrupt memory before `fopen` rejects them. Short ordinary paths are unaffected.

- **Impact**: Stack corruption or process termination on rank0, sufficient to disrupt the run.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_compute_neutrino_luminosities_global_sum_and_output_to_file.cc:68–71`; `NRPyLeakageET/param.ccl:129–133`.

- **Repository evidence and checks**: Direct buffer-size proof survives fresh review. The original exact-statement probe with a 600-character name and out_dir="." produced an AddressSanitizer stack-buffer-overflow for a 603-byte write and exited1. That receipt was not rerun during false-positive review and did not execute Carpet output.

- **Confidence**: High for the configured long-string mechanism.

- **Desired postcondition and fix direction**: Construct the path dynamically, or use bounded formatting and reject truncation before opening it.

- **Verification**: Exercise boundary-length, overflowing, and valid long nested paths under a memory checker. Require either complete safe construction or an explicit error without opening a truncated path.

</details>

<details>
<summary><strong>Issue 11 (medium, checkpoint recovery):</strong> Recovery omits leakage's constrained-variable MoL registration. [F09]</summary>

<br>

- **Description**: Both constrained-group registrations occur only inside the optical-depth initial-data routine at CCTK_INITIAL. The inspected Carpet recovery branch skips initial-data calls while MoL constructs fresh registration bookkeeping.

- **Production trigger and likelihood**: Checkpoint recovery in a new process using that lifecycle. Fresh-run initial data executes the current registration route; exact restart damage was not measured.

- **Impact**: Recovered leakage groups lack the registration-dependent constrained time-level/validity handling available on a fresh run. Corrupted recovered data is not demonstrated.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_utils.cc:283–285`; `NRPyLeakageET/schedule.ccl:4–20`; the framework recovery/registration lifecycle.

- **Repository evidence and checks**: The [MoL contract](https://einsteintoolkit.org/thornguide/CactusNumerical/MoL/documentation.html) and [schedule](https://bitbucket.org/cactuscode/cactusnumerical/raw/master/MoL/schedule.ccl) place registration at MoL_Register/WRAGH. [Carpet initialization](https://raw.githubusercontent.com/EinsteinToolkit/Carpet/master/Carpet/src/Initialise.cc) skips initial data on recovery; [MoL InitialCopy](https://bitbucket.org/cactuscode/cactusnumerical/raw/master/MoL/src/InitialCopy.c) uses constrained registration for old-to-current copies. Original and fresh reviews traced source; neither ran a checkpoint recovery.

- **Confidence**: High for the lifecycle omission; its precise recovered-data consequence remains unverified.

- **Desired postcondition and fix direction**: Move checked registration into a dedicated recovery-capable MoL_Register routine. Keep initial optical-depth calculation at CCTK_INITIAL without rerunning initial data during recovery.

- **Verification**: Compare fresh-start and recovered-process registration state, constrained copies, and validity/time-level behavior. Verify POLR's interior-only updates receive the intended lifecycle support.

</details>

<details>
<summary><strong>Issue 12 (medium, symmetry-dependent):</strong> Physical boundary writes can zero optical depth on symmetry faces. [F10]</summary>

<br>

- **Description**: The boundary routine zeroes tau/kappa on flagged coarsest bbox faces after iteration0 without excluding registered reflection/rotation faces. Scalar parity registration does not prevent this direct write; initialization also lacks target-owned symmetry application/Boundary selection.

- **Production trigger and likelihood**: A flagged face is a symmetry face. A false escape path requires POLR to consume the zeroed values before any external/driver reflection refill. Persistence depends on the actual boundary and PreSync schedule.

- **Impact**: Certain symmetry-incompatible writes; conditionally, an artificial transparent boundary and biased optical depths. Sustained evolution damage is not established.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_BoundaryConditions.c:25–26,39–62`; `NRPyLeakageET/src/NRPyLeakageET_InitSym.c:10–12`; `NRPyLeakageET/src/NRPyLeakageET_utils.cc:214–228`; `NRPyLeakageET/schedule.ccl:63–80`.

- **Repository evidence and checks**: The [SymBase contract](https://einsteintoolkit.org/thornguide/CactusBase/SymBase/documentation.html) requires checking symmetry handles in addition to bbox. An original unchanged-boundary probe set constant tau/kappa7 to0 at z-min where even reflection requires7. Fresh review confirmed the write and distinguished it from later refill; no reflection-engine run is claimed. [CartGrid3D's implementation](https://bitbucket.org/cactuscode/cactusbase/raw/master/CartGrid3D/src/Symmetry.c) shows the role of variable selection.

- **Confidence**: High for the overwrite; conditional for persistent numerical effects and missing initialization ghost fills.

- **Desired postcondition and fix direction**: Separate physical and symmetry faces through established SymBase/Boundary APIs. Ensure even-scalar ghost conditions are applied during each required initialization iteration and evolution stage.

- **Verification**: Run reflection/rotation cases with nonzero symmetric scalar data and inspect ghost values immediately before each POLR use, including initialization and relevant PreSync modes.

</details>

<details>
<summary><strong>Issue 13 (medium, truncated AMR solve):</strong> Global convergence masking can hide every changing solved point. [F11]</summary>

<br>

- **Description**: Squared differences are updated only on the selected initialization levels, but the GLOBAL reduction uses the full hierarchy's coverage mask. Finer active grids excluded from the solve can mask the changing coarse values.

- **Production trigger and likelihood**: A nonempty truncated solve window, for example maxInitRefLevel=1 with active level2 covering every changing emitting coarse point. The maximum0 sentinel selects all levels and avoids this excluded-finer-level trigger.

- **Impact**: The reported norm can be zero after the first sweep while the selected coarse solution is still changing, causing premature convergence.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_utils.cc:96–148,196–197,232–242`; `NRPyLeakageET/src/NRPyLeakageET_compute_optical_depth_change.c:44–95`.

- **Repository evidence and checks**: Current [CarpetReduce reduction](https://raw.githubusercontent.com/EinsteinToolkit/Carpet/master/CarpetReduce/src/reduce.cc) retains coverage weights for isum; it removes relative volume weighting, not masking. [Mask construction](https://raw.githubusercontent.com/EinsteinToolkit/Carpet/master/CarpetReduce/src/mask_carpet.cc) suppresses covered coarse points. Excluded fine differences remain initialized zero. Independent reviews retained this conditional source counterexample; no AMR execution was performed.

- **Confidence**: High for the stated source-derived masked-domain counterexample; numerical AMR reproduction remains pending.

- **Desired postcondition and fix direction**: Reduce over a mask/domain matching the selected solve, populate differences on the matching composite hierarchy, or reject unsupported truncated windows. Replacing isum with sum alone retains the masking error.

- **Verification**: Use at least three levels with all changing coarse emitters inside excluded fine coverage. Compare the reported norm with an independent norm over the solved domain; include the all-level case.

</details>

<details>
<summary><strong>Issue 14 (medium, accepted configuration):</strong> Inverted initialization bounds silently skip the optical-depth solve. [F12]</summary>

<br>

- **Description**: Minimum and maximum refinement bounds are independently accepted/clipped without checking their resolved order. Initialization then treats an empty solve interval as converged zero data.

- **Production trigger and likelihood**: At least three levels with minInitRefLevel=2 and maxInitRefLevel=1, or another accepted pair resolving to start>end. This requires a contradictory configuration, not invalid individual parameter values.

- **Impact**: Requested POLR initialization performs no solve and leaves zero optical depths, including in opaque material.

- **Locations**: `NRPyLeakageET/param.ccl:21–30`; `NRPyLeakageET/src/NRPyLeakageET_utils.cc:196–197,214–242,263–268`.

- **Repository evidence and checks**: Source tracing shows initial zeroing followed by empty POLR/difference loops and reduction of auxiliary zeros. Available-level clipping does not order reversed bounds. Independent review distinguishes this from explicitly selecting zero initialization and from F11's nonempty masked solve. Evidence is static.

- **Confidence**: High for the accepted reversed-bound control flow.

- **Desired postcondition and fix direction**: Validate resolved start<=end before solving and reject inconsistent inputs. Document maximum0 as the finest-level sentinel and explain clipping.

- **Verification**: Check reversed, equal, clipped, and sentinel bounds against their resolved intervals. Require an explicit error for an empty requested POLR interval and actual execution for valid intervals.

</details>

<details>
<summary><strong>Issue 15 (medium, multiple maps or multigrid):</strong> Carpet ownership is indexed by multigrid level instead of map. [F13]</summary>

<br>

- **Description**: The ownership helper queries `vhh.AT(Carpet::mglevel)->processor(reflevel,component)` inside map loops, but vhh is indexed by map. With mglevel0, every map consults map0.

- **Production trigger and likelihood**: Maps have different component ownership, or mglevel exceeds the available map index. Ordinary single-map/mglevel0 configurations hide the error; no multimap/MPI frequency estimate is available.

- **Impact**: Locally owned data can be skipped, nonlocal data claimed, or an out-of-range map accessed. Initialization, opacity, convergence, and luminosity helpers share this check.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET_utils.cc:20–24`, and its initialization/opacity/convergence/luminosity callers.

- **Repository evidence and checks**: Official [Carpet declarations](https://raw.githubusercontent.com/EinsteinToolkit/Carpet/master/Carpet/src/variables.hh) distinguish map, mglevel, and local_component and declare vhh by map. Target traversal and query were compared directly. No multimap runtime was executed. MPI_COMM_WORLD deserves a separate communicator check, but an ordinary-run rank mismatch is not established.

- **Confidence**: High for the wrong index; exact distributed consequences remain conditional and unexecuted.

- **Desired postcondition and fix direction**: Use driver-local ownership state or the current map's owner query with the appropriate driver communicator.

- **Verification**: Assign different owners on multiple maps, vary multigrid level, and check each component exactly against driver ownership. Inspect communicator-split behavior without assuming world ranks equal driver ranks.

</details>

<details>
<summary><strong>Issue 16 (medium, driver validity contract):</strong> Scheduling metadata omits helper accesses and dynamic RHS updates. [F14]</summary>

<br>

- **Description**: Initialization omits off-diagonal metric reads, opacity time-level writes, and auxiliary writes from its actual helper access closure. Five runtime string-selected hydro RHS arrays are read/updated without declarations or dynamic validity/modification notifications.

- **Production trigger and likelihood**: A driver/PreSync configuration enforcing those access contracts. Surrounding thorn metadata and stage ordering determine the first consequence; universal failure with enforcement disabled is not claimed.

- **Impact**: Missing synchronization, stale/poisoned reads, or inconsistent validity bookkeeping under enforcement.

- **Locations**: `NRPyLeakageET/schedule.ccl:8–19,31–39`; `NRPyLeakageET/src/NRPyLeakageET_optical_depths_initialization_routines.c:16–42,78–105`; `NRPyLeakageET/src/NRPyLeakageET_compute_neutrino_opacities.c:42–78,90–96`; `NRPyLeakageET/src/NRPyLeakageET_compute_optical_depth_change.c:59–95`; `NRPyLeakageET/src/NRPyLeakageET_compute_neutrino_opacities_and_add_source_terms_to_MHD_rhss.c:22–26,185–223`.

- **Repository evidence and checks**: Complete helper/source access tracing establishes the metadata mismatch. The [PreSync contract](https://docs.einsteintoolkit.org/et-docs/PreSync) supports static metadata and Driver_RequireValidData/Driver_NotifyDataModified for runtime-selected arrays. Fresh review checked [Carpet's public wrappers](https://raw.githubusercontent.com/EinsteinToolkit/Carpet/master/Carpet/src/PreSync.cc), including harmless disabled-mode returns. Unsuffixed READS and repeated _p names are legal. No enforced ET traversal was run.

- **Confidence**: High for incomplete declarations/notifications; runtime effects depend on the actual driver contract.

- **Desired postcondition and fix direction**: Declare actual regions, all relevant time levels, and the full helper access closure. Use supported driver APIs with correct arguments around dynamic RHS reads/writes.

- **Verification**: Trace initialization and source stages with validity enforcement/poisoning enabled and check synchronization and modification state. Confirm disabled-mode compatibility and every dynamically selected variable's time level/region.

</details>

<details>
<summary><strong>Issue 17 (medium, conditional runtime steering):</strong> Updating the method selector alone does not switch the installed optical-depth routine. [F16]</summary>

<br>

- **Description**: `optical_depth_evolution_type` is STEERABLE ALWAYS, but an outer schedule conditional selects static or POLR at initialization. The thorn supplies no runtime dispatcher or rebinding callback.

- **Production trigger and likelihood**: Live selector changes through a parameter update that does not reconstruct the schedule. This leaves the previously installed method active. Behavior of a particular external rebinding provider was not verified. Recovery in a new process is different: recovered parameter values are loaded before schedule initialization and can select a new branch then.

- **Impact**: Configured method and executed method diverge under the unchanged-schedule condition.

- **Locations**: `NRPyLeakageET/param.ccl:71–75`; `NRPyLeakageET/schedule.ccl:41–61`; the installed static/POLR step routines.

- **Repository evidence and checks**: Source inspection finds no current-selector dispatch in either installed routine or implementing notification registration. The binding/initializer/native-setter trace and conflicting guide passages described in F15 support only the narrowed path, not every possible steering interface. Fresh false-positive review required and accepted this qualification; no provider execution was performed.

- **Confidence**: High for absent thorn-owned dispatch and the stated unchanged-schedule path; arbitrary provider behavior remains unverified.

- **Desired postcondition and fix direction**: Schedule a runtime dispatcher, implement verified runtime rebinding support, or restrict the selector to startup use. Document when initialization-only selectors take effect.

- **Verification**: With a named provider, change static to POLR and back while observing actual invoked routines. Separately recover into a new process with changed selector values and verify branch selection before evolution.

</details>

<details>
<summary><strong>Issue 18 (medium, scientific configuration):</strong> The advertised physical-constants selector has no consumer. [F18]</summary>

<br>

- **Description**: `constants_key` advertises NRPy versus HARM3D+NUC constants and runtime steering, but no target code reads it. Current rate signatures have no selection argument and kernels use fixed macros.

- **Production trigger and likelihood**: A user requests the alternative key1 convention, including through steering. Key0 and key1 reach the same current implementation.

- **Impact**: Scientific configuration and comparison labels can claim a constants model that was never selected. Neither historical set is alleged to be intrinsically invalid.

- **Locations**: `NRPyLeakageET/param.ccl:14–18`; the three rate call sites in F01; `GRHayL-main/GRHayL/include/ghl_nrpyleakage.h:9–55` and current signatures.

- **Repository evidence and checks**: Complete target parameter-use search finds only the declaration; direct API/kernel tracing shows fixed constants. Alternative unused ZL macros do not implement selection. Both audit and false-positive review established this statically; no model-comparison execution is claimed.

- **Confidence**: High for the unconsumed selector and unchanged model path.

- **Desired postcondition and fix direction**: Remove/deprecate or explicitly reject the unsupported choice and document the actual convention, or expose it only through a supported implemented library mechanism. Do not invent an argument to the current API.

- **Verification**: Require each accepted choice to reach its documented constants model and produce independently expected rate differences, or require explicit rejection of unsupported values. Check existing parameter-file compatibility.

</details>

<details>
<summary><strong>Issue 19 (low, declaration defect):</strong> Post-boundary SYNC repeats optical depths and omits opacities. [F19]</summary>

<br>

- **Description**: The boundary routine writes depths and opacities, but its final SYNC lists `NRPyLeakageET_optical_depths` twice and omits the opacity group.

- **Production trigger and likelihood**: The declaration is present in the ordinary boundary schedule. Whether missing opacity synchronization affects a particular decomposition depends on ownership and later synchronization.

- **Impact**: Definite declaration drift; a separate MPI opacity failure is not demonstrated.

- **Locations**: `NRPyLeakageET/schedule.ccl:80`; `NRPyLeakageET/src/NRPyLeakageET_BoundaryConditions.c`.

- **Repository evidence and checks**: Literal declaration and boundary writes were compared by the original and fresh reviewers. No distributed boundary execution established whether additional opacity communication is necessary or already supplied later.

- **Confidence**: High for the duplicate/omission, without a measured MPI consequence.

- **Desired postcondition and fix direction**: Name the intended opacity group if post-boundary synchronization is required, or deliberately remove/document unnecessary synchronization after checking ownership and the surrounding schedule.

- **Verification**: Trace boundary ownership and subsequent communication on a decomposed grid. Confirm the final declarations match the intended contract for both groups.

</details>

### Additional documentation, validation, and maintenance gaps

These retain their original gap/dormant classifications. They are not three additional demonstrated evolution failures.

<details>
<summary><strong>Gap G01 (documentation):</strong> The thorn has no owning usage or integration manual.</summary>

<br>

- **Description**: The supplied thorn has no README or `doc/documentation.tex`. Parameters and scattered comments do not establish a complete EOS, units, source-projection, driver, or lifecycle contract.

- **Production trigger and likelihood**: Users configure or integrate the thorn from its local interface. The workspace README directs them to individual thorn manuals, but none is supplied here.

- **Impact**: Missing guidance for scientifically meaningful and supported configurations; absence of prose alone does not prove a new numerical error.

- **Locations**: `NRPyLeakageET/`; `NRPyLeakageET/param.ccl`; `README.md:25–31`.

- **Repository evidence and checks**: Target inventory and workspace guidance were inspected. The psi6 cutoff is coordinate-dependent and cannot generally locate a horizon. Misleading method/constants promises are separately captured by F16/F18.

- **Confidence**: High for the documentation gap in the supplied tree.

- **Desired postcondition and fix direction**: Document compatible tabulated EOS/units, conservative RHS and Q/R/tau conventions, initialization sentinels/clipping, supported drivers, symmetry/AMR/recovery assumptions, errors, output semantics, and steering limits. Describe psi6 as a configured heuristic.

- **Verification**: Reconcile each documented contract with current source and show a coherent supported configuration without implying unperformed runtime validation.

</details>

<details>
<summary><strong>Gap G02 (validation):</strong> Checked-in library evidence does not validate the ET wrapper.</summary>

<br>

- **Description**: The target has no checked-in tests or target parfiles. Upstream kernel tests do not establish this wrapper's stencil packing, projections, allocated indexing, boundaries, scheduling, or recovery.

- **Production trigger and likelihood**: A wrapper/integration regression occurs outside existing library coverage. No regression rate or new defect follows merely from absent coverage.

- **Impact**: Reduced ability to establish and retain wrapper/framework correctness; this is an evidence gap.

- **Locations**: `NRPyLeakageET/`; `GRHayL-main/Unit_Tests/unit_test_nrpyleakage_physics.c`; wrapper call sites and schedules identified above.

- **Repository evidence and checks**: Original and fresh reviewers inspected the target inventory and distinguished library oracles from wrapper checks. Disposable probes expose local errors; neither those probes nor historical original-thorn parfiles prove current ET execution. No complete upstream suite was run.

- **Confidence**: High for the stated coverage gap, without an additional numerical-defect claim.

- **Desired postcondition and fix direction**: Establish authorized wrapper and ET validation with independent face-cost/source-projection oracles, successful/failing statuses, padded indexing, and lifecycle/boundary/AMR coverage. Maintained test changes require separate explicit authorization.

- **Verification**: Demonstrate actual evolution, recovery, and named-provider steering in supported configurations; distinguish kernel, shim, and full-framework results. No maintained tests, assertions, fixtures, generators, or CI runs were added here.

</details>

<details>
<summary><strong>Gap G03 (dormant maintenance and extreme parameters):</strong> Internal declarations, debug indexing, and the convergence exit need bounded cleanup.</summary>

<br>

- **Description**: An internal header declares an old RHS wrapper name; an uncalled print helper indexes diagonal positions 0..11 regardless of grid size; a boundary comment refers to A_mu. Convergence exits by assigning `i=max_iterations*10`, which can overflow a 32-bit signed expression for legal very large settings.

- **Production trigger and likelihood**: Declaration/comment drift is dormant. The print helper has no target caller or schedule. The overflow requires both a 32-bit integer expression and sufficiently large legal max_iterations reaching convergence; ordinary defaults are not alleged to fail.

- **Impact**: Maintenance ambiguity and latent helper bounds risk; an extreme-configuration arithmetic hazard. No ordinary link failure or reachable print crash is established.

- **Locations**: `NRPyLeakageET/src/NRPyLeakageET.h:26`; `NRPyLeakageET/src/NRPyLeakageET_utils.cc:163–175,241`; `NRPyLeakageET/src/NRPyLeakageET_BoundaryConditions.c`.

- **Repository evidence and checks**: Target symbol/caller search and direct indexing/comment inspection establish the dormant scope. The exit multiplication admits signed overflow under the stated integer/range condition. No runtime reproduction of this extreme case or print call is claimed.

- **Confidence**: High for the inspected drift and bounded arithmetic mechanism; dormant code is not presented as a demonstrated production failure.

- **Desired postcondition and fix direction**: Align/remove the stale declaration, correct the comment, and use dimension-aware indexing if the print helper is retained or made callable. Exit the convergence loop with structured control flow rather than overflow-prone multiplication.

- **Verification**: Compare declarations, definitions, and callers. Check small local dimensions if the helper is used, and inspect/run the convergence exit with large supported settings under integer-overflow detection.

</details>
