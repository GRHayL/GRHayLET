# NRPyLeakageET

This Carpet thorn computes neutrino opacities, optical depths, conservative
leakage sources, and optional integrated luminosities using GRHayL's NRPyLeakage
API. It requires Carpet, MPI, HDF5, GRHayLib, ADMBase, HydroBase, Boundary,
CartGrid3D scalar parity support, SymBase, SpaceMask, and MoL. CarpetX is not implemented:
the global initialization and analysis helpers use Carpet's hierarchy API.

## EOS, units, and sources

Use a compatible HDF5-backed tabulated GRHayL EOS with density, electron fraction,
and temperature within its supported domain. HydroBase `rho` uses GRHayL
geometric density units; `temperature` and EOS chemical potentials use MeV.
Opacity moments use inverse GRHayL length units and optical depths are
dimensionless. The GRHayL NRPyLeakage header owns the fixed constants and unit
conversions. `constants_key=0` is retained for parameter-file compatibility;
the previously advertised choice `1` is rejected because no alternative model
is implemented by this API.

For Eulerian velocity `V^i`, the wrapper constructs transport velocity
`v^i=alpha V^i-beta^i` and a normalized four-velocity. Finite subluminal states
above `W_max` are rescaled using their spatial metric norm. Superluminal,
nonfinite, or invalid geometry inputs abort. `W_max=1` is supported.

The RHS names must identify timelevel-zero arrays with the host definitions
`Ye_star=sqrt(gamma) rho W Ye`,
`tau_tilde=alpha^2 sqrt(gamma) T^{00}-rho_star`, and
`S_tilde_i=alpha sqrt(gamma) T^0_i`. Given GRHayL matter sources `R` and `Q`,
the increments are

```text
Ye_star_rhs += alpha sqrt(gamma) R
 tau_rhs   += alpha^2 sqrt(gamma) Q u^0
 S_i_rhs   += alpha sqrt(gamma) Q u_i
```

The array-name parameters select storage, not a different projection or units.
Custom hosts must supply these same conventions and current HydroBase state.
The IllinoisGRMHD and GRHayLHD hooks refresh HydroBase velocity and Lorentz
factor before leakage RHS and luminosity analysis independently of optional
hydro diagnostic cadence. Dynamic RHS reads and modifications are reported to
Carpet through its PreSync APIs at timelevel zero in the interior.

Every non-success rate status aborts with cell/level and thermodynamic state,
including statuses that return sanitized finite outputs. Such fallback outputs
are not silently adopted. Successful calls must also return finite outputs;
source increments and proposed accumulated RHS values are checked individually
before that cell is written. This policy does not roll back other cells before
termination. It requires the current status-returning GRHayL API and its
`robust_isfinite` classifier; on binary64 this preserves classification under
fast-math. Other floating representations follow the classifier's documented
portable limitations.

## Initialization, boundaries, and recovery

Initialization writes all three stored opacity/depth timelevels. It rejects a
Carpet initialization context that hides past arrays before any helper writes.
Use `InitBase::initial_data_setup_method="init_some_levels"` with
`Carpet::init_each_timelevel=no` to keep past arrays accessible. The
`init_single_level` path with past access disabled is unsupported.

`initial_optical_depth` chooses zero initialization or iterative POLR at startup.
The POLR stencil uses axis-matched metric neighbors and minus/plus opacity/depth
pairs. `maxInitRefLevel=0` means the finest active level. Bounds are clipped to
available levels and reversed resolved intervals are rejected. POLR must include
the finest active level: truncated windows are rejected because the composite
Carpet reduction masks covered coarse points. A nonzero minimum can omit coarse
solve levels; their material opacities are computed with zero initial depths and
copied to all donor timelevels before the first sweep, so opacity prolongation
has computed coarse data. They do not run POLR sweeps; the optical-depth result
is restricted to coarser levels after convergence.

The stopping value is the square root of the composite sum of six squared
relative optical-depth changes on admitted material cells. `tauChangeThreshold`
is positive. A maximum-iteration exit warns when the threshold was not met;
it does not certify convergence. Refinement interpolation/restriction and
composite masks remain Carpet operations, not independent numerical oracles.

Scalar symmetry handling is selected through Driver/Boundary APIs. Initialization
uses Carpet's component-aware `ApplyBCs` traversal during the iterative solve;
evolution applies it after physical boundary writes. Zero initialization is
published across all levels and timelevels before solving. After final restriction,
current ghosts and symmetries are refilled coarse to fine before copying and
publishing past timelevels. Registered symmetry faces are excluded from the zero physical
outer boundary on the coarsest level. Physical boundary zeroing is skipped at
iteration zero. Evolution and initial-data ghost values must be checked with the
actual reflection/rotation, AMR, and PreSync configuration before scientific use.

Constrained opacity/depth groups are registered in `MoL_Register`, including a
new process recovering a checkpoint. Initial optical-depth computation remains
in `CCTK_INITIAL`; recovery must restore the stored timelevels rather than rerun
initial data. The change does not supply a checkpoint/runtime oracle.

Density thresholds select cells admitted for leakage. `psi6_threshold` is a
coordinate-dependent suppression heuristic, not a horizon finder. Neither cutoff
establishes EOS validity or a successful rate evaluation.

## Configuration and output

For an existing tabulated-EOS IllinoisGRMHD configuration, add NRPyLeakageET and
its required thorns to `ActiveThorns`, then select, for example:

```text
InitBase::initial_data_setup_method = "init_some_levels"
Carpet::init_each_timelevel = no
NRPyLeakageET::initial_optical_depth = "PathOfLeastResistance"
NRPyLeakageET::optical_depth_evolution_type = "PathOfLeastResistance"
NRPyLeakageET::constants_key = 0
NRPyLeakageET::minInitRefLevel = 0
NRPyLeakageET::maxInitRefLevel = 0
NRPyLeakageET::compute_luminosities_every = 10
```

This fragment augments host/EOS/driver setup; it is not a complete tested parfile.
GRHayLHD hosts additionally set all five `GFstring_*_rhs` parameters to the
corresponding `GRHayLHD::Ye_star_rhs`, `tau_rhs`, and `Stildex/y/z_rhs` names.

Initialization controls, constants compatibility key, evolution method, and
luminosity cadence are startup-only. Changing these values in a running process
is unsupported. New-process recovery selects its schedule/storage from recovered
startup parameters. Other parameters marked `STEERABLE ALWAYS` retain their
individual runtime semantics; changing array names requires compatible host data.

A positive `compute_luminosities_every` allocates and schedules luminosities;
zero or negative disables them. Output appends to `IO::out_dir/luminosities_outfile`
using a dynamically sized path. Columns are iteration, coordinate time, electron
neutrino luminosity, electron antineutrino luminosity, and one heavy-lepton species.
Multiply the last value by four for that species family and geometric luminosity
by `c^5/G` for erg/s. Global integration uses Carpet's composite reduction and
coordinate cell volume. Recovery appends to the existing file; no duplicate-time
suppression or automatic file rotation is supplied.

## Validation limits and required framework checks

Disposable wrapper probes can verify ordered face costs, conservative source
projections, allocation strides, error propagation, finite-output rejection,
face exclusion, and diagnostic cadence control. Controlled rates and framework
shims do not validate a real EOS, distributed ownership, or framework scheduling.
No checked-in NRPyLeakageET test/parfile suite or full Toolkit runtime result is
supplied. The owning issue report keeps this validation gap open.

A configured ET validation should check fresh starts and checkpoint recovery,
multiple maps with different component owners, padded arrays, three-level AMR
with rejected truncated/reversed intervals, reflection/rotation scalar ghosts
immediately before POLR, PreSync enforcement/poisoning on and off, hydro cadences
0/1/greater than one, and startup-only selector rejection through a named steering
provider. Use independent face-length and covariant-source oracles plus successful,
early failing, and sanitized-output rate statuses. Kernel-only tests cannot
replace these wrapper and lifecycle checks.
