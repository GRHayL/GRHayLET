# GRHayLib Contract

> Page status: reviewed · Last reviewed: 10-04-2026
> Up: [Integration](index.md)

## Scope and Non-Scope

This page owns the local declarations and visible dataflow described below.
Only GRHayLID sources are domain evidence. Framework execution, library
semantics, production-table provenance, and numerical validation remain
out of scope.

## Summary

Local CCL shares EOS_type from GRHayLib, includes GRHayLib.h, and requires
HDF5. The header rejects non-eight-byte Cactus real configurations. All
output-pointer EOS calls use local scalar temporaries, typed CCTK_REAL in the
gas and sphere initializers and double in the one-dimensional hydro, beta, and
entropy code. Beta uses interpolated
base potentials at each actual rho/T rather than the cached-root APIs.

## Mode Applicability

| Applicability | Local surface |
| --- | --- |
| Common | External initialized ghl_eos and EOS dispatch remain prerequisites. |
| HydroTest1D | Simple metadata or Hybrid cold/energy helpers. |
| BetaEquilibrium | Base-potential and pressure/energy APIs. |
| Entropy/Hybrid | Native entropy helper. |
| Entropy/Tabulated | Bounds and pressure/energy/entropy helpers. |

## Claim-Evidence

| Claim ID | Claim | Status | Evidence | Typed locator |
| --- | --- | --- | --- | --- |
| `GRHAYLIB-CONTRACT-01` | EOS_type is used from GRHayLib rather than declared locally. | declared | Named local source | `ccl:GRHayLID/param.ccl#parameter=EOS_type` |
| `GRHAYLIB-CONTRACT-02` | Header includes GRHayLib.h and enforces CCTK_REAL_PRECISION_8. | visible-implementation | Named local source | `macro:GRHayLID/src/GRHayLID.h#include=GRHayLib.h` |
| `GRHAYLIB-CONTRACT-03` | HDF5 is a build requirement. | declared | Named local source | `ccl:GRHayLID/configuration.ccl#requirement=HDF5` |
| `GRHAYLIB-CONTRACT-04` | Beta calls the base-potential API and checks its status/finiteness. | visible-implementation | Named local source | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_beta_residual` |
| `GRHAYLIB-CONTRACT-05` | Tabulated entropy uses double temporaries and checks the returned status. | visible-implementation | Named local source | `c:GRHayLID/src/ComputeEntropy.c#symbol=GRHayLID_compute_entropy_tabulated` |
| `GRHAYLIB-CONTRACT-06` | External EOS handle metadata, initialization, ABI, and table semantics are unverified locally. | out-of-scope | Named local source | `c:GRHayLID/src/BetaEquilibrium.c#symbol=GRHayLID_BetaEquilibrium` |

## Details

| Local API family | Visible names |
| --- | --- |
| Hybrid cold/energy | ghl_hybrid_compute_P_cold_and_eps_cold; ghl_hybrid_compute_epsilon |
| Native entropy | ghl_hybrid_compute_entropy_function |
| Tabulated thermodynamics | ghl_tabulated_compute_P_eps_from_T; ghl_tabulated_compute_P_eps_S_from_T |
| Base chemical potentials | ghl_tabulated_compute_P_eps_muhat_mue_mup_mun_from_T |
| Effective bounds | ghl_tabulated_enforce_bounds_rho_Ye_T |

Local code checks statuses against ghl_success and stages all EOS output
pointers through local scalar temporaries typed CCTK_REAL or double. The
header requires `CCTK_REAL` to be Cactus's eight-byte real (`double`); that
requirement also bounds the GRHayLib parameter-array interface and does not
establish a complete Cactus build. No pointer cast is used as a representation
bridge.

The local beta solver reads effective rho/Ye/T bounds, atmosphere metadata,
N_Ye, and table_Y_e. It no longer calls the cached-root builder or its density
interpolator and does not modify shared beta cache arrays. The selected library
must initialize the required function pointers, table nodes, and bounds.
Physical provenance and correctness of base potentials are external obligations.

## Caveats

Storage and schedule declarations do not prove allocation or execution.
Local guards and calls do not establish external error or interpolation
semantics. No checked-in GRHayLID test/parfile/oracle validates this path.

## Sources

- [param.ccl](../../../GRHayLID/param.ccl)
- [GRHayLID.h](../../../GRHayLID/src/GRHayLID.h)
- [configuration.ccl](../../../GRHayLID/configuration.ccl)
- [BetaEquilibrium.c](../../../GRHayLID/src/BetaEquilibrium.c)
- [ComputeEntropy.c](../../../GRHayLID/src/ComputeEntropy.c)

## Related Pages

- [Keyword Extensions](hydrobase-keyword-extensions.md)
- [Parameters](parameters-and-configurations.md)
- [Beta](../initial-data/beta-equilibrium.md)
- [Entropy](../initial-data/entropy-computation.md)
