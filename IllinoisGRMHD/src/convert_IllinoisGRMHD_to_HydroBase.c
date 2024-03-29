#include "IllinoisGRMHD.h"

void convert_IllinoisGRMHD_to_HydroBase(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_convert_IllinoisGRMHD_to_HydroBase;
  DECLARE_CCTK_PARAMETERS;

  // if/else for backward compatibility
  if(CCTK_IsThornActive("Convert_to_HydroBase")) {
    int partype;
    void const *const parptr = CCTK_ParameterGet("Convert_to_HydroBase_every", "Convert_to_HydroBase", &partype);
    const int old_Convert_to_HydroBase_every = *(CCTK_INT const *)parptr;
    if(old_Convert_to_HydroBase_every==0) return;
    if(cctk_iteration%old_Convert_to_HydroBase_every!=0) return;
  } else {
    // Generally, we only need the HydroBase variables for diagnostic purposes, so we run the below loop only at iterations in which diagnostics are run.
    if(cctk_iteration%Convert_to_HydroBase_every!=0) return;
  }

  const CCTK_REAL mag_factor = rescale_magnetics ? sqrt(4.0*M_PI) : 1;

#pragma omp parallel for
  for(int k=0; k<cctk_lsh[2]; k++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int i=0; i<cctk_lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        const int index4D0 = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,0);
        const int index4D1 = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,1);
        const int index4D2 = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,2);

        // IllinoisGRMHD defines v^i = u^i/u^0.

        // Meanwhile, the ET/HydroBase formalism, called the Valencia
        // formalism, splits the 4 velocity into a purely spatial part
        // and a part that is normal to the spatial hypersurface:
        // u^a = G (n^a + U^a), (Eq. 14 of arXiv:1304.5544; G=W, U^a=v^a)
        // where n^a is the unit normal vector to the spatial hypersurface,
        // n_a = {-\alpha,0,0,0}, and U^a is the purely spatial part, which
        // is defined in HydroBase as the vel[] vector gridfunction.
        // Then u^a n_a = - \alpha u^0 = G n^a n_a = -G, and
        // of course \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma,
        // the standard Lorentz factor.

        // Note that n^i = - \beta^i / \alpha, so
        // u^a = \Gamma (n^a + U^a)
        // -> u^i = \Gamma ( U^i - \beta^i / \alpha )
        // which implies
        // v^i = u^i/u^0
        //     = \Gamma/u^0 ( U^i - \beta^i / \alpha ) <- \Gamma = \alpha u^0
        //     = \alpha ( U^i - \beta^i / \alpha )
        //     = \alpha U^i - \beta^i
        const CCTK_REAL lapseL = alp[index];
        const CCTK_REAL lapseL_inv = 1.0/lapseL;
        const CCTK_REAL utU[3] = {vx[index] + betax[index],
                               vy[index] + betay[index],
                               vz[index] + betaz[index]};

        vel[index4D0] = utU[0]*lapseL_inv;
        vel[index4D1] = utU[1]*lapseL_inv;
        vel[index4D2] = utU[2]*lapseL_inv;

        // \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma = w_lorentz
        // First compute u^0:
        // Derivation of first equation:
        // \gamma_{ij} (v^i + \beta^i)(v^j + \beta^j)/(\alpha)^2
        //   = \gamma_{ij} 1/(u^0)^2 ( \gamma^{ik} u_k \gamma^{jl} u_l /(\alpha)^2 <- Using Eq. 53 of arXiv:astro-ph/0503420
        //   = 1/(u^0 \alpha)^2 u_j u_l \gamma^{jl}  <- Since \gamma_{ij} \gamma^{ik} = \delta^k_j
        //   = 1/(u^0 \alpha)^2 ( (u^0 \alpha)^2 - 1 ) <- Using Eq. 56 of arXiv:astro-ph/0503420
        //   = 1 - 1/(u^0 \alpha)^2 <= 1
        const CCTK_REAL gxxL = gxx[index];
        const CCTK_REAL gxyL = gxy[index];
        const CCTK_REAL gxzL = gxz[index];
        const CCTK_REAL gyyL = gyy[index];
        const CCTK_REAL gyzL = gyz[index];
        const CCTK_REAL gzzL = gzz[index];

        const CCTK_REAL one_minus_invW_squared =
              (gxxL* SQR(utU[0]) +
               gyyL* SQR(utU[1]) + 
               gzzL* SQR(utU[2]) +
               2.0*(gxyL*(utU[0])*(utU[1]) +
                    gxzL*(utU[0])*(utU[2]) +
                    gyzL*(utU[1])*(utU[2])
               ))*SQR(lapseL_inv);
        /*** Check for superluminal velocity ***/
        //FIXME: Instead of >1.0, should be one_minus_one_over_alpha_u0_squared > ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED, for consistency with conserv_to_prims routines

        if(one_minus_invW_squared > 1.0) {
          CCTK_VINFO("convert_from_IllinoisGRMHD_to_HydroBase WARNING: Found superluminal velocity. This should have been caught by IllinoisGRMHD.");
        }

        const CCTK_REAL W = 1.0/sqrt(1.0-one_minus_invW_squared);
        if(isnan(W*lapseL_inv)) CCTK_VINFO("BAD FOUND NAN ALPHAU0 CALC: %.15e %.15e %.15e\n", W, lapseL_inv, one_minus_invW_squared);

        w_lorentz[index] = W;

        Bvec[index4D0] = Bx_center[index]*mag_factor;
        Bvec[index4D1] = By_center[index]*mag_factor;
        Bvec[index4D2] = Bz_center[index]*mag_factor;
      }
    }
  }
}
