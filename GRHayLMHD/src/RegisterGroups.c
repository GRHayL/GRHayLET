#include "GRHayLMHD.h"

#define NUM_EVOLVED_GROUPS     (1)
#define NUM_CONSTRAINED_GROUPS (9)
#define NUM_RESTORED_GROUPS    (4)

typedef struct {
    const char *var, *rhs;
} group_t;

void GRHayLMHD_RegisterGroups(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const group_t evolved_groups[NUM_EVOLVED_GROUPS] = {
        { "GRHayLMHD::hydro_conservatives", "GRHayLMHD::hydro_rhss" },
        // { "GRHayLMHD::Phi_tilde", "GRHayLMHD::Phi_rhs" },
        // { "GRHayLMHD::A_x_tilde", "GRHayLMHD::A_x_rhs" },
        // { "GRHayLMHD::A_y_tilde", "GRHayLMHD::A_y_rhs" },
        // { "GRHayLMHD::A_z_tilde", "GRHayLMHD::A_z_rhs" },
    };

    for(int n = 0; n < NUM_EVOLVED_GROUPS; n++) {
        const char *var = evolved_groups[n].var;
        const char *rhs = evolved_groups[n].rhs;
        if(MoLRegisterEvolvedGroup(CCTK_GroupIndex(var), CCTK_GroupIndex(rhs))) {
            CCTK_VERROR("Problem registering evolved group (%s, %s)", var, rhs);
        }
    }

    const char *constrained_groups[NUM_CONSTRAINED_GROUPS] = {
        "HydroBase::rho",
        "HydroBase::press",
        "HydroBase::eps",
        "HydroBase::Y_e",
        "HydroBase::temperature",
        "HydroBase::entropy",
        "TmunuBase::stress_energy_scalar",
        "TmunuBase::stress_energy_vector",
        "TmunuBase::stress_energy_tensor",
    };

    for(int n = 0; n < NUM_CONSTRAINED_GROUPS; n++) {
        if(MoLRegisterConstrainedGroup(CCTK_GroupIndex(constrained_groups[n]))) {
            CCTK_VERROR("Problem registering constrained group %s", constrained_groups[n]);
        }
    }

    const char *restored_groups[NUM_RESTORED_GROUPS] = {
        "ADMBase::lapse",
        "ADMBase::shift",
        "ADMBase::metric",
        "ADMBase::curv",
    };

    for(int n = 0; n < NUM_RESTORED_GROUPS; n++) {
        if(MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex(restored_groups[n]))) {
            CCTK_VERROR("Problem registering restored group %s", restored_groups[n]);
        }
    }
}
