#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "NRPyLeakageET.h"

#ifdef CCTK_MPI
#  include <mpi.h>
#endif

#define DISABLE_PROLONGATIONS EnableProlongating(0);
#define  ENABLE_PROLONGATIONS EnableProlongating(1);

using namespace Carpet;

extern "C"
int NRPyLeakageET_ProcessOwnsData() {
  // Carpet sets this from the current map/component using its communicator.
  return Carpet::local_component >= 0;
}

extern "C"
void NRPyLeakageET_Prolongate(CCTK_ARGUMENTS, const char *varname) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(verbosity_level>1) CCTK_VINFO("Prolongating %s from ref. lvl. %d to ref. lvl. %d...", varname, GetRefinementLevel(cctkGH)-1, GetRefinementLevel(cctkGH));

  int ml=Carpet::mglevel;
  int ll=Carpet::reflevel;
  for (comm_state state; not state.done(); state.step()) {
    const int g = CCTK_GroupIndex(varname);
    const int active_tl = CCTK_ActiveTimeLevelsGI (cctkGH, g);
    assert (active_tl>=0);
    const int tl = active_tl > 1 ? timelevel : 0;
    for (int m=0; m<(int)arrdata.AT(g).size(); ++m) {
      for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
        ggf *const gv = arrdata.AT(g).AT(m).data.AT(v);
        gv->ref_prolongate_all (state, tl, ll, ml, cctk_time);
        gv->ref_bnd_prolongate_all (state, tl, ll, ml, cctk_time);
      }
    }
  } // for state
}

extern "C"
void NRPyLeakageET_Restrict(CCTK_ARGUMENTS, const char *varname) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(verbosity_level>1) CCTK_VINFO("Restricting %s from ref. lvl. %d to ref. lvl. %d...", varname, GetRefinementLevel(cctkGH)+1, GetRefinementLevel(cctkGH));

  int ml=Carpet::mglevel;
  int ll=Carpet::reflevel;
  for (comm_state state; not state.done(); state.step()) {
    const int g = CCTK_GroupIndex(varname);
    const int active_tl = CCTK_ActiveTimeLevelsGI (cctkGH, g);
    assert (active_tl>=0);
    const int tl = active_tl > 1 ? timelevel : 0;
    for (int m=0; m<(int)arrdata.AT(g).size(); ++m) {
      for (int v = 0; v < (int)arrdata.AT(g).AT(m).data.size(); ++v) {
        ggf *const gv = arrdata.AT(g).AT(m).data.AT(v);
        gv->ref_restrict_all (state, tl, ll, ml);
      }
    }
  } // for state
}

extern "C"
void NRPyLeakageET_SyncOpticalDepths(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( verbosity_level > 1 ) CCTK_VINFO("Synchronizing auxiliary optical depths at ref. lev. %d",GetRefinementLevel(cctkGH));

  DISABLE_PROLONGATIONS;
  const int status = CCTK_SyncGroup(cctkGH,"NRPyLeakageET::NRPyLeakageET_optical_depths");
  if( status < 0 ) CCTK_VERROR("Could not synchronize NRPyLeakageET::NRPyLeakageET_optical_depths.");
  ENABLE_PROLONGATIONS;

  if( verbosity_level > 1 ) CCTK_VINFO("Finished synchronizing auxiliary optical depths at ref. lev. %d",GetRefinementLevel(cctkGH));
}

extern "C"
void NRPyLeakageET_SelectBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  const char *groups[2] = {"NRPyLeakageET::NRPyLeakageET_optical_depths",
                          "NRPyLeakageET::NRPyLeakageET_opacities"};
  for(int g=0;g<2;g++)
    if(Boundary_SelectGroupForBC(cctkGH,CCTK_ALL_FACES,1,-1,groups[g],"none") < 0)
      CCTK_VERROR("Could not select scalar symmetry boundaries for %s",groups[g]);
}

// Called in GLOBAL mode after the unscheduled helpers write the whole hierarchy.
static void NotifyInitializedLeakageData(CCTK_ARGUMENTS) {
  const char *groups[3] = {"NRPyLeakageET::NRPyLeakageET_optical_depths",
                          "NRPyLeakageET::NRPyLeakageET_opacities",
                          "NRPyLeakageET::NRPyLeakageET_auxiliary_optical_depths"};
  const int timelevels[3] = {3,3,1};
  CCTK_INT vars[42], tls[42], everywhere[42];
  int n=0;
  for(int g=0;g<3;g++) {
    const int group = CCTK_GroupIndex(groups[g]);
    if(group < 0 || CCTK_NumVarsInGroupI(group) != 6)
      CCTK_VERROR("Invalid leakage initialization group %s",groups[g]);
    const int first = CCTK_FirstVarIndexI(group);
    if(first < 0) CCTK_VERROR("Invalid first variable for %s",groups[g]);
    for(int tl=0;tl<timelevels[g];tl++) {
      for(int v=0;v<6;v++) {
        vars[n]=first+v; tls[n]=tl; everywhere[n]=CCTK_VALID_EVERYWHERE;
        n++;
      }
    }
  }
  if(Driver_NotifyDataModified(cctkGH,vars,tls,n,everywhere) != 0)
    CCTK_ERROR("Could not publish leakage initialization writes");
}

// LEVEL mode permits Carpet's ApplyBCs traversal to visit every local component.
static void ApplyLeakageSymmetries(CCTK_ARGUMENTS) {
  // These nested initialization helpers have no schedule entry of their own.
  // Fill refinement/component ghosts and scalar symmetries explicitly in every
  // PreSync mode, then publish the completed writes before checking validity.
  CCTK_INT vars[12], tls[12], interior[12], everywhere[12];
  const char *groups[2] = {"NRPyLeakageET::NRPyLeakageET_optical_depths",
                          "NRPyLeakageET::NRPyLeakageET_opacities"};
  for(int g=0;g<2;g++) {
    const int group = CCTK_GroupIndex(groups[g]);
    const int first = CCTK_FirstVarIndexI(group);
    if(group < 0 || first < 0 || CCTK_NumVarsInGroupI(group) != 6)
      CCTK_VERROR("Invalid leakage boundary group %s",groups[g]);
    for(int v=0;v<6;v++) {
      const int n=6*g+v;
      vars[n]=first+v; tls[n]=0;
      interior[n]=CCTK_VALID_INTERIOR; everywhere[n]=CCTK_VALID_EVERYWHERE;
    }
  }
  if(Driver_NotifyDataModified(cctkGH,vars,tls,12,interior) != 0)
    CCTK_ERROR("Could not publish leakage initialization interior writes");
  NRPyLeakageET_SelectBoundaries(CCTK_PASS_CTOC);
  ENABLE_PROLONGATIONS;
  for(int g=0;g<2;g++)
    if(CCTK_SyncGroup(cctkGH,groups[g]) < 0)
      CCTK_VERROR("Could not fill leakage initialization ghosts for %s",groups[g]);
  if(CCTK_ScheduleTraverse("ApplyBCs",const_cast<cGH *>(cctkGH),Carpet::CallFunction) != 0)
    CCTK_ERROR("Could not apply leakage scalar symmetry boundaries");
  if(Driver_NotifyDataModified(cctkGH,vars,tls,12,everywhere) != 0 ||
     Driver_RequireValidData(cctkGH,vars,tls,12,everywhere) != 0)
    CCTK_ERROR("Could not validate leakage initialization ghosts and symmetry zones");
}

extern "C"
void NRPyLeakageET_RegisterVars(CCTK_ARGUMENTS) {
  const char *groups[2] = {"NRPyLeakageET::NRPyLeakageET_opacities",
                          "NRPyLeakageET::NRPyLeakageET_optical_depths"};
  for(int g=0;g<2;g++) {
    const int index = CCTK_GroupIndex(groups[g]);
    if(index < 0 || MoLRegisterConstrainedGroup(index) != 0)
      CCTK_VERROR("Could not register constrained group %s",groups[g]);
  }
}

extern "C"
void NRPyLeakageET_getGlobalL2norm(CCTK_ARGUMENTS, const int startRefLev, const int endRefLev, const int it, CCTK_REAL *restrict l2norm_global) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Step 1: Loop over refinement levels, maps, and components, computing the
  // differences in the optical depths between two consecutive iterations
  for(int rl=startRefLev;rl<=endRefLev;rl++) {
    ENTER_LEVEL_MODE(cctkGH,rl) {
      BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
        BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          NRPyLeakageET_compute_optical_depth_change(CCTK_PASS_CTOC,it);
        } END_COMPONENT_LOOP;
      } END_MAP_LOOP;
    } LEAVE_LEVEL_MODE;
  }

  // Step 3: Get sum operation handle
  const int op_sum = CCTK_ReductionHandle("isum");

  // Step 4: Set variables we want to reduce
  const char varnames[6][30] = {"NRPyLeakageET::tau_0_nue_aux ",
                                "NRPyLeakageET::tau_1_nue_aux ",
                                "NRPyLeakageET::tau_0_anue_aux",
                                "NRPyLeakageET::tau_1_anue_aux",
                                "NRPyLeakageET::tau_0_nux_aux ",
                                "NRPyLeakageET::tau_1_nux_aux "};


  // Step 5: Get indices for the variables we want to reduce
  const int varindices[6] = {CCTK_VarIndex("NRPyLeakageET::tau_0_nue_aux"),
                             CCTK_VarIndex("NRPyLeakageET::tau_1_nue_aux"),
                             CCTK_VarIndex("NRPyLeakageET::tau_0_anue_aux"),
                             CCTK_VarIndex("NRPyLeakageET::tau_1_anue_aux"),
                             CCTK_VarIndex("NRPyLeakageET::tau_0_nux_aux"),
                             CCTK_VarIndex("NRPyLeakageET::tau_1_nux_aux")};

  // Step 6: Compute the global L2-norm:
  // l2norm = sqrt(dTau_0_nue + dTau_0_anue + dTau_0_nux
  //             + dTau_1_nue + dTau_1_anue + dTau_1_nux)
  CCTK_REAL l2norm = 0.0;
  for(int i=0;i<6;i++) {
    int varindex = varindices[i];
    assert(varindex>=0);
    CCTK_REAL l2normL = 0.0;
    const int ierr = CCTK_Reduce(cctkGH,
                                 -1, // target processors; -1 -> all
                                 op_sum,
                                 1,  // Number of outputs
                                 CCTK_VARIABLE_REAL,
                                 &l2normL,
                                 1,  // Number of inputs
                                 varindex);
    if( ierr ) CCTK_VERROR("Error in reduction of L2-norm of %s",varnames[i]);

    l2norm += l2normL;

    if( verbosity_level > 1 ) CCTK_VINFO("L2-norm of %s: %e",varnames[i],sqrt(l2normL));
  }
  *l2norm_global = sqrt(l2norm);

  // Step 7: Update auxiliary variables
  for(int rl=startRefLev;rl<=endRefLev;rl++) {
    ENTER_LEVEL_MODE(cctkGH,rl) {
      BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
        BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
          NRPyLeakageET_CopyOpticalDepthsToAux(CCTK_PASS_CTOC);
        } END_COMPONENT_LOOP;
      } END_MAP_LOOP;
    } LEAVE_LEVEL_MODE;
  }
}

extern "C"
void NRPyLeakageET_Print(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!NRPyLeakageET_ProcessOwnsData()) return;

  for(int i=0;i<12 && i<cctk_lsh[0] && i<cctk_lsh[1] && i<cctk_lsh[2];i++) {
    const int index = CCTK_GFINDEX3D(cctkGH,i,i,i);
    CCTK_VINFO("Ref. Lev. %d - %02d,%02d,%02d - %g %g %g %g %g %g",
               GetRefinementLevel(cctkGH),i,i,i,
               tau_0_nue[index],tau_0_anue[index],tau_0_nux[index],
               tau_1_nue[index],tau_1_anue[index],tau_1_nux[index]);
  }
}

extern "C"
void NRPyLeakageET_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // The helpers initialize all three stored timelevels, including past arrays.
  if(!Carpet::do_allow_past_timelevels)
    CCTK_ERROR("NRPyLeakageET initialization requires accessible past timelevels; use InitBase::initial_data_setup_method=init_some_levels with Carpet::init_each_timelevel=no");

  // Step 1: Initialize all optical depth gridfunctions to zero
  if(verbosity_level>0) CCTK_INFO("Initializing optical depths gridfunctions to zero...");
  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        NRPyLeakageET_optical_depths_initialize_to_zero(CCTK_PASS_CTOC);
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
  } END_REFLEVEL_LOOP;
  NotifyInitializedLeakageData(CCTK_PASS_CTOC);
  if(verbosity_level>0) CCTK_INFO("Initialized all optical depths gridfunctions to zero");

  if( CCTK_EQUALS(initial_optical_depth,"PathOfLeastResistance") ) {

    const int startRefLev = ghl_imin(ghl_imax(minInitRefLevel,0),Carpet::reflevels-1);
    const int endRefLev   = maxInitRefLevel == 0 ? Carpet::reflevels-1 : ghl_imin(maxInitRefLevel,Carpet::reflevels-1);

    if(startRefLev > endRefLev)
      CCTK_VERROR("Empty POLR initialization interval: resolved levels %d through %d",startRefLev,endRefLev);
    if(endRefLev != Carpet::reflevels-1)
      CCTK_ERROR("POLR initialization must include the finest active level; truncated windows are unsupported by the composite convergence reduction");

    if(verbosity_level>0) {
      CCTK_VINFO("Optical depths will be initialized on levels %d through %d with the path of least resistance algorithm",startRefLev,endRefLev);
      CCTK_VINFO("Number of iterations to be performed on each refinement level: %d", max_iterations);
    }

    // Step 2: Now perform iterations of the path of least resistance algorithm
    int counter = 0;
    int RemainingIterations = max_iterations;

    for(int iteration=0;iteration<max_iterations;iteration++) {
      // Step 2.a: First perform the iteration on refinement level startRefLev
      if( verbosity_level>1 ) CCTK_VINFO("Beginning iteration %d...",counter+1);

      // Step 2.b: Now loop over remaining refinement levels, performing the
      //           path of least resistance (POLR) algorithm
      for(int rl=startRefLev;rl<=endRefLev;rl++) {
        ENTER_LEVEL_MODE(cctkGH,rl) {
          // Step 2.b.i: Prolongate results from the previous refinement level
          if( rl > startRefLev ) NRPyLeakageET_Prolongate(CCTK_PASS_CTOC,"NRPyLeakageET::NRPyLeakageET_optical_depths");
          ApplyLeakageSymmetries(CCTK_PASS_CTOC);
          BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
            BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
              NRPyLeakageET_compute_neutrino_opacities(CCTK_PASS_CTOC);
            } END_COMPONENT_LOOP;
          } END_MAP_LOOP;
          ApplyLeakageSymmetries(CCTK_PASS_CTOC);
          BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
            BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
              NRPyLeakageET_optical_depths_PathOfLeastResistance(CCTK_PASS_CTOC);
            } END_COMPONENT_LOOP;
          } END_MAP_LOOP;
          // Step 2.b.iii: Synchronize optical depths gridfunctions
          //               at this refinement level
          NRPyLeakageET_SyncOpticalDepths(CCTK_PASS_CTOC);
          ApplyLeakageSymmetries(CCTK_PASS_CTOC);
        } LEAVE_LEVEL_MODE;
      }

      // Step 2.c: Perform a reduction of the L2-norm among all processors
      CCTK_REAL l2norm_global;
      NRPyLeakageET_getGlobalL2norm(CCTK_PASS_CTOC, startRefLev, endRefLev, iteration+1, &l2norm_global);

      // Step 2.d: Increment counter; check convergence
      counter++;
      if( l2norm_global < tauChangeThreshold ) {
        // Step 2.d.i: Converged!
        RemainingIterations = 0;
        if( verbosity_level>0 ) CCTK_VINFO("Algorithm converged! l2norm_tau_change = %e (threshold = %e)",l2norm_global,tauChangeThreshold);
        break;
      }
      else {
        // Step 2.d.ii: Have not converged. Restrict current solution from the
        // finest to the coarsest initialized refinement levels; Decrement
        // number of iterations remaining; continue.
        if( RemainingIterations > 1 ) {
          for(int rl=endRefLev-1;rl>=startRefLev;rl--) {
            ENTER_LEVEL_MODE(cctkGH,rl) {
              NRPyLeakageET_Restrict(CCTK_PASS_CTOC,"NRPyLeakageET::NRPyLeakageET_optical_depths");
            } LEAVE_LEVEL_MODE;
          }
          RemainingIterations--;
          if( verbosity_level>0 ) CCTK_VINFO("Completed iteration %d. Global change: %e. Remaining iterations: %4d",counter,l2norm_global,RemainingIterations);
        }
        else {
          CCTK_WARN(CCTK_WARN_ALERT,"Maximum iterations reached, but convergence threshold not met.");
        }
      }
    }

    // Step 3: All done! Now loop over the refinement levels backwards, restricting the result
    for(int rl=endRefLev-1;rl>=0;rl--) {
      ENTER_LEVEL_MODE(cctkGH,rl) {
        NRPyLeakageET_Restrict(CCTK_PASS_CTOC,"NRPyLeakageET::NRPyLeakageET_optical_depths");
      } LEAVE_LEVEL_MODE;
    }
    // Restriction invalidates refinement and symmetry ghosts on both levels.
    // Refill coarse to fine before adopting current data as past timelevels.
    for(int rl=0;rl<Carpet::reflevels;rl++) {
      ENTER_LEVEL_MODE(cctkGH,rl) {
        ApplyLeakageSymmetries(CCTK_PASS_CTOC);
      } LEAVE_LEVEL_MODE;
    }
    if(verbosity_level>0) CCTK_INFO("Completed path of least resistance algorithm!");
  }

  // Step 4: Now copy the optical depths and opacities to all time levels
  if(verbosity_level>0) CCTK_INFO("Copying initial data to all time levels...");
  BEGIN_REFLEVEL_LOOP(cctkGH) {
    BEGIN_MAP_LOOP(cctkGH,CCTK_GF) {
      BEGIN_COMPONENT_LOOP(cctkGH, CCTK_GF) {
        NRPyLeakageET_copy_opacities_and_optical_depths_to_previous_time_levels(CCTK_PASS_CTOC);
      } END_COMPONENT_LOOP;
    } END_MAP_LOOP;
  } END_REFLEVEL_LOOP;
  NotifyInitializedLeakageData(CCTK_PASS_CTOC);
  if(verbosity_level>0) CCTK_INFO("Finished copying initial data to all time levels");

}
