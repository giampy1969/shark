/*
 * MATLAB Compiler: 2.1
 * Date: Wed Oct 17 16:15:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-h" "-x" "-W" "mex"
 * "-L" "C" "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "vxdot.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "libmatlb.h"
#include "vxdot.h"
#include "rpy2J.h"
#include "rpy2R_eb.h"
#include "tau_cor.h"
#include "tau_damp.h"
#include "tau_rest.h"
#include "vp.h"
#include "a2clcd.h"
#include "a2clcdxc.h"

mxArray * veh = NULL;

static mexGlobalTableEntry global_table[1] = { { "veh", &veh } };

static mexFunctionTableEntry function_table[9]
  = { { "vxdot", mlxVxdot, 1, 1, &_local_function_table_vxdot },
      { "rpy2J", mlxRpy2J, 1, 1, &_local_function_table_rpy2J },
      { "rpy2R_eb", mlxRpy2R_eb, 1, 1, &_local_function_table_rpy2R_eb },
      { "tau_cor", mlxTau_cor, 3, 1, &_local_function_table_tau_cor },
      { "tau_damp", mlxTau_damp, 3, 1, &_local_function_table_tau_damp },
      { "tau_rest", mlxTau_rest, 2, 1, &_local_function_table_tau_rest },
      { "vp", mlxVp, 2, 1, &_local_function_table_vp },
      { "a2clcd", mlxA2clcd, 1, 2, &_local_function_table_a2clcd },
      { "a2clcdxc", mlxA2clcdxc, 1, 3, &_local_function_table_a2clcdxc } };

static _mexInitTermTableEntry init_term_table[9]
  = { { InitializeModule_vxdot, TerminateModule_vxdot },
      { InitializeModule_rpy2J, TerminateModule_rpy2J },
      { InitializeModule_rpy2R_eb, TerminateModule_rpy2R_eb },
      { InitializeModule_tau_cor, TerminateModule_tau_cor },
      { InitializeModule_tau_damp, TerminateModule_tau_damp },
      { InitializeModule_tau_rest, TerminateModule_tau_rest },
      { InitializeModule_vp, TerminateModule_vp },
      { InitializeModule_a2clcd, TerminateModule_a2clcd },
      { InitializeModule_a2clcdxc, TerminateModule_a2clcdxc } };

static _mex_information _mex_info
  = { 1, 9, function_table, 1, global_table, 0, NULL, 9, init_term_table };

/*
 * The function "mexLibrary" is a Compiler-generated mex wrapper, suitable for
 * building a MEX-function. It initializes any persistent variables as well as
 * a function table for use by the feval function. It then calls the function
 * "mlxVxdot". Finally, it clears the feval table and exits.
 */
mex_information mexLibrary(void) {
    return &_mex_info;
}
