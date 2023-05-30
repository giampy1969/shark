/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:47 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "matlab.h"
#include "a2clcdxc.h"
#include "a2clcd.h"
#include "vp.h"
#include "tau_rest.h"
#include "tau_damp.h"
#include "tau_cor.h"
#include "rpy2R_eb.h"
#include "rpy2J.h"
#include "vxdot.h"

mxArray * veh = NULL;

static mlfFunctionTableEntry function_table[9]
  = { { "a2clcdxc", mlxA2clcdxc, 1, 3 },
      { "a2clcd", mlxA2clcd, 1, 2 },
      { "vp", mlxVp, 2, 1 },
      { "tau_rest", mlxTau_rest, 2, 1 },
      { "tau_damp", mlxTau_damp, 3, 1 },
      { "tau_cor", mlxTau_cor, 3, 1 },
      { "rpy2R_eb", mlxRpy2R_eb, 1, 1 },
      { "rpy2J", mlxRpy2J, 1, 1 },
      { "vxdot", mlxVxdot, 1, 1 } };

static mclGlobalTableEntry global_table[1] = { { "veh", &veh } };

/*
 * The function "mexFunction" is a Compiler-generated mex wrapper, suitable for
 * building a MEX-function. It initializes any persistent variables as well as
 * a function table for use by the feval function. It then calls the function
 * "mlxVxdot". Finally, it clears the feval table and exits.
 */
void mexFunction(int nlhs, mxArray * * plhs, int nrhs, mxArray * * prhs) {
    mlfTry {
        mlfFunctionTableSetup(9, function_table);
        mclImportGlobal(1, global_table);
        mlxVxdot(nlhs, plhs, nrhs, prhs);
        mclExportGlobal(1, global_table);
        mlfFunctionTableTakedown(9, function_table);
    } mlfCatch {
        mclExportGlobal(1, global_table);
        mlfFunctionTableTakedown(9, function_table);
        mclMexError();
    } mlfEndCatch
}
