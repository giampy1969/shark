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

#ifndef __rpy2J_h
#define __rpy2J_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_rpy2J(void);
extern void TerminateModule_rpy2J(void);
extern _mexLocalFunctionTable _local_function_table_rpy2J;

extern mxArray * mlfRpy2J(mxArray * rpy);
extern void mlxRpy2J(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
