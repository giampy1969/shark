/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:42 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __vxdot_h
#define __vxdot_h 1

#include "matlab.h"

extern mxArray * mlfVxdot(mxArray * xu);
extern void mlxVxdot(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#endif
