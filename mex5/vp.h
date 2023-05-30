/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:45 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __vp_h
#define __vp_h 1

#include "matlab.h"

extern mxArray * mlfVp(mxArray * x, mxArray * y);
extern void mlxVp(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#endif
