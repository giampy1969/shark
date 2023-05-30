/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:42 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __rpy2J_h
#define __rpy2J_h 1

#include "matlab.h"

extern mxArray * mlfRpy2J(mxArray * rpy);
extern void mlxRpy2J(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#endif
