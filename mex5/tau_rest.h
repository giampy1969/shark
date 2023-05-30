/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:45 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __tau_rest_h
#define __tau_rest_h 1

#include "matlab.h"

extern mxArray * mlfTau_rest(mxArray * veh, mxArray * p);
extern void mlxTau_rest(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#endif
