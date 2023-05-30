/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:44 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __tau_damp_h
#define __tau_damp_h 1

#include "matlab.h"

extern mxArray * mlfTau_damp(mxArray * veh, mxArray * vr, mxArray * de);
extern void mlxTau_damp(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#endif
