/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:44 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __tau_cor_h
#define __tau_cor_h 1

#include "matlab.h"

extern mxArray * mlfTau_cor(mxArray * veh, mxArray * v, mxArray * vr);
extern void mlxTau_cor(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#endif
