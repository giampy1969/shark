/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:46 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __a2clcdxc_h
#define __a2clcdxc_h 1

#include "matlab.h"

extern mxArray * mlfA2clcdxc(mxArray * * Cd, mxArray * * xcp, mxArray * alfa);
extern void mlxA2clcdxc(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#endif
