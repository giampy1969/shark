/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:43 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */
#include "rpy2R_eb.h"

/*
 * The function "Mrpy2R_eb" is the implementation version of the "rpy2R_eb"
 * M-function from file "E:\RTW\NEW3\rpy2R_eb.m" (lines 1-16). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function R_eb=rpy2R_eb(rpy)
 */
static mxArray * Mrpy2R_eb(int nargout_, mxArray * rpy) {
    mxArray * R_eb = mclGetUninitializedArray();
    mxArray * cf = mclGetUninitializedArray();
    mxArray * cp = mclGetUninitializedArray();
    mxArray * ct = mclGetUninitializedArray();
    mxArray * sf = mclGetUninitializedArray();
    mxArray * sp = mclGetUninitializedArray();
    mxArray * st = mclGetUninitializedArray();
    mclValidateInputs("rpy2R_eb", 1, &rpy);
    /*
     * 
     * % R_eb=rpy2R_eb(rpy), computes the rotation matrix of e-frame wrt b-frame,
     * % given in input a vector containing roll, pitch and yaw angles.
     * 
     * sf = sin(rpy(1));
     */
    mlfAssign(&sf, mlfSin(mlfIndexRef(rpy, "(?)", mlfScalar(1.0))));
    /*
     * cf = cos(rpy(1));
     */
    mlfAssign(&cf, mlfCos(mlfIndexRef(rpy, "(?)", mlfScalar(1.0))));
    /*
     * st = sin(rpy(2));
     */
    mlfAssign(&st, mlfSin(mlfIndexRef(rpy, "(?)", mlfScalar(2.0))));
    /*
     * ct = cos(rpy(2));
     */
    mlfAssign(&ct, mlfCos(mlfIndexRef(rpy, "(?)", mlfScalar(2.0))));
    /*
     * sp = sin(rpy(3));
     */
    mlfAssign(&sp, mlfSin(mlfIndexRef(rpy, "(?)", mlfScalar(3.0))));
    /*
     * cp = cos(rpy(3));
     */
    mlfAssign(&cp, mlfCos(mlfIndexRef(rpy, "(?)", mlfScalar(3.0))));
    /*
     * 
     * R_eb = [         +ct*cp               +ct*sp         -st
     */
    mlfAssign(
      &R_eb,
      mlfVertcat(
        mlfHorzcat(
          mlfMtimes(mlfUplus(ct), cp),
          mlfMtimes(mlfUplus(ct), sp),
          mlfUminus(st),
          NULL),
        mlfHorzcat(
          mlfMinus(
            mlfMtimes(mlfMtimes(mlfUplus(sf), st), cp), mlfMtimes(cf, sp)),
          mlfPlus(
            mlfMtimes(mlfMtimes(mlfUplus(sf), st), sp), mlfMtimes(cf, cp)),
          mlfMtimes(mlfUplus(sf), ct),
          NULL),
        mlfHorzcat(
          mlfPlus(
            mlfMtimes(mlfMtimes(mlfUplus(cf), st), cp), mlfMtimes(sf, sp)),
          mlfMinus(
            mlfMtimes(mlfMtimes(mlfUplus(cf), st), sp), mlfMtimes(sf, cp)),
          mlfMtimes(mlfUplus(cf), ct),
          NULL),
        NULL));
    mclValidateOutputs("rpy2R_eb", 1, nargout_, &R_eb);
    mxDestroyArray(cf);
    mxDestroyArray(cp);
    mxDestroyArray(ct);
    mxDestroyArray(sf);
    mxDestroyArray(sp);
    mxDestroyArray(st);
    /*
     * +sf*st*cp-cf*sp      +sf*st*sp+cf*cp      +sf*ct
     * +cf*st*cp+sf*sp      +cf*st*sp-sf*cp      +cf*ct ];
     */
    return R_eb;
}

/*
 * The function "mlfRpy2R_eb" contains the normal interface for the "rpy2R_eb"
 * M-function from file "E:\RTW\NEW3\rpy2R_eb.m" (lines 1-16). This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
mxArray * mlfRpy2R_eb(mxArray * rpy) {
    int nargout = 1;
    mxArray * R_eb = mclGetUninitializedArray();
    mlfEnterNewContext(0, 1, rpy);
    R_eb = Mrpy2R_eb(nargout, rpy);
    mlfRestorePreviousContext(0, 1, rpy);
    return mlfReturnValue(R_eb);
}

/*
 * The function "mlxRpy2R_eb" contains the feval interface for the "rpy2R_eb"
 * M-function from file "E:\RTW\NEW3\rpy2R_eb.m" (lines 1-16). The feval
 * function calls the implementation version of rpy2R_eb through this function.
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxRpy2R_eb(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: rpy2R_eb Line: 1 Column:"
            " 0 The function \"rpy2R_eb\" was called with m"
            "ore than the declared number of outputs (1)"));
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: rpy2R_eb Line: 1 Column:"
            " 0 The function \"rpy2R_eb\" was called with m"
            "ore than the declared number of inputs (1)"));
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 1 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 1; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 1, mprhs[0]);
    mplhs[0] = Mrpy2R_eb(nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
}
