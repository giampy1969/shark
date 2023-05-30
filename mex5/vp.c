/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:46 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */
#include "vp.h"

/*
 * The function "Mvp" is the implementation version of the "vp" M-function from
 * file "E:\RTW\NEW3\VP.M" (lines 1-13). It contains the actual compiled code
 * for that M-function. It is a static function and must only be called from
 * one of the interface functions, appearing below.
 */
/*
 * function z=vp(x,y)
 */
static mxArray * Mvp(int nargout_, mxArray * x, mxArray * y) {
    mxArray * z = mclGetUninitializedArray();
    mxArray * nargin_ = mclGetUninitializedArray();
    mlfAssign(&nargin_, mlfNargin(0, x, y, NULL));
    mclValidateInputs("vp", 2, &x, &y);
    /*
     * 
     * % z=vp(x,y); z = 3d cross product of x and y
     * % vp(x) is the 3d cross product matrix : vp(x)*y=vp(x,y).
     * %
     * % by Giampiero Campa.  
     * 
     * z=[  0    -x(3)   x(2);
     */
    mlfAssign(
      &z,
      mlfVertcat(
        mlfHorzcat(
          mlfScalar(0.0),
          mlfUminus(mlfIndexRef(x, "(?)", mlfScalar(3.0))),
          mlfIndexRef(x, "(?)", mlfScalar(2.0)),
          NULL),
        mclCreateEmptyArray(),
        mlfHorzcat(
          mlfIndexRef(x, "(?)", mlfScalar(3.0)),
          mlfScalar(0.0),
          mlfUminus(mlfIndexRef(x, "(?)", mlfScalar(1.0))),
          NULL),
        mclCreateEmptyArray(),
        mlfHorzcat(
          mlfUminus(mlfIndexRef(x, "(?)", mlfScalar(2.0))),
          mlfIndexRef(x, "(?)", mlfScalar(1.0)),
          mlfScalar(0.0),
          NULL),
        NULL));
    /*
     * x(3)    0    -x(1);
     * -x(2)   x(1)    0   ];
     * 
     * if nargin>1, z=z*y; end
     */
    if (mlfTobool(mlfGt(nargin_, mlfScalar(1.0)))) {
        mlfAssign(&z, mlfMtimes(z, y));
    }
    mclValidateOutputs("vp", 1, nargout_, &z);
    mxDestroyArray(nargin_);
    return z;
}

/*
 * The function "mlfVp" contains the normal interface for the "vp" M-function
 * from file "E:\RTW\NEW3\VP.M" (lines 1-13). This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
mxArray * mlfVp(mxArray * x, mxArray * y) {
    int nargout = 1;
    mxArray * z = mclGetUninitializedArray();
    mlfEnterNewContext(0, 2, x, y);
    z = Mvp(nargout, x, y);
    mlfRestorePreviousContext(0, 2, x, y);
    return mlfReturnValue(z);
}

/*
 * The function "mlxVp" contains the feval interface for the "vp" M-function
 * from file "E:\RTW\NEW3\VP.M" (lines 1-13). The feval function calls the
 * implementation version of vp through this function. This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
void mlxVp(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: vp Line: 1 Column: 0 The function \"vp\" "
            "was called with more than the declared number of outputs (1)"));
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: vp Line: 1 Column: 0 The function \"vp\" "
            "was called with more than the declared number of inputs (2)"));
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 2 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 2; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 2, mprhs[0], mprhs[1]);
    mplhs[0] = Mvp(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}
