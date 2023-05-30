/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:43 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */
#include "rpy2J.h"
#include "rpy2R_eb.h"

/*
 * The function "Mrpy2J" is the implementation version of the "rpy2J"
 * M-function from file "E:\RTW\NEW3\rpy2J.m" (lines 1-14). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function J=rpy2J(rpy)
 */
static mxArray * Mrpy2J(int nargout_, mxArray * rpy) {
    mxArray * J = mclGetUninitializedArray();
    mxArray * cf = mclGetUninitializedArray();
    mxArray * ct = mclGetUninitializedArray();
    mxArray * sf = mclGetUninitializedArray();
    mxArray * tt = mclGetUninitializedArray();
    mclValidateInputs("rpy2J", 1, &rpy);
    /*
     * 
     * % J=rpy2J(rpy); computes generalised Jacobian matrix which 
     * % transforms ni into eta derivatives, given in input roll pitch
     * % and yaw angles.
     * 
     * sf = sin(rpy(1));
     */
    mlfAssign(&sf, mlfSin(mlfIndexRef(rpy, "(?)", mlfScalar(1.0))));
    /*
     * cf = cos(rpy(1));
     */
    mlfAssign(&cf, mlfCos(mlfIndexRef(rpy, "(?)", mlfScalar(1.0))));
    /*
     * tt = tan(rpy(2));
     */
    mlfAssign(&tt, mlfTan(mlfIndexRef(rpy, "(?)", mlfScalar(2.0))));
    /*
     * ct = cos(rpy(2));
     */
    mlfAssign(&ct, mlfCos(mlfIndexRef(rpy, "(?)", mlfScalar(2.0))));
    /*
     * 
     * J = [ rpy2R_eb(rpy)' zeros(3,3)
     */
    mlfAssign(
      &J,
      mlfVertcat(
        mlfHorzcat(
          mlfCtranspose(mlfRpy2R_eb(rpy)),
          mlfZeros(mlfScalar(3.0), mlfScalar(3.0), NULL),
          NULL),
        mlfHorzcat(
          mlfZeros(mlfScalar(3.0), mlfScalar(3.0), NULL),
          mlfVertcat(
            mlfHorzcat(
              mlfScalar(1.0), mlfMtimes(sf, tt), mlfMtimes(cf, tt), NULL),
            mlfHorzcat(mlfScalar(0.0), cf, mlfUminus(sf), NULL),
            mlfHorzcat(
              mlfScalar(0.0), mlfMrdivide(sf, ct), mlfMrdivide(cf, ct), NULL),
            NULL),
          NULL),
        NULL));
    mclValidateOutputs("rpy2J", 1, nargout_, &J);
    mxDestroyArray(cf);
    mxDestroyArray(ct);
    mxDestroyArray(sf);
    mxDestroyArray(tt);
    /*
     * zeros(3,3)     [1 sf*tt cf*tt; 0 cf -sf; 0 sf/ct cf/ct]  ];
     */
    return J;
}

/*
 * The function "mlfRpy2J" contains the normal interface for the "rpy2J"
 * M-function from file "E:\RTW\NEW3\rpy2J.m" (lines 1-14). This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
mxArray * mlfRpy2J(mxArray * rpy) {
    int nargout = 1;
    mxArray * J = mclGetUninitializedArray();
    mlfEnterNewContext(0, 1, rpy);
    J = Mrpy2J(nargout, rpy);
    mlfRestorePreviousContext(0, 1, rpy);
    return mlfReturnValue(J);
}

/*
 * The function "mlxRpy2J" contains the feval interface for the "rpy2J"
 * M-function from file "E:\RTW\NEW3\rpy2J.m" (lines 1-14). The feval function
 * calls the implementation version of rpy2J through this function. This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
void mlxRpy2J(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: rpy2J Line: 1 Column: 0 The function \"rpy2J"
            "\" was called with more than the declared number of outputs (1)"));
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: rpy2J Line: 1 Column: 0 The function \"rpy2J"
            "\" was called with more than the declared number of inputs (1)"));
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
    mplhs[0] = Mrpy2J(nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
}
