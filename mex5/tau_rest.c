/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:45 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */
#include "tau_rest.h"
#include "rpy2R_eb.h"
#include "vp.h"

/*
 * The function "Mtau_rest" is the implementation version of the "tau_rest"
 * M-function from file "E:\RTW\NEW3\TAU_REST.M" (lines 1-19). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function tr=tau_rest(veh,p)
 */
static mxArray * Mtau_rest(int nargout_, mxArray * veh, mxArray * p) {
    mxArray * tr = mclGetUninitializedArray();
    mxArray * FB_b = mclGetUninitializedArray();
    mxArray * FB_e = mclGetUninitializedArray();
    mxArray * FG_b = mclGetUninitializedArray();
    mxArray * FG_e = mclGetUninitializedArray();
    mxArray * MB_b = mclGetUninitializedArray();
    mxArray * MG_b = mclGetUninitializedArray();
    mxArray * tb = mclGetUninitializedArray();
    mxArray * tg = mclGetUninitializedArray();
    mclValidateInputs("tau_rest", 2, &veh, &p);
    /*
     * 
     * % tr=tau_rest(veh,p); calculates restoring forces from 
     * % vehicle variables and generalized position p
     * 
     * % Hydrostatic force and moment
     * FB_e=-veh.vol*veh.rho*veh.g_e;
     */
    mlfAssign(
      &FB_e,
      mlfFeval(
        mclValueVarargout(),
        mlxMtimes,
        mlfFeval(
          mclValueVarargout(),
          mlxMtimes,
          mlfFeval(
            mclValueVarargout(), mlxUminus, mlfIndexRef(veh, ".vol"), NULL),
          mlfIndexRef(veh, ".rho"),
          NULL),
        mlfIndexRef(veh, ".g_e"),
        NULL));
    /*
     * FB_b=rpy2R_eb(p(4:6))*FB_e;
     */
    mlfAssign(
      &FB_b,
      mlfMtimes(
        mlfRpy2R_eb(
          mlfIndexRef(
            p, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL))),
        FB_e));
    /*
     * MB_b=vp(veh.B_b,FB_b);
     */
    mlfAssign(
      &MB_b,
      mlfFeval(
        mclValueVarargout(), mlxVp, mlfIndexRef(veh, ".B_b"), FB_b, NULL));
    /*
     * tb=[FB_b;MB_b];
     */
    mlfAssign(
      &tb, mlfVertcat(mlfHorzcat(FB_b, NULL), mlfHorzcat(MB_b, NULL), NULL));
    /*
     * 
     * % Gravitational force and moment
     * FG_e=veh.m*veh.g_e;
     */
    mlfAssign(
      &FG_e,
      mlfFeval(
        mclValueVarargout(),
        mlxMtimes,
        mlfIndexRef(veh, ".m"),
        mlfIndexRef(veh, ".g_e"),
        NULL));
    /*
     * FG_b=rpy2R_eb(p(4:6))*FG_e;
     */
    mlfAssign(
      &FG_b,
      mlfMtimes(
        mlfRpy2R_eb(
          mlfIndexRef(
            p, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL))),
        FG_e));
    /*
     * MG_b=vp(veh.G_b,FG_b);
     */
    mlfAssign(
      &MG_b,
      mlfFeval(
        mclValueVarargout(), mlxVp, mlfIndexRef(veh, ".G_b"), FG_b, NULL));
    /*
     * tg=[FG_b;MG_b];
     */
    mlfAssign(
      &tg, mlfVertcat(mlfHorzcat(FG_b, NULL), mlfHorzcat(MG_b, NULL), NULL));
    /*
     * 
     * tr=tb+tg;
     */
    mlfAssign(&tr, mlfPlus(tb, tg));
    mclValidateOutputs("tau_rest", 1, nargout_, &tr);
    mxDestroyArray(FB_b);
    mxDestroyArray(FB_e);
    mxDestroyArray(FG_b);
    mxDestroyArray(FG_e);
    mxDestroyArray(MB_b);
    mxDestroyArray(MG_b);
    mxDestroyArray(tb);
    mxDestroyArray(tg);
    return tr;
}

/*
 * The function "mlfTau_rest" contains the normal interface for the "tau_rest"
 * M-function from file "E:\RTW\NEW3\TAU_REST.M" (lines 1-19). This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
mxArray * mlfTau_rest(mxArray * veh, mxArray * p) {
    int nargout = 1;
    mxArray * tr = mclGetUninitializedArray();
    mlfEnterNewContext(0, 2, veh, p);
    tr = Mtau_rest(nargout, veh, p);
    mlfRestorePreviousContext(0, 2, veh, p);
    return mlfReturnValue(tr);
}

/*
 * The function "mlxTau_rest" contains the feval interface for the "tau_rest"
 * M-function from file "E:\RTW\NEW3\TAU_REST.M" (lines 1-19). The feval
 * function calls the implementation version of tau_rest through this function.
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxTau_rest(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: tau_rest Line: 1 Column:"
            " 0 The function \"tau_rest\" was called with m"
            "ore than the declared number of outputs (1)"));
    }
    if (nrhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: tau_rest Line: 1 Column:"
            " 0 The function \"tau_rest\" was called with m"
            "ore than the declared number of inputs (2)"));
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
    mplhs[0] = Mtau_rest(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}
