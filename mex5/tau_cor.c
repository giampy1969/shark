/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:44 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */
#include "tau_cor.h"
#include "vp.h"

/*
 * The function "Mtau_cor" is the implementation version of the "tau_cor"
 * M-function from file "E:\RTW\NEW3\TAU_COR.M" (lines 1-14). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function tc=tau_cor(veh,v,vr)
 */
static mxArray * Mtau_cor(int nargout_,
                          mxArray * veh,
                          mxArray * v,
                          mxArray * vr) {
    mxArray * tc = mclGetUninitializedArray();
    mxArray * Ca = mclGetUninitializedArray();
    mxArray * Crb = mclGetUninitializedArray();
    mclValidateInputs("tau_cor", 3, &veh, &v, &vr);
    /*
     * 
     * % tc=tau_cor(veh,v,vr) calculates coriolis forces from 
     * % vehicle variables and generalized velocities v and vr
     * 
     * 
     * Crb=[ zeros(3,3),           -vp(veh.Mrb(1:3,:)*v);
     */
    mlfAssign(
      &Crb,
      mlfVertcat(
        mlfHorzcat(
          mlfZeros(mlfScalar(3.0), mlfScalar(3.0), NULL),
          mlfUminus(
            mlfVp(
              mlfMtimes(
                mlfIndexRef(
                  veh,
                  ".Mrb(?,?)",
                  mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL),
                  mlfCreateColonIndex()),
                v),
              NULL)),
          NULL),
        mclCreateEmptyArray(),
        mlfHorzcat(
          mlfUminus(
            mlfVp(
              mlfMtimes(
                mlfIndexRef(
                  veh,
                  ".Mrb(?,?)",
                  mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL),
                  mlfCreateColonIndex()),
                v),
              NULL)),
          mlfUminus(
            mlfVp(
              mlfMtimes(
                mlfIndexRef(
                  veh,
                  ".Mrb(?,?)",
                  mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL),
                  mlfCreateColonIndex()),
                v),
              NULL)),
          NULL),
        NULL));
    /*
     * -vp(veh.Mrb(1:3,:)*v), -vp(veh.Mrb(4:6,:)*v)];
     * 
     * Ca= [ zeros(3,3),           -vp(veh.Ma(1:3,:)*vr);
     */
    mlfAssign(
      &Ca,
      mlfVertcat(
        mlfHorzcat(
          mlfZeros(mlfScalar(3.0), mlfScalar(3.0), NULL),
          mlfUminus(
            mlfVp(
              mlfMtimes(
                mlfIndexRef(
                  veh,
                  ".Ma(?,?)",
                  mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL),
                  mlfCreateColonIndex()),
                vr),
              NULL)),
          NULL),
        mclCreateEmptyArray(),
        mlfHorzcat(
          mlfUminus(
            mlfVp(
              mlfMtimes(
                mlfIndexRef(
                  veh,
                  ".Ma(?,?)",
                  mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL),
                  mlfCreateColonIndex()),
                vr),
              NULL)),
          mlfUminus(
            mlfVp(
              mlfMtimes(
                mlfIndexRef(
                  veh,
                  ".Ma(?,?)",
                  mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL),
                  mlfCreateColonIndex()),
                vr),
              NULL)),
          NULL),
        NULL));
    /*
     * -vp(veh.Ma(1:3,:)*vr), -vp(veh.Ma(4:6,:)*vr)];
     * 
     * tc=Crb*v+Ca*vr;
     */
    mlfAssign(&tc, mlfPlus(mlfMtimes(Crb, v), mlfMtimes(Ca, vr)));
    mclValidateOutputs("tau_cor", 1, nargout_, &tc);
    mxDestroyArray(Ca);
    mxDestroyArray(Crb);
    return tc;
}

/*
 * The function "mlfTau_cor" contains the normal interface for the "tau_cor"
 * M-function from file "E:\RTW\NEW3\TAU_COR.M" (lines 1-14). This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
mxArray * mlfTau_cor(mxArray * veh, mxArray * v, mxArray * vr) {
    int nargout = 1;
    mxArray * tc = mclGetUninitializedArray();
    mlfEnterNewContext(0, 3, veh, v, vr);
    tc = Mtau_cor(nargout, veh, v, vr);
    mlfRestorePreviousContext(0, 3, veh, v, vr);
    return mlfReturnValue(tc);
}

/*
 * The function "mlxTau_cor" contains the feval interface for the "tau_cor"
 * M-function from file "E:\RTW\NEW3\TAU_COR.M" (lines 1-14). The feval
 * function calls the implementation version of tau_cor through this function.
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxTau_cor(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: tau_cor Line: 1 Column: "
            "0 The function \"tau_cor\" was called with mor"
            "e than the declared number of outputs (1)"));
    }
    if (nrhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: tau_cor Line: 1 Column:"
            " 0 The function \"tau_cor\" was called with m"
            "ore than the declared number of inputs (3)"));
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 3 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 3; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    mplhs[0] = Mtau_cor(nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    plhs[0] = mplhs[0];
}
