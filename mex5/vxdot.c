/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:42 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */
#include "vxdot.h"
#include "rpy2J.h"
#include "rpy2R_eb.h"
#include "tau_cor.h"
#include "tau_damp.h"
#include "tau_rest.h"
#include "vp.h"

extern mxArray * veh;

/*
 * The function "Mvxdot" is the implementation version of the "vxdot"
 * M-function from file "E:\RTW\NEW3\VXDOT.M" (lines 1-30). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function out=vxdot(xu)
 */
static mxArray * Mvxdot(int nargout_, mxArray * xu) {
    mxArray * out = mclGetUninitializedArray();
    mxArray * R_eb = mclGetUninitializedArray();
    mxArray * a_cee = mclGetUninitializedArray();
    mxArray * de = mclGetUninitializedArray();
    mxArray * p = mclGetUninitializedArray();
    mxArray * pdot = mclGetUninitializedArray();
    mxArray * tau_b = mclGetUninitializedArray();
    mxArray * tau_e = mclGetUninitializedArray();
    mxArray * v = mclGetUninitializedArray();
    mxArray * v_cee = mclGetUninitializedArray();
    mxArray * vc = mclGetUninitializedArray();
    mxArray * vcdot = mclGetUninitializedArray();
    mxArray * vdot = mclGetUninitializedArray();
    mclValidateInputs("vxdot", 1, &xu);
    /*
     * 
     * global veh;
     * 
     * % computes v6dof state derivatives
     * 
     * de=xu(12+[1:8]);              % fin angles
     */
    mlfAssign(
      &de,
      mlfIndexRef(
        xu,
        "(?)",
        mlfPlus(
          mlfScalar(12.0),
          mlfHorzcat(mlfColon(mlfScalar(1.0), mlfScalar(8.0), NULL), NULL))));
    /*
     * tau_b=xu(12+[9:14]);          % external force and moment wrt b
     */
    mlfAssign(
      &tau_b,
      mlfIndexRef(
        xu,
        "(?)",
        mlfPlus(
          mlfScalar(12.0),
          mlfHorzcat(mlfColon(mlfScalar(9.0), mlfScalar(14.0), NULL), NULL))));
    /*
     * tau_e=xu(12+[15:20]);         % external force and moment wrt e
     */
    mlfAssign(
      &tau_e,
      mlfIndexRef(
        xu,
        "(?)",
        mlfPlus(
          mlfScalar(12.0),
          mlfHorzcat(mlfColon(mlfScalar(15.0), mlfScalar(20.0), NULL), NULL))));
    /*
     * 
     * v_cee=xu(12+[21:23]);         % current velocity
     */
    mlfAssign(
      &v_cee,
      mlfIndexRef(
        xu,
        "(?)",
        mlfPlus(
          mlfScalar(12.0),
          mlfHorzcat(mlfColon(mlfScalar(21.0), mlfScalar(23.0), NULL), NULL))));
    /*
     * a_cee=xu(12+[24:26]);         % current acceleration
     */
    mlfAssign(
      &a_cee,
      mlfIndexRef(
        xu,
        "(?)",
        mlfPlus(
          mlfScalar(12.0),
          mlfHorzcat(mlfColon(mlfScalar(24.0), mlfScalar(26.0), NULL), NULL))));
    /*
     * 
     * p=xu(1:6);               % generalised position (eta)
     */
    mlfAssign(
      &p,
      mlfIndexRef(xu, "(?)", mlfColon(mlfScalar(1.0), mlfScalar(6.0), NULL)));
    /*
     * v=xu(7:12);              % generalised velocity (ni)
     */
    mlfAssign(
      &v,
      mlfIndexRef(xu, "(?)", mlfColon(mlfScalar(7.0), mlfScalar(12.0), NULL)));
    /*
     * 
     * % rotation matrix
     * R_eb=rpy2R_eb(p(4:6));
     */
    mlfAssign(
      &R_eb,
      mlfRpy2R_eb(
        mlfIndexRef(p, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL))));
    /*
     * 
     * % vc and vcdot
     * vc=[R_eb*v_cee; zeros(3,1)];
     */
    mlfAssign(
      &vc,
      mlfVertcat(
        mlfHorzcat(mlfMtimes(R_eb, v_cee), NULL),
        mlfHorzcat(mlfZeros(mlfScalar(3.0), mlfScalar(1.0), NULL), NULL),
        NULL));
    /*
     * vcdot=[R_eb*a_cee-vp(v(4:6),R_eb*v_cee); zeros(3,1)];
     */
    mlfAssign(
      &vcdot,
      mlfVertcat(
        mlfHorzcat(
          mlfMinus(
            mlfMtimes(R_eb, a_cee),
            mlfVp(
              mlfIndexRef(
                v, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL)),
              mlfMtimes(R_eb, v_cee))),
          NULL),
        mlfHorzcat(mlfZeros(mlfScalar(3.0), mlfScalar(1.0), NULL), NULL),
        NULL));
    /*
     * 
     * % state derivative
     * pdot=rpy2J(p(4:6))*v;
     */
    mlfAssign(
      &pdot,
      mlfMtimes(
        mlfRpy2J(
          mlfIndexRef(
            p, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL))),
        v));
    /*
     * vdot=vcdot+veh.iM*(tau_cor(veh,v,v-vc)+tau_damp(veh,v-vc,de)+...
     */
    mlfAssign(
      &vdot,
      mlfPlus(
        vcdot,
        mlfFeval(
          mclValueVarargout(),
          mlxMtimes,
          mlfIndexRef(veh, ".iM"),
          mlfPlus(
            mlfPlus(
              mlfPlus(
                mlfPlus(
                  mlfTau_cor(veh, v, mlfMinus(v, vc)),
                  mlfTau_damp(veh, mlfMinus(v, vc), de)),
                mlfTau_rest(veh, p)),
              tau_b),
            mlfVertcat(
              mlfHorzcat(
                mlfMtimes(
                  R_eb,
                  mlfIndexRef(
                    tau_e,
                    "(?)",
                    mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL))),
                NULL),
              mlfHorzcat(
                mlfMtimes(
                  R_eb,
                  mlfIndexRef(
                    tau_e,
                    "(?)",
                    mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL))),
                NULL),
              NULL)),
          NULL)));
    /*
     * tau_rest(veh,p)+tau_b+[R_eb*tau_e(1:3);R_eb*tau_e(4:6)]);
     * 
     * out=[pdot;vdot];        % final result
     */
    mlfAssign(
      &out, mlfVertcat(mlfHorzcat(pdot, NULL), mlfHorzcat(vdot, NULL), NULL));
    mclValidateOutputs("vxdot", 1, nargout_, &out);
    mxDestroyArray(R_eb);
    mxDestroyArray(a_cee);
    mxDestroyArray(de);
    mxDestroyArray(p);
    mxDestroyArray(pdot);
    mxDestroyArray(tau_b);
    mxDestroyArray(tau_e);
    mxDestroyArray(v);
    mxDestroyArray(v_cee);
    mxDestroyArray(vc);
    mxDestroyArray(vcdot);
    mxDestroyArray(vdot);
    return out;
}

/*
 * The function "mlfVxdot" contains the normal interface for the "vxdot"
 * M-function from file "E:\RTW\NEW3\VXDOT.M" (lines 1-30). This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
mxArray * mlfVxdot(mxArray * xu) {
    int nargout = 1;
    mxArray * out = mclGetUninitializedArray();
    mlfEnterNewContext(0, 1, xu);
    out = Mvxdot(nargout, xu);
    mlfRestorePreviousContext(0, 1, xu);
    return mlfReturnValue(out);
}

/*
 * The function "mlxVxdot" contains the feval interface for the "vxdot"
 * M-function from file "E:\RTW\NEW3\VXDOT.M" (lines 1-30). The feval function
 * calls the implementation version of vxdot through this function. This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
void mlxVxdot(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: vxdot Line: 1 Column: 0 The function \"vxdot"
            "\" was called with more than the declared number of outputs (1)"));
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: vxdot Line: 1 Column: 0 The function \"vxdot"
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
    mplhs[0] = Mvxdot(nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
}
