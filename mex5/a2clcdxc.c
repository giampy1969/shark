/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:47 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */
#include "a2clcdxc.h"

/*
 * The function "Ma2clcdxc" is the implementation version of the "a2clcdxc"
 * M-function from file "E:\RTW\NEW3\A2CLCDXC.M" (lines 1-35). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function [Cl,Cd,xcp]=a2clcdxc(alfa)
 */
static mxArray * Ma2clcdxc(mxArray * * Cd,
                           mxArray * * xcp,
                           int nargout_,
                           mxArray * alfa) {
    mxArray * Cl = mclGetUninitializedArray();
    mxArray * ALFA1 = mclGetUninitializedArray();
    mxArray * ALFA2 = mclGetUninitializedArray();
    mxArray * C1 = mclGetUninitializedArray();
    mxArray * C2 = mclGetUninitializedArray();
    mxArray * C3 = mclGetUninitializedArray();
    mxArray * C4 = mclGetUninitializedArray();
    mxArray * C5 = mclGetUninitializedArray();
    mxArray * C6 = mclGetUninitializedArray();
    mxArray * CD0 = mclGetUninitializedArray();
    mxArray * CD90 = mclGetUninitializedArray();
    mxArray * K1 = mclGetUninitializedArray();
    mxArray * K2 = mclGetUninitializedArray();
    mxArray * mod_alfa = mclGetUninitializedArray();
    mclValidateInputs("a2clcdxc", 1, &alfa);
    /*
     * 
     * % [Cl,Cd,xcp]=a2clcdxc(alfa) computes hydrodynamic coefficients Cl and Cd
     * % for the body, and distance between prow and force application point xcp
     * 
     * % Costants
     * CD0 = 0.185;
     */
    mlfAssign(&CD0, mlfScalar(0.185));
    /*
     * CD90 = 4.773;
     */
    mlfAssign(&CD90, mlfScalar(4.773));
    /*
     * ALFA1 =0.5236;
     */
    mlfAssign(&ALFA1, mlfScalar(0.5236));
    /*
     * ALFA2 =1.047;
     */
    mlfAssign(&ALFA2, mlfScalar(1.047));
    /*
     * 
     * C1 =1.64;             % interpolating coefficients
     */
    mlfAssign(&C1, mlfScalar(1.64));
    /*
     * C2 = 2.387;
     */
    mlfAssign(&C2, mlfScalar(2.387));
    /*
     * C3 =1.971;
     */
    mlfAssign(&C3, mlfScalar(1.971));
    /*
     * C4 =0.481;
     */
    mlfAssign(&C4, mlfScalar(0.481));
    /*
     * C5 =-4.861;
     */
    mlfAssign(&C5, mlfScalar(-4.861));
    /*
     * C6 =7.635;
     */
    mlfAssign(&C6, mlfScalar(7.635));
    /*
     * K1 = 0.124;
     */
    mlfAssign(&K1, mlfScalar(0.124));
    /*
     * K2 = 0.243;
     */
    mlfAssign(&K2, mlfScalar(0.243));
    /*
     * 
     * % sign correction
     * mod_alfa=abs(alfa-sign(alfa)*(abs(alfa)>pi/2)*pi);
     */
    mlfAssign(
      &mod_alfa,
      mlfAbs(
        mlfMinus(
          alfa,
          mlfMtimes(
            mlfMtimes(
              mlfSign(alfa),
              mlfGt(mlfAbs(alfa), mlfMrdivide(mlfPi(), mlfScalar(2.0)))),
            mlfPi()))));
    /*
     * 
     * Cd = CD0 + (CD90 - CD0) * sin(mod_alfa)^3;
     */
    mlfAssign(
      Cd,
      mlfPlus(
        CD0,
        mlfMtimes(
          mlfMinus(CD90, CD0), mlfMpower(mlfSin(mod_alfa), mlfScalar(3.0)))));
    /*
     * xcp = K1*mod_alfa + K2*mod_alfa^0.5;
     */
    mlfAssign(
      xcp,
      mlfPlus(
        mlfMtimes(K1, mod_alfa),
        mlfMtimes(K2, mlfMpower(mod_alfa, mlfScalar(0.5)))));
    /*
     * 
     * if mod_alfa < ALFA1    Cl = C1*mod_alfa +C2*mod_alfa^2;
     */
    if (mlfTobool(mlfLt(mod_alfa, ALFA1))) {
        mlfAssign(
          &Cl,
          mlfPlus(
            mlfMtimes(C1, mod_alfa),
            mlfMtimes(C2, mlfMpower(mod_alfa, mlfScalar(2.0)))));
    /*
     * elseif mod_alfa < ALFA2     Cl = C3*mod_alfa +C4;
     */
    } else if (mlfTobool(mlfLt(mod_alfa, ALFA2))) {
        mlfAssign(&Cl, mlfPlus(mlfMtimes(C3, mod_alfa), C4));
    /*
     * else    Cl = C5*mod_alfa +C6;
     */
    } else {
        mlfAssign(&Cl, mlfPlus(mlfMtimes(C5, mod_alfa), C6));
    /*
     * end
     */
    }
    /*
     * 
     * % sign correction
     * Cl = Cl*sign(sin(2*alfa));
     */
    mlfAssign(
      &Cl, mlfMtimes(Cl, mlfSign(mlfSin(mlfMtimes(mlfScalar(2.0), alfa)))));
    /*
     * xcp =  xcp*sign(cos(alfa)) + (abs(alfa)>pi/2);
     */
    mlfAssign(
      xcp,
      mlfPlus(
        mlfMtimes(*xcp, mlfSign(mlfCos(alfa))),
        mlfGt(mlfAbs(alfa), mlfMrdivide(mlfPi(), mlfScalar(2.0)))));
    mclValidateOutputs("a2clcdxc", 3, nargout_, &Cl, Cd, xcp);
    mxDestroyArray(ALFA1);
    mxDestroyArray(ALFA2);
    mxDestroyArray(C1);
    mxDestroyArray(C2);
    mxDestroyArray(C3);
    mxDestroyArray(C4);
    mxDestroyArray(C5);
    mxDestroyArray(C6);
    mxDestroyArray(CD0);
    mxDestroyArray(CD90);
    mxDestroyArray(K1);
    mxDestroyArray(K2);
    mxDestroyArray(mod_alfa);
    return Cl;
}

/*
 * The function "mlfA2clcdxc" contains the normal interface for the "a2clcdxc"
 * M-function from file "E:\RTW\NEW3\A2CLCDXC.M" (lines 1-35). This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
mxArray * mlfA2clcdxc(mxArray * * Cd, mxArray * * xcp, mxArray * alfa) {
    int nargout = 1;
    mxArray * Cl = mclGetUninitializedArray();
    mxArray * Cd__ = mclGetUninitializedArray();
    mxArray * xcp__ = mclGetUninitializedArray();
    mlfEnterNewContext(2, 1, Cd, xcp, alfa);
    if (Cd != NULL) {
        ++nargout;
    }
    if (xcp != NULL) {
        ++nargout;
    }
    Cl = Ma2clcdxc(&Cd__, &xcp__, nargout, alfa);
    mlfRestorePreviousContext(2, 1, Cd, xcp, alfa);
    if (Cd != NULL) {
        mclCopyOutputArg(Cd, Cd__);
    } else {
        mxDestroyArray(Cd__);
    }
    if (xcp != NULL) {
        mclCopyOutputArg(xcp, xcp__);
    } else {
        mxDestroyArray(xcp__);
    }
    return mlfReturnValue(Cl);
}

/*
 * The function "mlxA2clcdxc" contains the feval interface for the "a2clcdxc"
 * M-function from file "E:\RTW\NEW3\A2CLCDXC.M" (lines 1-35). The feval
 * function calls the implementation version of a2clcdxc through this function.
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxA2clcdxc(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[3];
    int i;
    if (nlhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: a2clcdxc Line: 1 Column:"
            " 0 The function \"a2clcdxc\" was called with m"
            "ore than the declared number of outputs (3)"));
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: a2clcdxc Line: 1 Column:"
            " 0 The function \"a2clcdxc\" was called with m"
            "ore than the declared number of inputs (1)"));
    }
    for (i = 0; i < 3; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 1 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 1; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 1, mprhs[0]);
    mplhs[0] = Ma2clcdxc(&mplhs[1], &mplhs[2], nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 3 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 3; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}
