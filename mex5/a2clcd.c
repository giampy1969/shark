/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:46 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */
#include "a2clcd.h"

/*
 * The function "Ma2clcd" is the implementation version of the "a2clcd"
 * M-function from file "E:\RTW\NEW3\A2CLCD.M" (lines 1-39). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function [Cl, Cd] = a2clcd(alfa)
 */
static mxArray * Ma2clcd(mxArray * * Cd, int nargout_, mxArray * alfa) {
    mxArray * Cl = mclGetUninitializedArray();
    mxArray * ALFA1 = mclGetUninitializedArray();
    mxArray * ALFA2 = mclGetUninitializedArray();
    mxArray * C1 = mclGetUninitializedArray();
    mxArray * C2 = mclGetUninitializedArray();
    mxArray * C3 = mclGetUninitializedArray();
    mxArray * C4 = mclGetUninitializedArray();
    mxArray * CDmin = mclGetUninitializedArray();
    mxArray * CLa = mclGetUninitializedArray();
    mxArray * CT = mclGetUninitializedArray();
    mxArray * K = mclGetUninitializedArray();
    mxArray * mod_alfa = mclGetUninitializedArray();
    mclValidateInputs("a2clcd", 1, &alfa);
    /*
     * 
     * % [Cl, Cd] = a2clcd(alfa) computes hydrodynamic coefficients Cl and Cd
     * % for the fins
     * 
     * % Costants
     * CLa = 2.865;               % CLalfa [rad^-1]
     */
    mlfAssign(&CLa, mlfScalar(2.865));
    /*
     * CDmin = 0.0115;            % CDmin
     */
    mlfAssign(&CDmin, mlfScalar(0.0115));
    /*
     * K = 0.1309;                % K
     */
    mlfAssign(&K, mlfScalar(0.1309));
    /*
     * ALFA1 = 0.419;             % alfa_stall [rad]
     */
    mlfAssign(&ALFA1, mlfScalar(0.419));
    /*
     * ALFA2 = 0.7854;            % alfa45 [rad]
     */
    mlfAssign(&ALFA2, mlfScalar(0.7854));
    /*
     * 
     * C1 = -1.0572;              % interpolating coefficients
     */
    mlfAssign(&C1, mlfScalar(-1.0572));
    /*
     * C2 = 1.6434;               % between zone 1 and 3 
     */
    mlfAssign(&C2, mlfScalar(1.6434));
    /*
     * C3 = 1.6759;
     */
    mlfAssign(&C3, mlfScalar(1.6759));
    /*
     * C4 = -0.5021;
     */
    mlfAssign(&C4, mlfScalar(-0.5021));
    /*
     * 
     * CT = 1.15;                 
     */
    mlfAssign(&CT, mlfScalar(1.15));
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
     * if mod_alfa < ALFA1
     */
    if (mlfTobool(mlfLt(mod_alfa, ALFA1))) {
        /*
         * % zone 1
         * Cl = CLa * mod_alfa ;
         */
        mlfAssign(&Cl, mlfMtimes(CLa, mod_alfa));
        /*
         * Cd = CDmin + K * Cl^2 ;
         */
        mlfAssign(
          Cd, mlfPlus(CDmin, mlfMtimes(K, mlfMpower(Cl, mlfScalar(2.0)))));
    /*
     * elseif  mod_alfa < ALFA2
     */
    } else if (mlfTobool(mlfLt(mod_alfa, ALFA2))) {
        /*
         * % zone2 
         * Cl = C1 * mod_alfa + C2;
         */
        mlfAssign(&Cl, mlfPlus(mlfMtimes(C1, mod_alfa), C2));
        /*
         * Cd = C3 * mod_alfa + C4;
         */
        mlfAssign(Cd, mlfPlus(mlfMtimes(C3, mod_alfa), C4));
    /*
     * else
     */
    } else {
        /*
         * % zone 3 , piastra
         * Cl = CT*cos(mod_alfa);
         */
        mlfAssign(&Cl, mlfMtimes(CT, mlfCos(mod_alfa)));
        /*
         * Cd = CT*sin(mod_alfa);
         */
        mlfAssign(Cd, mlfMtimes(CT, mlfSin(mod_alfa)));
    /*
     * end;
     */
    }
    /*
     * 
     * % sign correction
     * Cl = Cl*sign(sin(2*alfa));
     */
    mlfAssign(
      &Cl, mlfMtimes(Cl, mlfSign(mlfSin(mlfMtimes(mlfScalar(2.0), alfa)))));
    mclValidateOutputs("a2clcd", 2, nargout_, &Cl, Cd);
    mxDestroyArray(ALFA1);
    mxDestroyArray(ALFA2);
    mxDestroyArray(C1);
    mxDestroyArray(C2);
    mxDestroyArray(C3);
    mxDestroyArray(C4);
    mxDestroyArray(CDmin);
    mxDestroyArray(CLa);
    mxDestroyArray(CT);
    mxDestroyArray(K);
    mxDestroyArray(mod_alfa);
    return Cl;
}

/*
 * The function "mlfA2clcd" contains the normal interface for the "a2clcd"
 * M-function from file "E:\RTW\NEW3\A2CLCD.M" (lines 1-39). This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
mxArray * mlfA2clcd(mxArray * * Cd, mxArray * alfa) {
    int nargout = 1;
    mxArray * Cl = mclGetUninitializedArray();
    mxArray * Cd__ = mclGetUninitializedArray();
    mlfEnterNewContext(1, 1, Cd, alfa);
    if (Cd != NULL) {
        ++nargout;
    }
    Cl = Ma2clcd(&Cd__, nargout, alfa);
    mlfRestorePreviousContext(1, 1, Cd, alfa);
    if (Cd != NULL) {
        mclCopyOutputArg(Cd, Cd__);
    } else {
        mxDestroyArray(Cd__);
    }
    return mlfReturnValue(Cl);
}

/*
 * The function "mlxA2clcd" contains the feval interface for the "a2clcd"
 * M-function from file "E:\RTW\NEW3\A2CLCD.M" (lines 1-39). The feval function
 * calls the implementation version of a2clcd through this function. This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
void mlxA2clcd(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[2];
    int i;
    if (nlhs > 2) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: a2clcd Line: 1 Column: "
            "0 The function \"a2clcd\" was called with mor"
            "e than the declared number of outputs (2)"));
    }
    if (nrhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: a2clcd Line: 1 Column: "
            "0 The function \"a2clcd\" was called with mor"
            "e than the declared number of inputs (1)"));
    }
    for (i = 0; i < 2; ++i) {
        mplhs[i] = NULL;
    }
    for (i = 0; i < 1 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 1; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(0, 1, mprhs[0]);
    mplhs[0] = Ma2clcd(&mplhs[1], nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
    for (i = 1; i < 2 && i < nlhs; ++i) {
        plhs[i] = mplhs[i];
    }
    for (; i < 2; ++i) {
        mxDestroyArray(mplhs[i]);
    }
}
