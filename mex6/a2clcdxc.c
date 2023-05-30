/*
 * MATLAB Compiler: 2.1
 * Date: Wed Oct 17 16:15:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-h" "-x" "-W" "mex"
 * "-L" "C" "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "vxdot.m" 
 */
#include "a2clcdxc.h"
#include "libmatlbm.h"

static mxChar _array1_[134] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'a', '2', 'c', 'l', 'c',
                                'd', 'x', 'c', ' ', 'L', 'i', 'n', 'e', ':',
                                ' ', '1', ' ', 'C', 'o', 'l', 'u', 'm', 'n',
                                ':', ' ', '1', ' ', 'T', 'h', 'e', ' ', 'f',
                                'u', 'n', 'c', 't', 'i', 'o', 'n', ' ', '"',
                                'a', '2', 'c', 'l', 'c', 'd', 'x', 'c', '"',
                                ' ', 'w', 'a', 's', ' ', 'c', 'a', 'l', 'l',
                                'e', 'd', ' ', 'w', 'i', 't', 'h', ' ', 'm',
                                'o', 'r', 'e', ' ', 't', 'h', 'a', 'n', ' ',
                                't', 'h', 'e', ' ', 'd', 'e', 'c', 'l', 'a',
                                'r', 'e', 'd', ' ', 'n', 'u', 'm', 'b', 'e',
                                'r', ' ', 'o', 'f', ' ', 'o', 'u', 't', 'p',
                                'u', 't', 's', ' ', '(', '3', ')', '.' };
static mxArray * _mxarray0_;

static mxChar _array3_[133] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'a', '2', 'c', 'l', 'c',
                                'd', 'x', 'c', ' ', 'L', 'i', 'n', 'e', ':',
                                ' ', '1', ' ', 'C', 'o', 'l', 'u', 'm', 'n',
                                ':', ' ', '1', ' ', 'T', 'h', 'e', ' ', 'f',
                                'u', 'n', 'c', 't', 'i', 'o', 'n', ' ', '"',
                                'a', '2', 'c', 'l', 'c', 'd', 'x', 'c', '"',
                                ' ', 'w', 'a', 's', ' ', 'c', 'a', 'l', 'l',
                                'e', 'd', ' ', 'w', 'i', 't', 'h', ' ', 'm',
                                'o', 'r', 'e', ' ', 't', 'h', 'a', 'n', ' ',
                                't', 'h', 'e', ' ', 'd', 'e', 'c', 'l', 'a',
                                'r', 'e', 'd', ' ', 'n', 'u', 'm', 'b', 'e',
                                'r', ' ', 'o', 'f', ' ', 'i', 'n', 'p', 'u',
                                't', 's', ' ', '(', '1', ')', '.' };
static mxArray * _mxarray2_;
static mxArray * _mxarray4_;
static mxArray * _mxarray5_;
static mxArray * _mxarray6_;
static mxArray * _mxarray7_;
static mxArray * _mxarray8_;
static mxArray * _mxarray9_;
static mxArray * _mxarray10_;
static mxArray * _mxarray11_;
static mxArray * _mxarray12_;
static mxArray * _mxarray13_;
static mxArray * _mxarray14_;
static mxArray * _mxarray15_;
static mxArray * _mxarray16_;
static mxArray * _mxarray17_;
static mxArray * _mxarray18_;
static mxArray * _mxarray19_;
static mxArray * _mxarray20_;

void InitializeModule_a2clcdxc(void) {
    _mxarray0_ = mclInitializeString(134, _array1_);
    _mxarray2_ = mclInitializeString(133, _array3_);
    _mxarray4_ = mclInitializeDouble(.185);
    _mxarray5_ = mclInitializeDouble(4.773);
    _mxarray6_ = mclInitializeDouble(.5236);
    _mxarray7_ = mclInitializeDouble(1.047);
    _mxarray8_ = mclInitializeDouble(1.64);
    _mxarray9_ = mclInitializeDouble(2.387);
    _mxarray10_ = mclInitializeDouble(1.971);
    _mxarray11_ = mclInitializeDouble(.481);
    _mxarray12_ = mclInitializeDouble(-4.861);
    _mxarray13_ = mclInitializeDouble(7.635);
    _mxarray14_ = mclInitializeDouble(.124);
    _mxarray15_ = mclInitializeDouble(.243);
    _mxarray16_ = mclInitializeDouble(1.5707963267948966);
    _mxarray17_ = mclInitializeDouble(3.141592653589793);
    _mxarray18_ = mclInitializeDouble(3.0);
    _mxarray19_ = mclInitializeDouble(.5);
    _mxarray20_ = mclInitializeDouble(2.0);
}

void TerminateModule_a2clcdxc(void) {
    mxDestroyArray(_mxarray20_);
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray18_);
    mxDestroyArray(_mxarray17_);
    mxDestroyArray(_mxarray16_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray14_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray10_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Ma2clcdxc(mxArray * * Cd,
                           mxArray * * xcp,
                           int nargout_,
                           mxArray * alfa);

_mexLocalFunctionTable _local_function_table_a2clcdxc
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfA2clcdxc" contains the normal interface for the "a2clcdxc"
 * M-function from file "D:\RTW\new3\a2clcdxc.m" (lines 1-35). This function
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
 * M-function from file "D:\RTW\new3\a2clcdxc.m" (lines 1-35). The feval
 * function calls the implementation version of a2clcdxc through this function.
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxA2clcdxc(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[3];
    int i;
    if (nlhs > 3) {
        mlfError(_mxarray0_);
    }
    if (nrhs > 1) {
        mlfError(_mxarray2_);
    }
    for (i = 0; i < 3; ++i) {
        mplhs[i] = mclGetUninitializedArray();
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

/*
 * The function "Ma2clcdxc" is the implementation version of the "a2clcdxc"
 * M-function from file "D:\RTW\new3\a2clcdxc.m" (lines 1-35). It contains the
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
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_a2clcdxc);
    mxArray * Cl = mclGetUninitializedArray();
    mxArray * mod_alfa = mclGetUninitializedArray();
    mxArray * K2 = mclGetUninitializedArray();
    mxArray * K1 = mclGetUninitializedArray();
    mxArray * C6 = mclGetUninitializedArray();
    mxArray * C5 = mclGetUninitializedArray();
    mxArray * C4 = mclGetUninitializedArray();
    mxArray * C3 = mclGetUninitializedArray();
    mxArray * C2 = mclGetUninitializedArray();
    mxArray * C1 = mclGetUninitializedArray();
    mxArray * ALFA2 = mclGetUninitializedArray();
    mxArray * ALFA1 = mclGetUninitializedArray();
    mxArray * CD90 = mclGetUninitializedArray();
    mxArray * CD0 = mclGetUninitializedArray();
    mclCopyArray(&alfa);
    /*
     * 
     * % [Cl,Cd,xcp]=a2clcdxc(alfa) computes hydrodynamic coefficients Cl and Cd
     * % for the body, and distance between prow and force application point xcp
     * 
     * % Costants
     * CD0 = 0.185;
     */
    mlfAssign(&CD0, _mxarray4_);
    /*
     * CD90 = 4.773;
     */
    mlfAssign(&CD90, _mxarray5_);
    /*
     * ALFA1 =0.5236;
     */
    mlfAssign(&ALFA1, _mxarray6_);
    /*
     * ALFA2 =1.047;
     */
    mlfAssign(&ALFA2, _mxarray7_);
    /*
     * 
     * C1 =1.64;             % interpolating coefficients
     */
    mlfAssign(&C1, _mxarray8_);
    /*
     * C2 = 2.387;
     */
    mlfAssign(&C2, _mxarray9_);
    /*
     * C3 =1.971;
     */
    mlfAssign(&C3, _mxarray10_);
    /*
     * C4 =0.481;
     */
    mlfAssign(&C4, _mxarray11_);
    /*
     * C5 =-4.861;
     */
    mlfAssign(&C5, _mxarray12_);
    /*
     * C6 =7.635;
     */
    mlfAssign(&C6, _mxarray13_);
    /*
     * K1 = 0.124;
     */
    mlfAssign(&K1, _mxarray14_);
    /*
     * K2 = 0.243;
     */
    mlfAssign(&K2, _mxarray15_);
    /*
     * 
     * % sign correction
     * mod_alfa=abs(alfa-sign(alfa)*(abs(alfa)>pi/2)*pi);
     */
    mlfAssign(
      &mod_alfa,
      mlfAbs(
        mclMinus(
          mclVa(alfa, "alfa"),
          mclMtimes(
            mclMtimes(
              mclVe(mlfSign(mclVa(alfa, "alfa"))),
              mclGt(mclVe(mlfAbs(mclVa(alfa, "alfa"))), _mxarray16_)),
            _mxarray17_))));
    /*
     * 
     * Cd = CD0 + (CD90 - CD0) * sin(mod_alfa)^3;
     */
    mlfAssign(
      Cd,
      mclPlus(
        mclVv(CD0, "CD0"),
        mclMtimes(
          mclMinus(mclVv(CD90, "CD90"), mclVv(CD0, "CD0")),
          mclMpower(mclVe(mlfSin(mclVv(mod_alfa, "mod_alfa"))), _mxarray18_))));
    /*
     * xcp = K1*mod_alfa + K2*mod_alfa^0.5;
     */
    mlfAssign(
      xcp,
      mclPlus(
        mclMtimes(mclVv(K1, "K1"), mclVv(mod_alfa, "mod_alfa")),
        mclMtimes(
          mclVv(K2, "K2"),
          mclMpower(mclVv(mod_alfa, "mod_alfa"), _mxarray19_))));
    /*
     * 
     * if mod_alfa < ALFA1    Cl = C1*mod_alfa +C2*mod_alfa^2;
     */
    if (mclLtBool(mclVv(mod_alfa, "mod_alfa"), mclVv(ALFA1, "ALFA1"))) {
        mlfAssign(
          &Cl,
          mclPlus(
            mclMtimes(mclVv(C1, "C1"), mclVv(mod_alfa, "mod_alfa")),
            mclMtimes(
              mclVv(C2, "C2"),
              mclMpower(mclVv(mod_alfa, "mod_alfa"), _mxarray20_))));
    /*
     * elseif mod_alfa < ALFA2     Cl = C3*mod_alfa +C4;
     */
    } else if (mclLtBool(mclVv(mod_alfa, "mod_alfa"), mclVv(ALFA2, "ALFA2"))) {
        mlfAssign(
          &Cl,
          mclPlus(
            mclMtimes(mclVv(C3, "C3"), mclVv(mod_alfa, "mod_alfa")),
            mclVv(C4, "C4")));
    /*
     * else    Cl = C5*mod_alfa +C6;
     */
    } else {
        mlfAssign(
          &Cl,
          mclPlus(
            mclMtimes(mclVv(C5, "C5"), mclVv(mod_alfa, "mod_alfa")),
            mclVv(C6, "C6")));
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
      &Cl,
      mclMtimes(
        mclVv(Cl, "Cl"),
        mclVe(
          mlfSign(
            mclVe(mlfSin(mclMtimes(_mxarray20_, mclVa(alfa, "alfa"))))))));
    /*
     * xcp =  xcp*sign(cos(alfa)) + (abs(alfa)>pi/2);
     */
    mlfAssign(
      xcp,
      mclPlus(
        mclMtimes(
          mclVv(*xcp, "xcp"),
          mclVe(mlfSign(mclVe(mlfCos(mclVa(alfa, "alfa")))))),
        mclGt(mclVe(mlfAbs(mclVa(alfa, "alfa"))), _mxarray16_)));
    mclValidateOutput(Cl, 1, nargout_, "Cl", "a2clcdxc");
    mclValidateOutput(*Cd, 2, nargout_, "Cd", "a2clcdxc");
    mclValidateOutput(*xcp, 3, nargout_, "xcp", "a2clcdxc");
    mxDestroyArray(CD0);
    mxDestroyArray(CD90);
    mxDestroyArray(ALFA1);
    mxDestroyArray(ALFA2);
    mxDestroyArray(C1);
    mxDestroyArray(C2);
    mxDestroyArray(C3);
    mxDestroyArray(C4);
    mxDestroyArray(C5);
    mxDestroyArray(C6);
    mxDestroyArray(K1);
    mxDestroyArray(K2);
    mxDestroyArray(mod_alfa);
    mxDestroyArray(alfa);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return Cl;
}
