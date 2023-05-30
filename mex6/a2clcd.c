/*
 * MATLAB Compiler: 2.1
 * Date: Wed Oct 17 16:15:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-h" "-x" "-W" "mex"
 * "-L" "C" "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "vxdot.m" 
 */
#include "a2clcd.h"
#include "libmatlbm.h"

static mxChar _array1_[130] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'a', '2', 'c', 'l', 'c',
                                'd', ' ', 'L', 'i', 'n', 'e', ':', ' ', '1',
                                ' ', 'C', 'o', 'l', 'u', 'm', 'n', ':', ' ',
                                '1', ' ', 'T', 'h', 'e', ' ', 'f', 'u', 'n',
                                'c', 't', 'i', 'o', 'n', ' ', '"', 'a', '2',
                                'c', 'l', 'c', 'd', '"', ' ', 'w', 'a', 's',
                                ' ', 'c', 'a', 'l', 'l', 'e', 'd', ' ', 'w',
                                'i', 't', 'h', ' ', 'm', 'o', 'r', 'e', ' ',
                                't', 'h', 'a', 'n', ' ', 't', 'h', 'e', ' ',
                                'd', 'e', 'c', 'l', 'a', 'r', 'e', 'd', ' ',
                                'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o', 'f',
                                ' ', 'o', 'u', 't', 'p', 'u', 't', 's', ' ',
                                '(', '2', ')', '.' };
static mxArray * _mxarray0_;

static mxChar _array3_[129] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'a', '2', 'c', 'l', 'c',
                                'd', ' ', 'L', 'i', 'n', 'e', ':', ' ', '1',
                                ' ', 'C', 'o', 'l', 'u', 'm', 'n', ':', ' ',
                                '1', ' ', 'T', 'h', 'e', ' ', 'f', 'u', 'n',
                                'c', 't', 'i', 'o', 'n', ' ', '"', 'a', '2',
                                'c', 'l', 'c', 'd', '"', ' ', 'w', 'a', 's',
                                ' ', 'c', 'a', 'l', 'l', 'e', 'd', ' ', 'w',
                                'i', 't', 'h', ' ', 'm', 'o', 'r', 'e', ' ',
                                't', 'h', 'a', 'n', ' ', 't', 'h', 'e', ' ',
                                'd', 'e', 'c', 'l', 'a', 'r', 'e', 'd', ' ',
                                'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o', 'f',
                                ' ', 'i', 'n', 'p', 'u', 't', 's', ' ', '(',
                                '1', ')', '.' };
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

void InitializeModule_a2clcd(void) {
    _mxarray0_ = mclInitializeString(130, _array1_);
    _mxarray2_ = mclInitializeString(129, _array3_);
    _mxarray4_ = mclInitializeDouble(2.865);
    _mxarray5_ = mclInitializeDouble(.0115);
    _mxarray6_ = mclInitializeDouble(.1309);
    _mxarray7_ = mclInitializeDouble(.419);
    _mxarray8_ = mclInitializeDouble(.7854);
    _mxarray9_ = mclInitializeDouble(-1.0572);
    _mxarray10_ = mclInitializeDouble(1.6434);
    _mxarray11_ = mclInitializeDouble(1.6759);
    _mxarray12_ = mclInitializeDouble(-.5021);
    _mxarray13_ = mclInitializeDouble(1.15);
    _mxarray14_ = mclInitializeDouble(1.5707963267948966);
    _mxarray15_ = mclInitializeDouble(3.141592653589793);
    _mxarray16_ = mclInitializeDouble(2.0);
}

void TerminateModule_a2clcd(void) {
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

static mxArray * Ma2clcd(mxArray * * Cd, int nargout_, mxArray * alfa);

_mexLocalFunctionTable _local_function_table_a2clcd
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfA2clcd" contains the normal interface for the "a2clcd"
 * M-function from file "D:\RTW\new3\a2clcd.m" (lines 1-39). This function
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
 * M-function from file "D:\RTW\new3\a2clcd.m" (lines 1-39). The feval function
 * calls the implementation version of a2clcd through this function. This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
void mlxA2clcd(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[2];
    int i;
    if (nlhs > 2) {
        mlfError(_mxarray0_);
    }
    if (nrhs > 1) {
        mlfError(_mxarray2_);
    }
    for (i = 0; i < 2; ++i) {
        mplhs[i] = mclGetUninitializedArray();
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

/*
 * The function "Ma2clcd" is the implementation version of the "a2clcd"
 * M-function from file "D:\RTW\new3\a2clcd.m" (lines 1-39). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function [Cl, Cd] = a2clcd(alfa)
 */
static mxArray * Ma2clcd(mxArray * * Cd, int nargout_, mxArray * alfa) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_a2clcd);
    mxArray * Cl = mclGetUninitializedArray();
    mxArray * mod_alfa = mclGetUninitializedArray();
    mxArray * CT = mclGetUninitializedArray();
    mxArray * C4 = mclGetUninitializedArray();
    mxArray * C3 = mclGetUninitializedArray();
    mxArray * C2 = mclGetUninitializedArray();
    mxArray * C1 = mclGetUninitializedArray();
    mxArray * ALFA2 = mclGetUninitializedArray();
    mxArray * ALFA1 = mclGetUninitializedArray();
    mxArray * K = mclGetUninitializedArray();
    mxArray * CDmin = mclGetUninitializedArray();
    mxArray * CLa = mclGetUninitializedArray();
    mclCopyArray(&alfa);
    /*
     * 
     * % [Cl, Cd] = a2clcd(alfa) computes hydrodynamic coefficients Cl and Cd
     * % for the fins
     * 
     * % Costants
     * CLa = 2.865;               % CLalfa [rad^-1]
     */
    mlfAssign(&CLa, _mxarray4_);
    /*
     * CDmin = 0.0115;            % CDmin
     */
    mlfAssign(&CDmin, _mxarray5_);
    /*
     * K = 0.1309;                % K
     */
    mlfAssign(&K, _mxarray6_);
    /*
     * ALFA1 = 0.419;             % alfa_stall [rad]
     */
    mlfAssign(&ALFA1, _mxarray7_);
    /*
     * ALFA2 = 0.7854;            % alfa45 [rad]
     */
    mlfAssign(&ALFA2, _mxarray8_);
    /*
     * 
     * C1 = -1.0572;              % interpolating coefficients
     */
    mlfAssign(&C1, _mxarray9_);
    /*
     * C2 = 1.6434;               % between zone 1 and 3 
     */
    mlfAssign(&C2, _mxarray10_);
    /*
     * C3 = 1.6759;
     */
    mlfAssign(&C3, _mxarray11_);
    /*
     * C4 = -0.5021;
     */
    mlfAssign(&C4, _mxarray12_);
    /*
     * 
     * CT = 1.15;                 
     */
    mlfAssign(&CT, _mxarray13_);
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
              mclGt(mclVe(mlfAbs(mclVa(alfa, "alfa"))), _mxarray14_)),
            _mxarray15_))));
    /*
     * 
     * if mod_alfa < ALFA1
     */
    if (mclLtBool(mclVv(mod_alfa, "mod_alfa"), mclVv(ALFA1, "ALFA1"))) {
        /*
         * % zone 1
         * Cl = CLa * mod_alfa ;
         */
        mlfAssign(
          &Cl, mclMtimes(mclVv(CLa, "CLa"), mclVv(mod_alfa, "mod_alfa")));
        /*
         * Cd = CDmin + K * Cl^2 ;
         */
        mlfAssign(
          Cd,
          mclPlus(
            mclVv(CDmin, "CDmin"),
            mclMtimes(mclVv(K, "K"), mclMpower(mclVv(Cl, "Cl"), _mxarray16_))));
    /*
     * elseif  mod_alfa < ALFA2
     */
    } else if (mclLtBool(mclVv(mod_alfa, "mod_alfa"), mclVv(ALFA2, "ALFA2"))) {
        /*
         * % zone2 
         * Cl = C1 * mod_alfa + C2;
         */
        mlfAssign(
          &Cl,
          mclPlus(
            mclMtimes(mclVv(C1, "C1"), mclVv(mod_alfa, "mod_alfa")),
            mclVv(C2, "C2")));
        /*
         * Cd = C3 * mod_alfa + C4;
         */
        mlfAssign(
          Cd,
          mclPlus(
            mclMtimes(mclVv(C3, "C3"), mclVv(mod_alfa, "mod_alfa")),
            mclVv(C4, "C4")));
    /*
     * else
     */
    } else {
        /*
         * % zone 3 , piastra
         * Cl = CT*cos(mod_alfa);
         */
        mlfAssign(
          &Cl,
          mclMtimes(
            mclVv(CT, "CT"), mclVe(mlfCos(mclVv(mod_alfa, "mod_alfa")))));
        /*
         * Cd = CT*sin(mod_alfa);
         */
        mlfAssign(
          Cd,
          mclMtimes(
            mclVv(CT, "CT"), mclVe(mlfSin(mclVv(mod_alfa, "mod_alfa")))));
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
      &Cl,
      mclMtimes(
        mclVv(Cl, "Cl"),
        mclVe(
          mlfSign(
            mclVe(mlfSin(mclMtimes(_mxarray16_, mclVa(alfa, "alfa"))))))));
    mclValidateOutput(Cl, 1, nargout_, "Cl", "a2clcd");
    mclValidateOutput(*Cd, 2, nargout_, "Cd", "a2clcd");
    mxDestroyArray(CLa);
    mxDestroyArray(CDmin);
    mxDestroyArray(K);
    mxDestroyArray(ALFA1);
    mxDestroyArray(ALFA2);
    mxDestroyArray(C1);
    mxDestroyArray(C2);
    mxDestroyArray(C3);
    mxDestroyArray(C4);
    mxDestroyArray(CT);
    mxDestroyArray(mod_alfa);
    mxDestroyArray(alfa);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return Cl;
}
