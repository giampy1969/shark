/*
 * MATLAB Compiler: 2.1
 * Date: Wed Oct 17 16:15:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-h" "-x" "-W" "mex"
 * "-L" "C" "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "vxdot.m" 
 */
#include "tau_rest.h"
#include "libmatlbm.h"
#include "rpy2R_eb.h"
#include "vp.h"

static mxChar _array1_[134] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 't', 'a', 'u', '_', 'r',
                                'e', 's', 't', ' ', 'L', 'i', 'n', 'e', ':',
                                ' ', '1', ' ', 'C', 'o', 'l', 'u', 'm', 'n',
                                ':', ' ', '1', ' ', 'T', 'h', 'e', ' ', 'f',
                                'u', 'n', 'c', 't', 'i', 'o', 'n', ' ', '"',
                                't', 'a', 'u', '_', 'r', 'e', 's', 't', '"',
                                ' ', 'w', 'a', 's', ' ', 'c', 'a', 'l', 'l',
                                'e', 'd', ' ', 'w', 'i', 't', 'h', ' ', 'm',
                                'o', 'r', 'e', ' ', 't', 'h', 'a', 'n', ' ',
                                't', 'h', 'e', ' ', 'd', 'e', 'c', 'l', 'a',
                                'r', 'e', 'd', ' ', 'n', 'u', 'm', 'b', 'e',
                                'r', ' ', 'o', 'f', ' ', 'o', 'u', 't', 'p',
                                'u', 't', 's', ' ', '(', '1', ')', '.' };
static mxArray * _mxarray0_;

static mxChar _array3_[133] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 't', 'a', 'u', '_', 'r',
                                'e', 's', 't', ' ', 'L', 'i', 'n', 'e', ':',
                                ' ', '1', ' ', 'C', 'o', 'l', 'u', 'm', 'n',
                                ':', ' ', '1', ' ', 'T', 'h', 'e', ' ', 'f',
                                'u', 'n', 'c', 't', 'i', 'o', 'n', ' ', '"',
                                't', 'a', 'u', '_', 'r', 'e', 's', 't', '"',
                                ' ', 'w', 'a', 's', ' ', 'c', 'a', 'l', 'l',
                                'e', 'd', ' ', 'w', 'i', 't', 'h', ' ', 'm',
                                'o', 'r', 'e', ' ', 't', 'h', 'a', 'n', ' ',
                                't', 'h', 'e', ' ', 'd', 'e', 'c', 'l', 'a',
                                'r', 'e', 'd', ' ', 'n', 'u', 'm', 'b', 'e',
                                'r', ' ', 'o', 'f', ' ', 'i', 'n', 'p', 'u',
                                't', 's', ' ', '(', '2', ')', '.' };
static mxArray * _mxarray2_;
static mxArray * _mxarray4_;
static mxArray * _mxarray5_;

void InitializeModule_tau_rest(void) {
    _mxarray0_ = mclInitializeString(134, _array1_);
    _mxarray2_ = mclInitializeString(133, _array3_);
    _mxarray4_ = mclInitializeDouble(4.0);
    _mxarray5_ = mclInitializeDouble(6.0);
}

void TerminateModule_tau_rest(void) {
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mtau_rest(int nargout_, mxArray * veh, mxArray * p);

_mexLocalFunctionTable _local_function_table_tau_rest
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfTau_rest" contains the normal interface for the "tau_rest"
 * M-function from file "D:\RTW\new3\tau_rest.m" (lines 1-19). This function
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
 * M-function from file "D:\RTW\new3\tau_rest.m" (lines 1-19). The feval
 * function calls the implementation version of tau_rest through this function.
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxTau_rest(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[2];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(_mxarray0_);
    }
    if (nrhs > 2) {
        mlfError(_mxarray2_);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = mclGetUninitializedArray();
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

/*
 * The function "Mtau_rest" is the implementation version of the "tau_rest"
 * M-function from file "D:\RTW\new3\tau_rest.m" (lines 1-19). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function tr=tau_rest(veh,p)
 */
static mxArray * Mtau_rest(int nargout_, mxArray * veh, mxArray * p) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_tau_rest);
    mxArray * tr = mclGetUninitializedArray();
    mxArray * tg = mclGetUninitializedArray();
    mxArray * MG_b = mclGetUninitializedArray();
    mxArray * FG_b = mclGetUninitializedArray();
    mxArray * FG_e = mclGetUninitializedArray();
    mxArray * tb = mclGetUninitializedArray();
    mxArray * MB_b = mclGetUninitializedArray();
    mxArray * FB_b = mclGetUninitializedArray();
    mxArray * FB_e = mclGetUninitializedArray();
    mclCopyArray(&veh);
    mclCopyArray(&p);
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
      mclFeval(
        mclValueVarargout(),
        mlxMtimes,
        mclFeval(
          mclValueVarargout(),
          mlxMtimes,
          mclFeval(
            mclValueVarargout(),
            mlxUminus,
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".vol")),
            NULL),
          mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".rho")),
          NULL),
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".g_e")),
        NULL));
    /*
     * FB_b=rpy2R_eb(p(4:6))*FB_e;
     */
    mlfAssign(
      &FB_b,
      mclMtimes(
        mclVe(
          mlfRpy2R_eb(
            mclVe(
              mclArrayRef1(
                mclVsa(p, "p"), mlfColon(_mxarray4_, _mxarray5_, NULL))))),
        mclVv(FB_e, "FB_e")));
    /*
     * MB_b=vp(veh.B_b,FB_b);
     */
    mlfAssign(
      &MB_b,
      mclFeval(
        mclValueVarargout(),
        mlxVp,
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".B_b")),
        mclVv(FB_b, "FB_b"),
        NULL));
    /*
     * tb=[FB_b;MB_b];
     */
    mlfAssign(&tb, mlfVertcat(mclVv(FB_b, "FB_b"), mclVv(MB_b, "MB_b"), NULL));
    /*
     * 
     * % Gravitational force and moment
     * FG_e=veh.m*veh.g_e;
     */
    mlfAssign(
      &FG_e,
      mclFeval(
        mclValueVarargout(),
        mlxMtimes,
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".m")),
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".g_e")),
        NULL));
    /*
     * FG_b=rpy2R_eb(p(4:6))*FG_e;
     */
    mlfAssign(
      &FG_b,
      mclMtimes(
        mclVe(
          mlfRpy2R_eb(
            mclVe(
              mclArrayRef1(
                mclVsa(p, "p"), mlfColon(_mxarray4_, _mxarray5_, NULL))))),
        mclVv(FG_e, "FG_e")));
    /*
     * MG_b=vp(veh.G_b,FG_b);
     */
    mlfAssign(
      &MG_b,
      mclFeval(
        mclValueVarargout(),
        mlxVp,
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".G_b")),
        mclVv(FG_b, "FG_b"),
        NULL));
    /*
     * tg=[FG_b;MG_b];
     */
    mlfAssign(&tg, mlfVertcat(mclVv(FG_b, "FG_b"), mclVv(MG_b, "MG_b"), NULL));
    /*
     * 
     * tr=tb+tg;
     */
    mlfAssign(&tr, mclPlus(mclVv(tb, "tb"), mclVv(tg, "tg")));
    mclValidateOutput(tr, 1, nargout_, "tr", "tau_rest");
    mxDestroyArray(FB_e);
    mxDestroyArray(FB_b);
    mxDestroyArray(MB_b);
    mxDestroyArray(tb);
    mxDestroyArray(FG_e);
    mxDestroyArray(FG_b);
    mxDestroyArray(MG_b);
    mxDestroyArray(tg);
    mxDestroyArray(p);
    mxDestroyArray(veh);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return tr;
}
