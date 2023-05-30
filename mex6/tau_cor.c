/*
 * MATLAB Compiler: 2.1
 * Date: Wed Oct 17 16:15:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-h" "-x" "-W" "mex"
 * "-L" "C" "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "vxdot.m" 
 */
#include "tau_cor.h"
#include "libmatlbm.h"
#include "vp.h"

static mxChar _array1_[132] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 't', 'a', 'u', '_', 'c',
                                'o', 'r', ' ', 'L', 'i', 'n', 'e', ':', ' ',
                                '1', ' ', 'C', 'o', 'l', 'u', 'm', 'n', ':',
                                ' ', '1', ' ', 'T', 'h', 'e', ' ', 'f', 'u',
                                'n', 'c', 't', 'i', 'o', 'n', ' ', '"', 't',
                                'a', 'u', '_', 'c', 'o', 'r', '"', ' ', 'w',
                                'a', 's', ' ', 'c', 'a', 'l', 'l', 'e', 'd',
                                ' ', 'w', 'i', 't', 'h', ' ', 'm', 'o', 'r',
                                'e', ' ', 't', 'h', 'a', 'n', ' ', 't', 'h',
                                'e', ' ', 'd', 'e', 'c', 'l', 'a', 'r', 'e',
                                'd', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ',
                                'o', 'f', ' ', 'o', 'u', 't', 'p', 'u', 't',
                                's', ' ', '(', '1', ')', '.' };
static mxArray * _mxarray0_;

static mxChar _array3_[131] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 't', 'a', 'u', '_', 'c',
                                'o', 'r', ' ', 'L', 'i', 'n', 'e', ':', ' ',
                                '1', ' ', 'C', 'o', 'l', 'u', 'm', 'n', ':',
                                ' ', '1', ' ', 'T', 'h', 'e', ' ', 'f', 'u',
                                'n', 'c', 't', 'i', 'o', 'n', ' ', '"', 't',
                                'a', 'u', '_', 'c', 'o', 'r', '"', ' ', 'w',
                                'a', 's', ' ', 'c', 'a', 'l', 'l', 'e', 'd',
                                ' ', 'w', 'i', 't', 'h', ' ', 'm', 'o', 'r',
                                'e', ' ', 't', 'h', 'a', 'n', ' ', 't', 'h',
                                'e', ' ', 'd', 'e', 'c', 'l', 'a', 'r', 'e',
                                'd', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ',
                                'o', 'f', ' ', 'i', 'n', 'p', 'u', 't', 's',
                                ' ', '(', '3', ')', '.' };
static mxArray * _mxarray2_;
static mxArray * _mxarray4_;
static mxArray * _mxarray5_;
static mxArray * _mxarray6_;
static mxArray * _mxarray7_;

void InitializeModule_tau_cor(void) {
    _mxarray0_ = mclInitializeString(132, _array1_);
    _mxarray2_ = mclInitializeString(131, _array3_);
    _mxarray4_ = mclInitializeDouble(3.0);
    _mxarray5_ = mclInitializeDouble(1.0);
    _mxarray6_ = mclInitializeDouble(4.0);
    _mxarray7_ = mclInitializeDouble(6.0);
}

void TerminateModule_tau_cor(void) {
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mtau_cor(int nargout_,
                          mxArray * veh,
                          mxArray * v,
                          mxArray * vr);

_mexLocalFunctionTable _local_function_table_tau_cor
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfTau_cor" contains the normal interface for the "tau_cor"
 * M-function from file "D:\RTW\new3\tau_cor.m" (lines 1-14). This function
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
 * M-function from file "D:\RTW\new3\tau_cor.m" (lines 1-14). The feval
 * function calls the implementation version of tau_cor through this function.
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxTau_cor(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(_mxarray0_);
    }
    if (nrhs > 3) {
        mlfError(_mxarray2_);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = mclGetUninitializedArray();
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

/*
 * The function "Mtau_cor" is the implementation version of the "tau_cor"
 * M-function from file "D:\RTW\new3\tau_cor.m" (lines 1-14). It contains the
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
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_tau_cor);
    mxArray * tc = mclGetUninitializedArray();
    mxArray * Ca = mclGetUninitializedArray();
    mxArray * Crb = mclGetUninitializedArray();
    mclCopyArray(&veh);
    mclCopyArray(&v);
    mclCopyArray(&vr);
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
          mclVe(mlfZeros(_mxarray4_, _mxarray4_, NULL)),
          mclUminus(
            mclVe(
              mlfVp(
                mclFeval(
                  mclValueVarargout(),
                  mlxMtimes,
                  mclVe(
                    mlfIndexRef(
                      mclVsa(veh, "veh"),
                      ".Mrb(?,?)",
                      mlfColon(_mxarray5_, _mxarray4_, NULL),
                      mlfCreateColonIndex())),
                  mclVa(v, "v"),
                  NULL),
                NULL))),
          NULL),
        mlfHorzcat(
          mclUminus(
            mclVe(
              mlfVp(
                mclFeval(
                  mclValueVarargout(),
                  mlxMtimes,
                  mclVe(
                    mlfIndexRef(
                      mclVsa(veh, "veh"),
                      ".Mrb(?,?)",
                      mlfColon(_mxarray5_, _mxarray4_, NULL),
                      mlfCreateColonIndex())),
                  mclVa(v, "v"),
                  NULL),
                NULL))),
          mclUminus(
            mclVe(
              mlfVp(
                mclFeval(
                  mclValueVarargout(),
                  mlxMtimes,
                  mclVe(
                    mlfIndexRef(
                      mclVsa(veh, "veh"),
                      ".Mrb(?,?)",
                      mlfColon(_mxarray6_, _mxarray7_, NULL),
                      mlfCreateColonIndex())),
                  mclVa(v, "v"),
                  NULL),
                NULL))),
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
          mclVe(mlfZeros(_mxarray4_, _mxarray4_, NULL)),
          mclUminus(
            mclVe(
              mlfVp(
                mclFeval(
                  mclValueVarargout(),
                  mlxMtimes,
                  mclVe(
                    mlfIndexRef(
                      mclVsa(veh, "veh"),
                      ".Ma(?,?)",
                      mlfColon(_mxarray5_, _mxarray4_, NULL),
                      mlfCreateColonIndex())),
                  mclVa(vr, "vr"),
                  NULL),
                NULL))),
          NULL),
        mlfHorzcat(
          mclUminus(
            mclVe(
              mlfVp(
                mclFeval(
                  mclValueVarargout(),
                  mlxMtimes,
                  mclVe(
                    mlfIndexRef(
                      mclVsa(veh, "veh"),
                      ".Ma(?,?)",
                      mlfColon(_mxarray5_, _mxarray4_, NULL),
                      mlfCreateColonIndex())),
                  mclVa(vr, "vr"),
                  NULL),
                NULL))),
          mclUminus(
            mclVe(
              mlfVp(
                mclFeval(
                  mclValueVarargout(),
                  mlxMtimes,
                  mclVe(
                    mlfIndexRef(
                      mclVsa(veh, "veh"),
                      ".Ma(?,?)",
                      mlfColon(_mxarray6_, _mxarray7_, NULL),
                      mlfCreateColonIndex())),
                  mclVa(vr, "vr"),
                  NULL),
                NULL))),
          NULL),
        NULL));
    /*
     * -vp(veh.Ma(1:3,:)*vr), -vp(veh.Ma(4:6,:)*vr)];
     * 
     * tc=Crb*v+Ca*vr;
     */
    mlfAssign(
      &tc,
      mclPlus(
        mclMtimes(mclVv(Crb, "Crb"), mclVa(v, "v")),
        mclMtimes(mclVv(Ca, "Ca"), mclVa(vr, "vr"))));
    mclValidateOutput(tc, 1, nargout_, "tc", "tau_cor");
    mxDestroyArray(Crb);
    mxDestroyArray(Ca);
    mxDestroyArray(vr);
    mxDestroyArray(v);
    mxDestroyArray(veh);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return tc;
}
