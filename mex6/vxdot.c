/*
 * MATLAB Compiler: 2.1
 * Date: Wed Oct 17 16:15:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-h" "-x" "-W" "mex"
 * "-L" "C" "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "vxdot.m" 
 */
#include "vxdot.h"
#include "libmatlbm.h"
#include "rpy2J.h"
#include "rpy2R_eb.h"
#include "tau_cor.h"
#include "tau_damp.h"
#include "tau_rest.h"
#include "vp.h"

extern mxArray * veh;

static mxChar _array1_[128] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'v', 'x', 'd', 'o', 't',
                                ' ', 'L', 'i', 'n', 'e', ':', ' ', '1', ' ',
                                'C', 'o', 'l', 'u', 'm', 'n', ':', ' ', '1',
                                ' ', 'T', 'h', 'e', ' ', 'f', 'u', 'n', 'c',
                                't', 'i', 'o', 'n', ' ', '"', 'v', 'x', 'd',
                                'o', 't', '"', ' ', 'w', 'a', 's', ' ', 'c',
                                'a', 'l', 'l', 'e', 'd', ' ', 'w', 'i', 't',
                                'h', ' ', 'm', 'o', 'r', 'e', ' ', 't', 'h',
                                'a', 'n', ' ', 't', 'h', 'e', ' ', 'd', 'e',
                                'c', 'l', 'a', 'r', 'e', 'd', ' ', 'n', 'u',
                                'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ', 'o',
                                'u', 't', 'p', 'u', 't', 's', ' ', '(', '1',
                                ')', '.' };
static mxArray * _mxarray0_;

static mxChar _array3_[127] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'v', 'x', 'd', 'o', 't',
                                ' ', 'L', 'i', 'n', 'e', ':', ' ', '1', ' ',
                                'C', 'o', 'l', 'u', 'm', 'n', ':', ' ', '1',
                                ' ', 'T', 'h', 'e', ' ', 'f', 'u', 'n', 'c',
                                't', 'i', 'o', 'n', ' ', '"', 'v', 'x', 'd',
                                'o', 't', '"', ' ', 'w', 'a', 's', ' ', 'c',
                                'a', 'l', 'l', 'e', 'd', ' ', 'w', 'i', 't',
                                'h', ' ', 'm', 'o', 'r', 'e', ' ', 't', 'h',
                                'a', 'n', ' ', 't', 'h', 'e', ' ', 'd', 'e',
                                'c', 'l', 'a', 'r', 'e', 'd', ' ', 'n', 'u',
                                'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ', 'i',
                                'n', 'p', 'u', 't', 's', ' ', '(', '1', ')',
                                '.' };
static mxArray * _mxarray2_;
static mxArray * _mxarray4_;

static double _array6_[8] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
static mxArray * _mxarray5_;

static double _array8_[6] = { 9.0, 10.0, 11.0, 12.0, 13.0, 14.0 };
static mxArray * _mxarray7_;

static double _array10_[6] = { 15.0, 16.0, 17.0, 18.0, 19.0, 20.0 };
static mxArray * _mxarray9_;

static double _array12_[3] = { 21.0, 22.0, 23.0 };
static mxArray * _mxarray11_;

static double _array14_[3] = { 24.0, 25.0, 26.0 };
static mxArray * _mxarray13_;
static mxArray * _mxarray15_;
static mxArray * _mxarray16_;
static mxArray * _mxarray17_;
static mxArray * _mxarray18_;
static mxArray * _mxarray19_;

void InitializeModule_vxdot(void) {
    _mxarray0_ = mclInitializeString(128, _array1_);
    _mxarray2_ = mclInitializeString(127, _array3_);
    _mxarray4_ = mclInitializeDouble(12.0);
    _mxarray5_ = mclInitializeDoubleVector(1, 8, _array6_);
    _mxarray7_ = mclInitializeDoubleVector(1, 6, _array8_);
    _mxarray9_ = mclInitializeDoubleVector(1, 6, _array10_);
    _mxarray11_ = mclInitializeDoubleVector(1, 3, _array12_);
    _mxarray13_ = mclInitializeDoubleVector(1, 3, _array14_);
    _mxarray15_ = mclInitializeDouble(1.0);
    _mxarray16_ = mclInitializeDouble(6.0);
    _mxarray17_ = mclInitializeDouble(7.0);
    _mxarray18_ = mclInitializeDouble(4.0);
    _mxarray19_ = mclInitializeDouble(3.0);
}

void TerminateModule_vxdot(void) {
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray18_);
    mxDestroyArray(_mxarray17_);
    mxDestroyArray(_mxarray16_);
    mxDestroyArray(_mxarray15_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray11_);
    mxDestroyArray(_mxarray9_);
    mxDestroyArray(_mxarray7_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mvxdot(int nargout_, mxArray * xu);

_mexLocalFunctionTable _local_function_table_vxdot
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfVxdot" contains the normal interface for the "vxdot"
 * M-function from file "D:\RTW\new3\vxdot.m" (lines 1-30). This function
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
 * M-function from file "D:\RTW\new3\vxdot.m" (lines 1-30). The feval function
 * calls the implementation version of vxdot through this function. This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
void mlxVxdot(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[1];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(_mxarray0_);
    }
    if (nrhs > 1) {
        mlfError(_mxarray2_);
    }
    for (i = 0; i < 1; ++i) {
        mplhs[i] = mclGetUninitializedArray();
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

/*
 * The function "Mvxdot" is the implementation version of the "vxdot"
 * M-function from file "D:\RTW\new3\vxdot.m" (lines 1-30). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function out=vxdot(xu)
 */
static mxArray * Mvxdot(int nargout_, mxArray * xu) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_vxdot);
    mxArray * out = mclGetUninitializedArray();
    mxArray * vdot = mclGetUninitializedArray();
    mxArray * pdot = mclGetUninitializedArray();
    mxArray * vcdot = mclGetUninitializedArray();
    mxArray * vc = mclGetUninitializedArray();
    mxArray * R_eb = mclGetUninitializedArray();
    mxArray * v = mclGetUninitializedArray();
    mxArray * p = mclGetUninitializedArray();
    mxArray * a_cee = mclGetUninitializedArray();
    mxArray * v_cee = mclGetUninitializedArray();
    mxArray * tau_e = mclGetUninitializedArray();
    mxArray * tau_b = mclGetUninitializedArray();
    mxArray * de = mclGetUninitializedArray();
    mxArray * ans = mclGetUninitializedArray();
    mclCopyArray(&xu);
    /*
     * 
     * global veh;
     * 
     * % computes v6dof state derivatives
     * 
     * de=xu(12+[1:8]);              % fin angles
     */
    mlfAssign(
      &de, mclArrayRef1(mclVsa(xu, "xu"), mclPlus(_mxarray4_, _mxarray5_)));
    /*
     * tau_b=xu(12+[9:14]);          % external force and moment wrt b
     */
    mlfAssign(
      &tau_b, mclArrayRef1(mclVsa(xu, "xu"), mclPlus(_mxarray4_, _mxarray7_)));
    /*
     * tau_e=xu(12+[15:20]);         % external force and moment wrt e
     */
    mlfAssign(
      &tau_e, mclArrayRef1(mclVsa(xu, "xu"), mclPlus(_mxarray4_, _mxarray9_)));
    /*
     * 
     * v_cee=xu(12+[21:23]);         % current velocity
     */
    mlfAssign(
      &v_cee, mclArrayRef1(mclVsa(xu, "xu"), mclPlus(_mxarray4_, _mxarray11_)));
    /*
     * a_cee=xu(12+[24:26]);         % current acceleration
     */
    mlfAssign(
      &a_cee, mclArrayRef1(mclVsa(xu, "xu"), mclPlus(_mxarray4_, _mxarray13_)));
    /*
     * 
     * p=xu(1:6);               % generalised position (eta)
     */
    mlfAssign(
      &p,
      mclArrayRef1(mclVsa(xu, "xu"), mlfColon(_mxarray15_, _mxarray16_, NULL)));
    /*
     * v=xu(7:12);              % generalised velocity (ni)
     */
    mlfAssign(
      &v,
      mclArrayRef1(mclVsa(xu, "xu"), mlfColon(_mxarray17_, _mxarray4_, NULL)));
    /*
     * 
     * % rotation matrix
     * R_eb=rpy2R_eb(p(4:6));
     */
    mlfAssign(
      &R_eb,
      mlfRpy2R_eb(
        mclVe(
          mclArrayRef1(
            mclVsv(p, "p"), mlfColon(_mxarray18_, _mxarray16_, NULL)))));
    /*
     * 
     * % vc and vcdot
     * vc=[R_eb*v_cee; zeros(3,1)];
     */
    mlfAssign(
      &vc,
      mlfVertcat(
        mclMtimes(mclVv(R_eb, "R_eb"), mclVv(v_cee, "v_cee")),
        mclVe(mlfZeros(_mxarray19_, _mxarray15_, NULL)),
        NULL));
    /*
     * vcdot=[R_eb*a_cee-vp(v(4:6),R_eb*v_cee); zeros(3,1)];
     */
    mlfAssign(
      &vcdot,
      mlfVertcat(
        mclMinus(
          mclMtimes(mclVv(R_eb, "R_eb"), mclVv(a_cee, "a_cee")),
          mclVe(
            mlfVp(
              mclVe(
                mclArrayRef1(
                  mclVsv(v, "v"), mlfColon(_mxarray18_, _mxarray16_, NULL))),
              mclMtimes(mclVv(R_eb, "R_eb"), mclVv(v_cee, "v_cee"))))),
        mclVe(mlfZeros(_mxarray19_, _mxarray15_, NULL)),
        NULL));
    /*
     * 
     * % state derivative
     * pdot=rpy2J(p(4:6))*v;
     */
    mlfAssign(
      &pdot,
      mclMtimes(
        mclVe(
          mlfRpy2J(
            mclVe(
              mclArrayRef1(
                mclVsv(p, "p"), mlfColon(_mxarray18_, _mxarray16_, NULL))))),
        mclVv(v, "v")));
    /*
     * vdot=vcdot+veh.iM*(tau_cor(veh,v,v-vc)+tau_damp(veh,v-vc,de)+...
     */
    mlfAssign(
      &vdot,
      mclPlus(
        mclVv(vcdot, "vcdot"),
        mclFeval(
          mclValueVarargout(),
          mlxMtimes,
          mclVe(mlfIndexRef(mclVg(&veh, "veh"), ".iM")),
          mclPlus(
            mclPlus(
              mclPlus(
                mclPlus(
                  mclVe(
                    mlfTau_cor(
                      mclVg(&veh, "veh"),
                      mclVv(v, "v"),
                      mclMinus(mclVv(v, "v"), mclVv(vc, "vc")))),
                  mclVe(
                    mlfTau_damp(
                      mclVg(&veh, "veh"),
                      mclMinus(mclVv(v, "v"), mclVv(vc, "vc")),
                      mclVv(de, "de")))),
                mclVe(mlfTau_rest(mclVg(&veh, "veh"), mclVv(p, "p")))),
              mclVv(tau_b, "tau_b")),
            mlfVertcat(
              mclMtimes(
                mclVv(R_eb, "R_eb"),
                mclVe(
                  mclArrayRef1(
                    mclVsv(tau_e, "tau_e"),
                    mlfColon(_mxarray15_, _mxarray19_, NULL)))),
              mclMtimes(
                mclVv(R_eb, "R_eb"),
                mclVe(
                  mclArrayRef1(
                    mclVsv(tau_e, "tau_e"),
                    mlfColon(_mxarray18_, _mxarray16_, NULL)))),
              NULL)),
          NULL)));
    /*
     * tau_rest(veh,p)+tau_b+[R_eb*tau_e(1:3);R_eb*tau_e(4:6)]);
     * 
     * out=[pdot;vdot];        % final result
     */
    mlfAssign(&out, mlfVertcat(mclVv(pdot, "pdot"), mclVv(vdot, "vdot"), NULL));
    mclValidateOutput(out, 1, nargout_, "out", "vxdot");
    mxDestroyArray(ans);
    mxDestroyArray(de);
    mxDestroyArray(tau_b);
    mxDestroyArray(tau_e);
    mxDestroyArray(v_cee);
    mxDestroyArray(a_cee);
    mxDestroyArray(p);
    mxDestroyArray(v);
    mxDestroyArray(R_eb);
    mxDestroyArray(vc);
    mxDestroyArray(vcdot);
    mxDestroyArray(pdot);
    mxDestroyArray(vdot);
    mxDestroyArray(xu);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return out;
}
