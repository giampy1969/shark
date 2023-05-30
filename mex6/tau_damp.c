/*
 * MATLAB Compiler: 2.1
 * Date: Wed Oct 17 16:15:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-h" "-x" "-W" "mex"
 * "-L" "C" "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "vxdot.m" 
 */
#include "tau_damp.h"
#include "a2clcd.h"
#include "a2clcdxc.h"
#include "libmatlbm.h"
#include "vp.h"

static mxChar _array1_[134] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 't', 'a', 'u', '_', 'd',
                                'a', 'm', 'p', ' ', 'L', 'i', 'n', 'e', ':',
                                ' ', '1', ' ', 'C', 'o', 'l', 'u', 'm', 'n',
                                ':', ' ', '1', ' ', 'T', 'h', 'e', ' ', 'f',
                                'u', 'n', 'c', 't', 'i', 'o', 'n', ' ', '"',
                                't', 'a', 'u', '_', 'd', 'a', 'm', 'p', '"',
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
                                'l', 'e', ':', ' ', 't', 'a', 'u', '_', 'd',
                                'a', 'm', 'p', ' ', 'L', 'i', 'n', 'e', ':',
                                ' ', '1', ' ', 'C', 'o', 'l', 'u', 'm', 'n',
                                ':', ' ', '1', ' ', 'T', 'h', 'e', ' ', 'f',
                                'u', 'n', 'c', 't', 'i', 'o', 'n', ' ', '"',
                                't', 'a', 'u', '_', 'd', 'a', 'm', 'p', '"',
                                ' ', 'w', 'a', 's', ' ', 'c', 'a', 'l', 'l',
                                'e', 'd', ' ', 'w', 'i', 't', 'h', ' ', 'm',
                                'o', 'r', 'e', ' ', 't', 'h', 'a', 'n', ' ',
                                't', 'h', 'e', ' ', 'd', 'e', 'c', 'l', 'a',
                                'r', 'e', 'd', ' ', 'n', 'u', 'm', 'b', 'e',
                                'r', ' ', 'o', 'f', ' ', 'i', 'n', 'p', 'u',
                                't', 's', ' ', '(', '3', ')', '.' };
static mxArray * _mxarray2_;
static mxArray * _mxarray4_;
static mxArray * _mxarray5_;
static mxArray * _mxarray6_;
static mxArray * _mxarray7_;
static mxArray * _mxarray8_;
static mxArray * _mxarray9_;

static double _array11_[3] = { 1.0, 0.0, 0.0 };
static mxArray * _mxarray10_;
static mxArray * _mxarray12_;
static mxArray * _mxarray13_;

static double _array15_[3] = { 0.0, 0.0, 1.0 };
static mxArray * _mxarray14_;

static double _array17_[3] = { 0.0, 1.0, 0.0 };
static mxArray * _mxarray16_;
static mxArray * _mxarray18_;

static double _array20_[3] = { 0.0, -1.0, 0.0 };
static mxArray * _mxarray19_;

void InitializeModule_tau_damp(void) {
    _mxarray0_ = mclInitializeString(134, _array1_);
    _mxarray2_ = mclInitializeString(133, _array3_);
    _mxarray4_ = mclInitializeDouble(.7853981633974483);
    _mxarray5_ = mclInitializeDouble(2.0);
    _mxarray6_ = mclInitializeDouble(1.0);
    _mxarray7_ = mclInitializeDouble(3.0);
    _mxarray8_ = mclInitializeDouble(4.0);
    _mxarray9_ = mclInitializeDouble(6.0);
    _mxarray10_ = mclInitializeDoubleVector(3, 1, _array11_);
    _mxarray12_ = mclInitializeDouble(0.0);
    _mxarray13_ = mclInitializeDouble(1e-12);
    _mxarray14_ = mclInitializeDoubleVector(3, 1, _array15_);
    _mxarray16_ = mclInitializeDoubleVector(1, 3, _array17_);
    _mxarray18_ = mclInitializeDouble(-.5);
    _mxarray19_ = mclInitializeDoubleVector(1, 3, _array20_);
}

void TerminateModule_tau_damp(void) {
    mxDestroyArray(_mxarray19_);
    mxDestroyArray(_mxarray18_);
    mxDestroyArray(_mxarray16_);
    mxDestroyArray(_mxarray14_);
    mxDestroyArray(_mxarray13_);
    mxDestroyArray(_mxarray12_);
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

static mxArray * Mtau_damp(int nargout_,
                           mxArray * veh,
                           mxArray * vr,
                           mxArray * de);

_mexLocalFunctionTable _local_function_table_tau_damp
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfTau_damp" contains the normal interface for the "tau_damp"
 * M-function from file "D:\RTW\new3\tau_damp.m" (lines 1-307). This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
mxArray * mlfTau_damp(mxArray * veh, mxArray * vr, mxArray * de) {
    int nargout = 1;
    mxArray * td = mclGetUninitializedArray();
    mlfEnterNewContext(0, 3, veh, vr, de);
    td = Mtau_damp(nargout, veh, vr, de);
    mlfRestorePreviousContext(0, 3, veh, vr, de);
    return mlfReturnValue(td);
}

/*
 * The function "mlxTau_damp" contains the feval interface for the "tau_damp"
 * M-function from file "D:\RTW\new3\tau_damp.m" (lines 1-307). The feval
 * function calls the implementation version of tau_damp through this function.
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxTau_damp(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
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
    mplhs[0] = Mtau_damp(nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mtau_damp" is the implementation version of the "tau_damp"
 * M-function from file "D:\RTW\new3\tau_damp.m" (lines 1-307). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function td=tau_damp(veh,vr,de)
 */
static mxArray * Mtau_damp(int nargout_,
                           mxArray * veh,
                           mxArray * vr,
                           mxArray * de) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_tau_damp);
    mxArray * td = mclGetUninitializedArray();
    mxArray * t8 = mclGetUninitializedArray();
    mxArray * M_8bb = mclGetUninitializedArray();
    mxArray * F_8b = mclGetUninitializedArray();
    mxArray * F_8w = mclGetUninitializedArray();
    mxArray * R_w8 = mclGetUninitializedArray();
    mxArray * a8 = mclGetUninitializedArray();
    mxArray * v_8c8 = mclGetUninitializedArray();
    mxArray * v_8cb = mclGetUninitializedArray();
    mxArray * R_8b = mclGetUninitializedArray();
    mxArray * t7 = mclGetUninitializedArray();
    mxArray * M_7bb = mclGetUninitializedArray();
    mxArray * F_7b = mclGetUninitializedArray();
    mxArray * F_7w = mclGetUninitializedArray();
    mxArray * R_w7 = mclGetUninitializedArray();
    mxArray * a7 = mclGetUninitializedArray();
    mxArray * v_7c7 = mclGetUninitializedArray();
    mxArray * v_7cb = mclGetUninitializedArray();
    mxArray * R_7b = mclGetUninitializedArray();
    mxArray * t6 = mclGetUninitializedArray();
    mxArray * M_6bb = mclGetUninitializedArray();
    mxArray * F_6b = mclGetUninitializedArray();
    mxArray * F_6w = mclGetUninitializedArray();
    mxArray * R_w6 = mclGetUninitializedArray();
    mxArray * a6 = mclGetUninitializedArray();
    mxArray * v_6c6 = mclGetUninitializedArray();
    mxArray * v_6cb = mclGetUninitializedArray();
    mxArray * R_6b = mclGetUninitializedArray();
    mxArray * t5 = mclGetUninitializedArray();
    mxArray * M_5bb = mclGetUninitializedArray();
    mxArray * F_5b = mclGetUninitializedArray();
    mxArray * F_5w = mclGetUninitializedArray();
    mxArray * R_w5 = mclGetUninitializedArray();
    mxArray * a5 = mclGetUninitializedArray();
    mxArray * v_5c5 = mclGetUninitializedArray();
    mxArray * v_5cb = mclGetUninitializedArray();
    mxArray * R_5b = mclGetUninitializedArray();
    mxArray * t4 = mclGetUninitializedArray();
    mxArray * M_4bb = mclGetUninitializedArray();
    mxArray * F_4b = mclGetUninitializedArray();
    mxArray * F_4w = mclGetUninitializedArray();
    mxArray * R_w4 = mclGetUninitializedArray();
    mxArray * a4 = mclGetUninitializedArray();
    mxArray * v_4c4 = mclGetUninitializedArray();
    mxArray * v_4cb = mclGetUninitializedArray();
    mxArray * R_4b = mclGetUninitializedArray();
    mxArray * t3 = mclGetUninitializedArray();
    mxArray * M_3bb = mclGetUninitializedArray();
    mxArray * F_3b = mclGetUninitializedArray();
    mxArray * F_3w = mclGetUninitializedArray();
    mxArray * R_w3 = mclGetUninitializedArray();
    mxArray * a3 = mclGetUninitializedArray();
    mxArray * v_3c3 = mclGetUninitializedArray();
    mxArray * v_3cb = mclGetUninitializedArray();
    mxArray * R_3b = mclGetUninitializedArray();
    mxArray * t2 = mclGetUninitializedArray();
    mxArray * M_2bb = mclGetUninitializedArray();
    mxArray * F_2b = mclGetUninitializedArray();
    mxArray * F_2w = mclGetUninitializedArray();
    mxArray * R_w2 = mclGetUninitializedArray();
    mxArray * a2 = mclGetUninitializedArray();
    mxArray * v_2c2 = mclGetUninitializedArray();
    mxArray * v_2cb = mclGetUninitializedArray();
    mxArray * R_2b = mclGetUninitializedArray();
    mxArray * t1 = mclGetUninitializedArray();
    mxArray * M_1bb = mclGetUninitializedArray();
    mxArray * F_1b = mclGetUninitializedArray();
    mxArray * F_1w = mclGetUninitializedArray();
    mxArray * R_w1 = mclGetUninitializedArray();
    mxArray * a1 = mclGetUninitializedArray();
    mxArray * v_1c1 = mclGetUninitializedArray();
    mxArray * v_1cb = mclGetUninitializedArray();
    mxArray * R_1b = mclGetUninitializedArray();
    mxArray * tf = mclGetUninitializedArray();
    mxArray * M_Bbb = mclGetUninitializedArray();
    mxArray * Pf_b = mclGetUninitializedArray();
    mxArray * F_Bb = mclGetUninitializedArray();
    mxArray * F_Bw = mclGetUninitializedArray();
    mxArray * xcp = mclGetUninitializedArray();
    mxArray * cd = mclGetUninitializedArray();
    mxArray * cl = mclGetUninitializedArray();
    mxArray * R_wf = mclGetUninitializedArray();
    mxArray * af = mclGetUninitializedArray();
    mxArray * v_Bcf = mclGetUninitializedArray();
    mxArray * R_fb = mclGetUninitializedArray();
    mxArray * k_fb = mclGetUninitializedArray();
    mxArray * i_fb = mclGetUninitializedArray();
    mxArray * v_Bcb = mclGetUninitializedArray();
    mxArray * sf = mclGetUninitializedArray();
    mclCopyArray(&veh);
    mclCopyArray(&vr);
    mclCopyArray(&de);
    /*
     * 
     * % td=tau_damp(veh,vr,de); calculates damping forces from 
     * % vehicle variables ,generalized velocity vr and delta angles de
     * 
     * % --------------------------------------------------------------
     * % forces on FUSELAGE 
     * 
     * % Fuselage reference surface
     * sf=pi/4*veh.d^2;
     */
    mlfAssign(
      &sf,
      mclMtimes(
        _mxarray4_,
        mclFeval(
          mclValueVarargout(),
          mlxMpower,
          mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".d")),
          _mxarray5_,
          NULL)));
    /*
     * 
     * % relative velocity of B_b wrt c in b
     * v_Bcb=vr(1:3)+vp(vr(4:6),veh.B_b);
     */
    mlfAssign(
      &v_Bcb,
      mclPlus(
        mclVe(
          mclArrayRef1(
            mclVsa(vr, "vr"), mlfColon(_mxarray6_, _mxarray7_, NULL))),
        mclVe(
          mclFeval(
            mclValueVarargout(),
            mlxVp,
            mclVe(
              mclArrayRef1(
                mclVsa(vr, "vr"), mlfColon(_mxarray8_, _mxarray9_, NULL))),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".B_b")),
            NULL))));
    /*
     * 
     * % fuselage rotation matrix
     * i_fb=[1; 0; 0];
     */
    mlfAssign(&i_fb, _mxarray10_);
    /*
     * if   norm([0; v_Bcb(2); v_Bcb(3)])<1e-12, k_fb=[0; 0; 1];
     */
    if (mclLtBool(
          mclVe(
            mlfNorm(
              mlfVertcat(
                _mxarray12_,
                mclVe(mclIntArrayRef1(mclVsv(v_Bcb, "v_Bcb"), 2)),
                mclVe(mclIntArrayRef1(mclVsv(v_Bcb, "v_Bcb"), 3)),
                NULL),
              NULL)),
          _mxarray13_)) {
        mlfAssign(&k_fb, _mxarray14_);
    /*
     * else k_fb=[0; v_Bcb(2); v_Bcb(3)]/norm([0; v_Bcb(2); v_Bcb(3)]);end
     */
    } else {
        mlfAssign(
          &k_fb,
          mclMrdivide(
            mlfVertcat(
              _mxarray12_,
              mclVe(mclIntArrayRef1(mclVsv(v_Bcb, "v_Bcb"), 2)),
              mclVe(mclIntArrayRef1(mclVsv(v_Bcb, "v_Bcb"), 3)),
              NULL),
            mclVe(
              mlfNorm(
                mlfVertcat(
                  _mxarray12_,
                  mclVe(mclIntArrayRef1(mclVsv(v_Bcb, "v_Bcb"), 2)),
                  mclVe(mclIntArrayRef1(mclVsv(v_Bcb, "v_Bcb"), 3)),
                  NULL),
                NULL))));
    }
    /*
     * R_fb=[i_fb, vp(k_fb,i_fb), k_fb];
     */
    mlfAssign(
      &R_fb,
      mlfHorzcat(
        mclVv(i_fb, "i_fb"),
        mclVe(mlfVp(mclVv(k_fb, "k_fb"), mclVv(i_fb, "i_fb"))),
        mclVv(k_fb, "k_fb"),
        NULL));
    /*
     * 
     * % relative velocity of B_b wrt c in f
     * v_Bcf=R_fb'*v_Bcb;
     */
    mlfAssign(
      &v_Bcf,
      mclMtimes(mlfCtranspose(mclVv(R_fb, "R_fb")), mclVv(v_Bcb, "v_Bcb")));
    /*
     * 
     * % attack angle
     * af=atan2(v_Bcf(3),v_Bcf(1));
     */
    mlfAssign(
      &af,
      mlfAtan2(
        mclVe(mclIntArrayRef1(mclVsv(v_Bcf, "v_Bcf"), 3)),
        mclVe(mclIntArrayRef1(mclVsv(v_Bcf, "v_Bcf"), 1))));
    /*
     * 
     * % wind frame rotation matrix
     * R_wf=[cos(af) 0 -sin(af); 0 1 0; sin(af) 0 cos(af)];
     */
    mlfAssign(
      &R_wf,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVv(af, "af"))),
          _mxarray12_,
          mclUminus(mclVe(mlfSin(mclVv(af, "af")))),
          NULL),
        _mxarray16_,
        mlfHorzcat(
          mclVe(mlfSin(mclVv(af, "af"))),
          _mxarray12_,
          mclVe(mlfCos(mclVv(af, "af"))),
          NULL),
        NULL));
    /*
     * 
     * % cl cd xcp computation
     * [cl,cd,xcp]=a2clcdxc(af);
     */
    mlfAssign(&cl, mlfA2clcdxc(&cd, &xcp, mclVv(af, "af")));
    /*
     * 
     * % damping forces on B wrt w 
     * F_Bw=-0.5*veh.rho*sf*v_Bcf'*v_Bcf*[cd; 0; cl];
     */
    mlfAssign(
      &F_Bw,
      mclMtimes(
        mclMtimes(
          mclMtimes(
            mclMtimes(
              mclFeval(
                mclValueVarargout(),
                mlxMtimes,
                _mxarray18_,
                mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".rho")),
                NULL),
              mclVv(sf, "sf")),
            mlfCtranspose(mclVv(v_Bcf, "v_Bcf"))),
          mclVv(v_Bcf, "v_Bcf")),
        mlfVertcat(mclVv(cd, "cd"), _mxarray12_, mclVv(cl, "cl"), NULL)));
    /*
     * 
     * % damping forces on B wrt b 
     * F_Bb=R_fb*R_wf*F_Bw;
     */
    mlfAssign(
      &F_Bb,
      mclMtimes(
        mclMtimes(mclVv(R_fb, "R_fb"), mclVv(R_wf, "R_wf")),
        mclVv(F_Bw, "F_Bw")));
    /*
     * 
     * % force application point
     * Pf_b=[-xcp*veh.l; 0; 0];
     */
    mlfAssign(
      &Pf_b,
      mlfVertcat(
        mclFeval(
          mclValueVarargout(),
          mlxMtimes,
          mclUminus(mclVv(xcp, "xcp")),
          mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".l")),
          NULL),
        _mxarray12_,
        _mxarray12_,
        NULL));
    /*
     * 
     * % moments on B with pole in b wrt b
     * M_Bbb=vp(Pf_b,F_Bb);
     */
    mlfAssign(&M_Bbb, mlfVp(mclVv(Pf_b, "Pf_b"), mclVv(F_Bb, "F_Bb")));
    /*
     * 
     * tf=[F_Bb;M_Bbb];
     */
    mlfAssign(
      &tf, mlfVertcat(mclVv(F_Bb, "F_Bb"), mclVv(M_Bbb, "M_Bbb"), NULL));
    /*
     * 
     * % --------------------------------------------------------------
     * % forces on FIN 1 
     * 
     * % fin1 rotation matrix
     * R_1b=[cos(de(1)) 0 sin(de(1)); 0 1 0; -sin(de(1)) 0 cos(de(1))];
     */
    mlfAssign(
      &R_1b,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 1)))),
          _mxarray12_,
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 1)))),
          NULL),
        _mxarray16_,
        mlfHorzcat(
          mclUminus(mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 1))))),
          _mxarray12_,
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 1)))),
          NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin1 middle point wrt c in b
     * v_1cb=vr(1:3)+vp(vr(4:6),veh.P1_b);
     */
    mlfAssign(
      &v_1cb,
      mclPlus(
        mclVe(
          mclArrayRef1(
            mclVsa(vr, "vr"), mlfColon(_mxarray6_, _mxarray7_, NULL))),
        mclVe(
          mclFeval(
            mclValueVarargout(),
            mlxVp,
            mclVe(
              mclArrayRef1(
                mclVsa(vr, "vr"), mlfColon(_mxarray8_, _mxarray9_, NULL))),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P1_b")),
            NULL))));
    /*
     * 
     * % relative velocity of fin1 middle point wrt c in 1
     * v_1c1=R_1b'*v_1cb;
     */
    mlfAssign(
      &v_1c1,
      mclMtimes(mlfCtranspose(mclVv(R_1b, "R_1b")), mclVv(v_1cb, "v_1cb")));
    /*
     * 
     * % attack angle
     * a1=atan2(v_1c1(3),v_1c1(1));
     */
    mlfAssign(
      &a1,
      mlfAtan2(
        mclVe(mclIntArrayRef1(mclVsv(v_1c1, "v_1c1"), 3)),
        mclVe(mclIntArrayRef1(mclVsv(v_1c1, "v_1c1"), 1))));
    /*
     * 
     * % fin1 wind frame rotation matrix
     * R_w1=[cos(a1) 0 -sin(a1); 0 1 0; sin(a1) 0 cos(a1)];
     */
    mlfAssign(
      &R_w1,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVv(a1, "a1"))),
          _mxarray12_,
          mclUminus(mclVe(mlfSin(mclVv(a1, "a1")))),
          NULL),
        _mxarray16_,
        mlfHorzcat(
          mclVe(mlfSin(mclVv(a1, "a1"))),
          _mxarray12_,
          mclVe(mlfCos(mclVv(a1, "a1"))),
          NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a1);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, mclVv(a1, "a1")));
    /*
     * 
     * % damping forces on 1 wrt w 
     * F_1w=-0.5*veh.rho*veh.sw*(v_1c1(1)^2+v_1c1(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_1w,
      mclMtimes(
        mclMtimes(
          mclFeval(
            mclValueVarargout(),
            mlxMtimes,
            mclFeval(
              mclValueVarargout(),
              mlxMtimes,
              _mxarray18_,
              mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".rho")),
              NULL),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".sw")),
            NULL),
          mclPlus(
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_1c1, "v_1c1"), 1)), _mxarray5_),
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_1c1, "v_1c1"), 3)), _mxarray5_))),
        mlfVertcat(mclVv(cd, "cd"), _mxarray12_, mclVv(cl, "cl"), NULL)));
    /*
     * 
     * % damping forces on 1 wrt b 
     * F_1b=R_1b*R_w1*F_1w;
     */
    mlfAssign(
      &F_1b,
      mclMtimes(
        mclMtimes(mclVv(R_1b, "R_1b"), mclVv(R_w1, "R_w1")),
        mclVv(F_1w, "F_1w")));
    /*
     * 
     * % moments on 1 with pole in b wrt b
     * M_1bb=vp(veh.P1_b,F_1b);
     */
    mlfAssign(
      &M_1bb,
      mclFeval(
        mclValueVarargout(),
        mlxVp,
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P1_b")),
        mclVv(F_1b, "F_1b"),
        NULL));
    /*
     * 
     * t1=[F_1b;M_1bb];
     */
    mlfAssign(
      &t1, mlfVertcat(mclVv(F_1b, "F_1b"), mclVv(M_1bb, "M_1bb"), NULL));
    /*
     * 
     * % --------------------------------------------------------------
     * % forces on FIN 2 
     * 
     * % fin2 rotation matrix
     * R_2b=[cos(de(2)) 0 sin(de(2)); sin(de(2)) 0 -cos(de(2)); 0 1 0];
     */
    mlfAssign(
      &R_2b,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 2)))),
          _mxarray12_,
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 2)))),
          NULL),
        mlfHorzcat(
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 2)))),
          _mxarray12_,
          mclUminus(mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 2))))),
          NULL),
        _mxarray16_,
        NULL));
    /*
     * 
     * % relative velocity of fin2 middle point wrt c in b
     * v_2cb=vr(1:3)+vp(vr(4:6),veh.P2_b);
     */
    mlfAssign(
      &v_2cb,
      mclPlus(
        mclVe(
          mclArrayRef1(
            mclVsa(vr, "vr"), mlfColon(_mxarray6_, _mxarray7_, NULL))),
        mclVe(
          mclFeval(
            mclValueVarargout(),
            mlxVp,
            mclVe(
              mclArrayRef1(
                mclVsa(vr, "vr"), mlfColon(_mxarray8_, _mxarray9_, NULL))),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P2_b")),
            NULL))));
    /*
     * 
     * % relative velocity of fin2 middle point wrt c in 2
     * v_2c2=R_2b'*v_2cb;
     */
    mlfAssign(
      &v_2c2,
      mclMtimes(mlfCtranspose(mclVv(R_2b, "R_2b")), mclVv(v_2cb, "v_2cb")));
    /*
     * 
     * % attack angle
     * a2=atan2(v_2c2(3),v_2c2(1));
     */
    mlfAssign(
      &a2,
      mlfAtan2(
        mclVe(mclIntArrayRef1(mclVsv(v_2c2, "v_2c2"), 3)),
        mclVe(mclIntArrayRef1(mclVsv(v_2c2, "v_2c2"), 1))));
    /*
     * 
     * % fin2 wind frame rotation matrix
     * R_w2=[cos(a2) 0 -sin(a2); 0 1 0; sin(a2) 0 cos(a2)];
     */
    mlfAssign(
      &R_w2,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVv(a2, "a2"))),
          _mxarray12_,
          mclUminus(mclVe(mlfSin(mclVv(a2, "a2")))),
          NULL),
        _mxarray16_,
        mlfHorzcat(
          mclVe(mlfSin(mclVv(a2, "a2"))),
          _mxarray12_,
          mclVe(mlfCos(mclVv(a2, "a2"))),
          NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a2);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, mclVv(a2, "a2")));
    /*
     * 
     * % damping forces on 2 wrt w 
     * F_2w=-0.5*veh.rho*veh.sw*(v_2c2(1)^2+v_2c2(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_2w,
      mclMtimes(
        mclMtimes(
          mclFeval(
            mclValueVarargout(),
            mlxMtimes,
            mclFeval(
              mclValueVarargout(),
              mlxMtimes,
              _mxarray18_,
              mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".rho")),
              NULL),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".sw")),
            NULL),
          mclPlus(
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_2c2, "v_2c2"), 1)), _mxarray5_),
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_2c2, "v_2c2"), 3)), _mxarray5_))),
        mlfVertcat(mclVv(cd, "cd"), _mxarray12_, mclVv(cl, "cl"), NULL)));
    /*
     * 
     * % damping forces on 2 wrt b 
     * F_2b=R_2b*R_w2*F_2w;
     */
    mlfAssign(
      &F_2b,
      mclMtimes(
        mclMtimes(mclVv(R_2b, "R_2b"), mclVv(R_w2, "R_w2")),
        mclVv(F_2w, "F_2w")));
    /*
     * 
     * % moments on 2 with pole in b wrt b
     * M_2bb=vp(veh.P2_b,F_2b);
     */
    mlfAssign(
      &M_2bb,
      mclFeval(
        mclValueVarargout(),
        mlxVp,
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P2_b")),
        mclVv(F_2b, "F_2b"),
        NULL));
    /*
     * 
     * t2=[F_2b;M_2bb];
     */
    mlfAssign(
      &t2, mlfVertcat(mclVv(F_2b, "F_2b"), mclVv(M_2bb, "M_2bb"), NULL));
    /*
     * 
     * % --------------------------------------------------------------
     * % forces on FIN 3 
     * 
     * % fin3 rotation matrix
     * R_3b=[cos(de(3)) 0 sin(de(3));  0 -1 0; sin(de(3)) 0 -cos(de(3))];
     */
    mlfAssign(
      &R_3b,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 3)))),
          _mxarray12_,
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 3)))),
          NULL),
        _mxarray19_,
        mlfHorzcat(
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 3)))),
          _mxarray12_,
          mclUminus(mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 3))))),
          NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin3 middle point wrt c in b
     * v_3cb=vr(1:3)+vp(vr(4:6),veh.P3_b);
     */
    mlfAssign(
      &v_3cb,
      mclPlus(
        mclVe(
          mclArrayRef1(
            mclVsa(vr, "vr"), mlfColon(_mxarray6_, _mxarray7_, NULL))),
        mclVe(
          mclFeval(
            mclValueVarargout(),
            mlxVp,
            mclVe(
              mclArrayRef1(
                mclVsa(vr, "vr"), mlfColon(_mxarray8_, _mxarray9_, NULL))),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P3_b")),
            NULL))));
    /*
     * 
     * % relative velocity of fin3 middle point wrt c in 3
     * v_3c3=R_3b'*v_3cb;
     */
    mlfAssign(
      &v_3c3,
      mclMtimes(mlfCtranspose(mclVv(R_3b, "R_3b")), mclVv(v_3cb, "v_3cb")));
    /*
     * 
     * % attack angle
     * a3=atan2(v_3c3(3),v_3c3(1));
     */
    mlfAssign(
      &a3,
      mlfAtan2(
        mclVe(mclIntArrayRef1(mclVsv(v_3c3, "v_3c3"), 3)),
        mclVe(mclIntArrayRef1(mclVsv(v_3c3, "v_3c3"), 1))));
    /*
     * 
     * % fin3 wind frame rotation matrix
     * R_w3=[cos(a3) 0 -sin(a3); 0 1 0; sin(a3) 0 cos(a3)];
     */
    mlfAssign(
      &R_w3,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVv(a3, "a3"))),
          _mxarray12_,
          mclUminus(mclVe(mlfSin(mclVv(a3, "a3")))),
          NULL),
        _mxarray16_,
        mlfHorzcat(
          mclVe(mlfSin(mclVv(a3, "a3"))),
          _mxarray12_,
          mclVe(mlfCos(mclVv(a3, "a3"))),
          NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a3);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, mclVv(a3, "a3")));
    /*
     * 
     * % damping forces on 3 wrt w 
     * F_3w=-0.5*veh.rho*veh.sw*(v_3c3(1)^2+v_3c3(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_3w,
      mclMtimes(
        mclMtimes(
          mclFeval(
            mclValueVarargout(),
            mlxMtimes,
            mclFeval(
              mclValueVarargout(),
              mlxMtimes,
              _mxarray18_,
              mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".rho")),
              NULL),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".sw")),
            NULL),
          mclPlus(
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_3c3, "v_3c3"), 1)), _mxarray5_),
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_3c3, "v_3c3"), 3)), _mxarray5_))),
        mlfVertcat(mclVv(cd, "cd"), _mxarray12_, mclVv(cl, "cl"), NULL)));
    /*
     * 
     * % damping forces on 3 wrt b 
     * F_3b=R_3b*R_w3*F_3w;
     */
    mlfAssign(
      &F_3b,
      mclMtimes(
        mclMtimes(mclVv(R_3b, "R_3b"), mclVv(R_w3, "R_w3")),
        mclVv(F_3w, "F_3w")));
    /*
     * 
     * % moments on 3 with pole in b wrt b
     * M_3bb=vp(veh.P3_b,F_3b);
     */
    mlfAssign(
      &M_3bb,
      mclFeval(
        mclValueVarargout(),
        mlxVp,
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P3_b")),
        mclVv(F_3b, "F_3b"),
        NULL));
    /*
     * 
     * t3=[F_3b;M_3bb];
     */
    mlfAssign(
      &t3, mlfVertcat(mclVv(F_3b, "F_3b"), mclVv(M_3bb, "M_3bb"), NULL));
    /*
     * 
     * % --------------------------------------------------------------
     * % forces on FIN 4 
     * 
     * % fin4 rotation matrix
     * R_4b=[cos(de(4)) 0 sin(de(4)); -sin(de(4)) 0 cos(de(4));  0 -1 0];
     */
    mlfAssign(
      &R_4b,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 4)))),
          _mxarray12_,
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 4)))),
          NULL),
        mlfHorzcat(
          mclUminus(mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 4))))),
          _mxarray12_,
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 4)))),
          NULL),
        _mxarray19_,
        NULL));
    /*
     * 
     * % relative velocity of fin4 middle point wrt c in b
     * v_4cb=vr(1:3)+vp(vr(4:6),veh.P4_b);
     */
    mlfAssign(
      &v_4cb,
      mclPlus(
        mclVe(
          mclArrayRef1(
            mclVsa(vr, "vr"), mlfColon(_mxarray6_, _mxarray7_, NULL))),
        mclVe(
          mclFeval(
            mclValueVarargout(),
            mlxVp,
            mclVe(
              mclArrayRef1(
                mclVsa(vr, "vr"), mlfColon(_mxarray8_, _mxarray9_, NULL))),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P4_b")),
            NULL))));
    /*
     * 
     * % relative velocity of fin4 middle point wrt c in 4
     * v_4c4=R_4b'*v_4cb;
     */
    mlfAssign(
      &v_4c4,
      mclMtimes(mlfCtranspose(mclVv(R_4b, "R_4b")), mclVv(v_4cb, "v_4cb")));
    /*
     * 
     * % attack angle
     * a4=atan2(v_4c4(3),v_4c4(1));
     */
    mlfAssign(
      &a4,
      mlfAtan2(
        mclVe(mclIntArrayRef1(mclVsv(v_4c4, "v_4c4"), 3)),
        mclVe(mclIntArrayRef1(mclVsv(v_4c4, "v_4c4"), 1))));
    /*
     * 
     * % fin4 wind frame rotation matrix
     * R_w4=[cos(a4) 0 -sin(a4); 0 1 0; sin(a4) 0 cos(a4)];
     */
    mlfAssign(
      &R_w4,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVv(a4, "a4"))),
          _mxarray12_,
          mclUminus(mclVe(mlfSin(mclVv(a4, "a4")))),
          NULL),
        _mxarray16_,
        mlfHorzcat(
          mclVe(mlfSin(mclVv(a4, "a4"))),
          _mxarray12_,
          mclVe(mlfCos(mclVv(a4, "a4"))),
          NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a4);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, mclVv(a4, "a4")));
    /*
     * 
     * % damping forces on 4 wrt w 
     * F_4w=-0.5*veh.rho*veh.sw*(v_4c4(1)^2+v_4c4(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_4w,
      mclMtimes(
        mclMtimes(
          mclFeval(
            mclValueVarargout(),
            mlxMtimes,
            mclFeval(
              mclValueVarargout(),
              mlxMtimes,
              _mxarray18_,
              mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".rho")),
              NULL),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".sw")),
            NULL),
          mclPlus(
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_4c4, "v_4c4"), 1)), _mxarray5_),
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_4c4, "v_4c4"), 3)), _mxarray5_))),
        mlfVertcat(mclVv(cd, "cd"), _mxarray12_, mclVv(cl, "cl"), NULL)));
    /*
     * 
     * % damping forces on 4 wrt b 
     * F_4b=R_4b*R_w4*F_4w;
     */
    mlfAssign(
      &F_4b,
      mclMtimes(
        mclMtimes(mclVv(R_4b, "R_4b"), mclVv(R_w4, "R_w4")),
        mclVv(F_4w, "F_4w")));
    /*
     * 
     * % moments on 4 with pole in b wrt b
     * M_4bb=vp(veh.P4_b,F_4b);
     */
    mlfAssign(
      &M_4bb,
      mclFeval(
        mclValueVarargout(),
        mlxVp,
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P4_b")),
        mclVv(F_4b, "F_4b"),
        NULL));
    /*
     * 
     * t4=[F_4b;M_4bb];
     */
    mlfAssign(
      &t4, mlfVertcat(mclVv(F_4b, "F_4b"), mclVv(M_4bb, "M_4bb"), NULL));
    /*
     * 
     * % --------------------------------------------------------------
     * % forces on FIN 5 
     * 
     * % fin5 rotation matrix
     * R_5b=[cos(de(5)) 0 sin(de(5)); 0 1 0; -sin(de(5)) 0 cos(de(5))];
     */
    mlfAssign(
      &R_5b,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 5)))),
          _mxarray12_,
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 5)))),
          NULL),
        _mxarray16_,
        mlfHorzcat(
          mclUminus(mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 5))))),
          _mxarray12_,
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 5)))),
          NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin5 middle point wrt c in b
     * v_5cb=vr(1:3)+vp(vr(4:6),veh.P5_b);
     */
    mlfAssign(
      &v_5cb,
      mclPlus(
        mclVe(
          mclArrayRef1(
            mclVsa(vr, "vr"), mlfColon(_mxarray6_, _mxarray7_, NULL))),
        mclVe(
          mclFeval(
            mclValueVarargout(),
            mlxVp,
            mclVe(
              mclArrayRef1(
                mclVsa(vr, "vr"), mlfColon(_mxarray8_, _mxarray9_, NULL))),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P5_b")),
            NULL))));
    /*
     * 
     * % relative velocity of fin5 middle point wrt c in 5
     * v_5c5=R_5b'*v_5cb;
     */
    mlfAssign(
      &v_5c5,
      mclMtimes(mlfCtranspose(mclVv(R_5b, "R_5b")), mclVv(v_5cb, "v_5cb")));
    /*
     * 
     * % attack angle
     * a5=atan2(v_5c5(3),v_5c5(1));
     */
    mlfAssign(
      &a5,
      mlfAtan2(
        mclVe(mclIntArrayRef1(mclVsv(v_5c5, "v_5c5"), 3)),
        mclVe(mclIntArrayRef1(mclVsv(v_5c5, "v_5c5"), 1))));
    /*
     * 
     * % fin5 wind frame rotation matrix
     * R_w5=[cos(a5) 0 -sin(a5); 0 1 0; sin(a5) 0 cos(a5)];
     */
    mlfAssign(
      &R_w5,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVv(a5, "a5"))),
          _mxarray12_,
          mclUminus(mclVe(mlfSin(mclVv(a5, "a5")))),
          NULL),
        _mxarray16_,
        mlfHorzcat(
          mclVe(mlfSin(mclVv(a5, "a5"))),
          _mxarray12_,
          mclVe(mlfCos(mclVv(a5, "a5"))),
          NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a5);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, mclVv(a5, "a5")));
    /*
     * 
     * % damping forces on 5 wrt w 
     * F_5w=-0.5*veh.rho*veh.st*(v_5c5(1)^2+v_5c5(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_5w,
      mclMtimes(
        mclMtimes(
          mclFeval(
            mclValueVarargout(),
            mlxMtimes,
            mclFeval(
              mclValueVarargout(),
              mlxMtimes,
              _mxarray18_,
              mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".rho")),
              NULL),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".st")),
            NULL),
          mclPlus(
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_5c5, "v_5c5"), 1)), _mxarray5_),
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_5c5, "v_5c5"), 3)), _mxarray5_))),
        mlfVertcat(mclVv(cd, "cd"), _mxarray12_, mclVv(cl, "cl"), NULL)));
    /*
     * 
     * % damping forces on 5 wrt b 
     * F_5b=R_5b*R_w5*F_5w;
     */
    mlfAssign(
      &F_5b,
      mclMtimes(
        mclMtimes(mclVv(R_5b, "R_5b"), mclVv(R_w5, "R_w5")),
        mclVv(F_5w, "F_5w")));
    /*
     * 
     * % moments on 5 with pole in b wrt b
     * M_5bb=vp(veh.P5_b,F_5b);
     */
    mlfAssign(
      &M_5bb,
      mclFeval(
        mclValueVarargout(),
        mlxVp,
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P5_b")),
        mclVv(F_5b, "F_5b"),
        NULL));
    /*
     * 
     * t5=[F_5b;M_5bb];
     */
    mlfAssign(
      &t5, mlfVertcat(mclVv(F_5b, "F_5b"), mclVv(M_5bb, "M_5bb"), NULL));
    /*
     * 
     * % --------------------------------------------------------------
     * % forces on FIN 6 
     * 
     * % fin6 rotation matrix
     * R_6b=[cos(de(6)) 0 sin(de(6)); sin(de(6)) 0 -cos(de(6)); 0 1 0];
     */
    mlfAssign(
      &R_6b,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 6)))),
          _mxarray12_,
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 6)))),
          NULL),
        mlfHorzcat(
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 6)))),
          _mxarray12_,
          mclUminus(mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 6))))),
          NULL),
        _mxarray16_,
        NULL));
    /*
     * 
     * % relative velocity of fin6 middle point wrt c in b
     * v_6cb=vr(1:3)+vp(vr(4:6),veh.P6_b);
     */
    mlfAssign(
      &v_6cb,
      mclPlus(
        mclVe(
          mclArrayRef1(
            mclVsa(vr, "vr"), mlfColon(_mxarray6_, _mxarray7_, NULL))),
        mclVe(
          mclFeval(
            mclValueVarargout(),
            mlxVp,
            mclVe(
              mclArrayRef1(
                mclVsa(vr, "vr"), mlfColon(_mxarray8_, _mxarray9_, NULL))),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P6_b")),
            NULL))));
    /*
     * 
     * % relative velocity of fin6 middle point wrt c in 6
     * v_6c6=R_6b'*v_6cb;
     */
    mlfAssign(
      &v_6c6,
      mclMtimes(mlfCtranspose(mclVv(R_6b, "R_6b")), mclVv(v_6cb, "v_6cb")));
    /*
     * 
     * % attack angle
     * a6=atan2(v_6c6(3),v_6c6(1));
     */
    mlfAssign(
      &a6,
      mlfAtan2(
        mclVe(mclIntArrayRef1(mclVsv(v_6c6, "v_6c6"), 3)),
        mclVe(mclIntArrayRef1(mclVsv(v_6c6, "v_6c6"), 1))));
    /*
     * 
     * % fin6 wind frame rotation matrix
     * R_w6=[cos(a6) 0 -sin(a6); 0 1 0; sin(a6) 0 cos(a6)];
     */
    mlfAssign(
      &R_w6,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVv(a6, "a6"))),
          _mxarray12_,
          mclUminus(mclVe(mlfSin(mclVv(a6, "a6")))),
          NULL),
        _mxarray16_,
        mlfHorzcat(
          mclVe(mlfSin(mclVv(a6, "a6"))),
          _mxarray12_,
          mclVe(mlfCos(mclVv(a6, "a6"))),
          NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a6);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, mclVv(a6, "a6")));
    /*
     * 
     * % damping forces on 6 wrt w 
     * F_6w=-0.5*veh.rho*veh.st*(v_6c6(1)^2+v_6c6(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_6w,
      mclMtimes(
        mclMtimes(
          mclFeval(
            mclValueVarargout(),
            mlxMtimes,
            mclFeval(
              mclValueVarargout(),
              mlxMtimes,
              _mxarray18_,
              mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".rho")),
              NULL),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".st")),
            NULL),
          mclPlus(
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_6c6, "v_6c6"), 1)), _mxarray5_),
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_6c6, "v_6c6"), 3)), _mxarray5_))),
        mlfVertcat(mclVv(cd, "cd"), _mxarray12_, mclVv(cl, "cl"), NULL)));
    /*
     * 
     * % damping forces on 6 wrt b 
     * F_6b=R_6b*R_w6*F_6w;
     */
    mlfAssign(
      &F_6b,
      mclMtimes(
        mclMtimes(mclVv(R_6b, "R_6b"), mclVv(R_w6, "R_w6")),
        mclVv(F_6w, "F_6w")));
    /*
     * 
     * % moments on 6 with pole in b wrt b
     * M_6bb=vp(veh.P6_b,F_6b);
     */
    mlfAssign(
      &M_6bb,
      mclFeval(
        mclValueVarargout(),
        mlxVp,
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P6_b")),
        mclVv(F_6b, "F_6b"),
        NULL));
    /*
     * 
     * t6=[F_6b;M_6bb];
     */
    mlfAssign(
      &t6, mlfVertcat(mclVv(F_6b, "F_6b"), mclVv(M_6bb, "M_6bb"), NULL));
    /*
     * 
     * % --------------------------------------------------------------
     * % forces on FIN 7 
     * 
     * % fin7 rotation matrix
     * R_7b=[cos(de(7)) 0 sin(de(7));  0 -1 0; sin(de(7)) 0 -cos(de(7))];
     */
    mlfAssign(
      &R_7b,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 7)))),
          _mxarray12_,
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 7)))),
          NULL),
        _mxarray19_,
        mlfHorzcat(
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 7)))),
          _mxarray12_,
          mclUminus(mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 7))))),
          NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin7 middle point wrt c in b
     * v_7cb=vr(1:3)+vp(vr(4:6),veh.P7_b);
     */
    mlfAssign(
      &v_7cb,
      mclPlus(
        mclVe(
          mclArrayRef1(
            mclVsa(vr, "vr"), mlfColon(_mxarray6_, _mxarray7_, NULL))),
        mclVe(
          mclFeval(
            mclValueVarargout(),
            mlxVp,
            mclVe(
              mclArrayRef1(
                mclVsa(vr, "vr"), mlfColon(_mxarray8_, _mxarray9_, NULL))),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P7_b")),
            NULL))));
    /*
     * 
     * % relative velocity of fin7 middle point wrt c in 7
     * v_7c7=R_7b'*v_7cb;
     */
    mlfAssign(
      &v_7c7,
      mclMtimes(mlfCtranspose(mclVv(R_7b, "R_7b")), mclVv(v_7cb, "v_7cb")));
    /*
     * 
     * % attack angle
     * a7=atan2(v_7c7(3),v_7c7(1));
     */
    mlfAssign(
      &a7,
      mlfAtan2(
        mclVe(mclIntArrayRef1(mclVsv(v_7c7, "v_7c7"), 3)),
        mclVe(mclIntArrayRef1(mclVsv(v_7c7, "v_7c7"), 1))));
    /*
     * 
     * % fin7 wind frame rotation matrix
     * R_w7=[cos(a7) 0 -sin(a7); 0 1 0; sin(a7) 0 cos(a7)];
     */
    mlfAssign(
      &R_w7,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVv(a7, "a7"))),
          _mxarray12_,
          mclUminus(mclVe(mlfSin(mclVv(a7, "a7")))),
          NULL),
        _mxarray16_,
        mlfHorzcat(
          mclVe(mlfSin(mclVv(a7, "a7"))),
          _mxarray12_,
          mclVe(mlfCos(mclVv(a7, "a7"))),
          NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a7);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, mclVv(a7, "a7")));
    /*
     * 
     * % damping forces on 7 wrt w 
     * F_7w=-0.5*veh.rho*veh.st*(v_7c7(1)^2+v_7c7(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_7w,
      mclMtimes(
        mclMtimes(
          mclFeval(
            mclValueVarargout(),
            mlxMtimes,
            mclFeval(
              mclValueVarargout(),
              mlxMtimes,
              _mxarray18_,
              mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".rho")),
              NULL),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".st")),
            NULL),
          mclPlus(
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_7c7, "v_7c7"), 1)), _mxarray5_),
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_7c7, "v_7c7"), 3)), _mxarray5_))),
        mlfVertcat(mclVv(cd, "cd"), _mxarray12_, mclVv(cl, "cl"), NULL)));
    /*
     * 
     * % damping forces on 7 wrt b 
     * F_7b=R_7b*R_w7*F_7w;
     */
    mlfAssign(
      &F_7b,
      mclMtimes(
        mclMtimes(mclVv(R_7b, "R_7b"), mclVv(R_w7, "R_w7")),
        mclVv(F_7w, "F_7w")));
    /*
     * 
     * % moments on 7 with pole in b wrt b
     * M_7bb=vp(veh.P7_b,F_7b);
     */
    mlfAssign(
      &M_7bb,
      mclFeval(
        mclValueVarargout(),
        mlxVp,
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P7_b")),
        mclVv(F_7b, "F_7b"),
        NULL));
    /*
     * 
     * t7=[F_7b;M_7bb];
     */
    mlfAssign(
      &t7, mlfVertcat(mclVv(F_7b, "F_7b"), mclVv(M_7bb, "M_7bb"), NULL));
    /*
     * 
     * % --------------------------------------------------------------
     * % forces on FIN 8 
     * 
     * % fin8 rotation matrix
     * R_8b=[cos(de(8)) 0 sin(de(8)); -sin(de(8)) 0 cos(de(8));  0 -1 0];
     */
    mlfAssign(
      &R_8b,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 8)))),
          _mxarray12_,
          mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 8)))),
          NULL),
        mlfHorzcat(
          mclUminus(mclVe(mlfSin(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 8))))),
          _mxarray12_,
          mclVe(mlfCos(mclVe(mclIntArrayRef1(mclVsa(de, "de"), 8)))),
          NULL),
        _mxarray19_,
        NULL));
    /*
     * 
     * % relative velocity of fin8 middle point wrt c in b
     * v_8cb=vr(1:3)+vp(vr(4:6),veh.P8_b);
     */
    mlfAssign(
      &v_8cb,
      mclPlus(
        mclVe(
          mclArrayRef1(
            mclVsa(vr, "vr"), mlfColon(_mxarray6_, _mxarray7_, NULL))),
        mclVe(
          mclFeval(
            mclValueVarargout(),
            mlxVp,
            mclVe(
              mclArrayRef1(
                mclVsa(vr, "vr"), mlfColon(_mxarray8_, _mxarray9_, NULL))),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P8_b")),
            NULL))));
    /*
     * 
     * % relative velocity of fin8 middle point wrt c in 8
     * v_8c8=R_8b'*v_8cb;
     */
    mlfAssign(
      &v_8c8,
      mclMtimes(mlfCtranspose(mclVv(R_8b, "R_8b")), mclVv(v_8cb, "v_8cb")));
    /*
     * 
     * % attack angle
     * a8=atan2(v_8c8(3),v_8c8(1));
     */
    mlfAssign(
      &a8,
      mlfAtan2(
        mclVe(mclIntArrayRef1(mclVsv(v_8c8, "v_8c8"), 3)),
        mclVe(mclIntArrayRef1(mclVsv(v_8c8, "v_8c8"), 1))));
    /*
     * 
     * % fin8 wind frame rotation matrix
     * R_w8=[cos(a8) 0 -sin(a8); 0 1 0; sin(a8) 0 cos(a8)];
     */
    mlfAssign(
      &R_w8,
      mlfVertcat(
        mlfHorzcat(
          mclVe(mlfCos(mclVv(a8, "a8"))),
          _mxarray12_,
          mclUminus(mclVe(mlfSin(mclVv(a8, "a8")))),
          NULL),
        _mxarray16_,
        mlfHorzcat(
          mclVe(mlfSin(mclVv(a8, "a8"))),
          _mxarray12_,
          mclVe(mlfCos(mclVv(a8, "a8"))),
          NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a8);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, mclVv(a8, "a8")));
    /*
     * 
     * % damping forces on 8 wrt w 
     * F_8w=-0.5*veh.rho*veh.st*(v_8c8(1)^2+v_8c8(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_8w,
      mclMtimes(
        mclMtimes(
          mclFeval(
            mclValueVarargout(),
            mlxMtimes,
            mclFeval(
              mclValueVarargout(),
              mlxMtimes,
              _mxarray18_,
              mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".rho")),
              NULL),
            mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".st")),
            NULL),
          mclPlus(
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_8c8, "v_8c8"), 1)), _mxarray5_),
            mclMpower(
              mclVe(mclIntArrayRef1(mclVsv(v_8c8, "v_8c8"), 3)), _mxarray5_))),
        mlfVertcat(mclVv(cd, "cd"), _mxarray12_, mclVv(cl, "cl"), NULL)));
    /*
     * 
     * % damping forces on 8 wrt b 
     * F_8b=R_8b*R_w8*F_8w;
     */
    mlfAssign(
      &F_8b,
      mclMtimes(
        mclMtimes(mclVv(R_8b, "R_8b"), mclVv(R_w8, "R_w8")),
        mclVv(F_8w, "F_8w")));
    /*
     * 
     * % moments on 8 with pole in b wrt b
     * M_8bb=vp(veh.P8_b,F_8b);
     */
    mlfAssign(
      &M_8bb,
      mclFeval(
        mclValueVarargout(),
        mlxVp,
        mclVe(mlfIndexRef(mclVsa(veh, "veh"), ".P8_b")),
        mclVv(F_8b, "F_8b"),
        NULL));
    /*
     * 
     * t8=[F_8b;M_8bb];
     */
    mlfAssign(
      &t8, mlfVertcat(mclVv(F_8b, "F_8b"), mclVv(M_8bb, "M_8bb"), NULL));
    /*
     * 
     * % --------------------------------------------------------------
     * % resulting hydrodynamic force and moment with pole in b wrt b
     * 
     * td=tf+t1+t2+t3+t4+t5+t6+t7+t8; 
     */
    mlfAssign(
      &td,
      mclPlus(
        mclPlus(
          mclPlus(
            mclPlus(
              mclPlus(
                mclPlus(
                  mclPlus(
                    mclPlus(mclVv(tf, "tf"), mclVv(t1, "t1")), mclVv(t2, "t2")),
                  mclVv(t3, "t3")),
                mclVv(t4, "t4")),
              mclVv(t5, "t5")),
            mclVv(t6, "t6")),
          mclVv(t7, "t7")),
        mclVv(t8, "t8")));
    mclValidateOutput(td, 1, nargout_, "td", "tau_damp");
    mxDestroyArray(sf);
    mxDestroyArray(v_Bcb);
    mxDestroyArray(i_fb);
    mxDestroyArray(k_fb);
    mxDestroyArray(R_fb);
    mxDestroyArray(v_Bcf);
    mxDestroyArray(af);
    mxDestroyArray(R_wf);
    mxDestroyArray(cl);
    mxDestroyArray(cd);
    mxDestroyArray(xcp);
    mxDestroyArray(F_Bw);
    mxDestroyArray(F_Bb);
    mxDestroyArray(Pf_b);
    mxDestroyArray(M_Bbb);
    mxDestroyArray(tf);
    mxDestroyArray(R_1b);
    mxDestroyArray(v_1cb);
    mxDestroyArray(v_1c1);
    mxDestroyArray(a1);
    mxDestroyArray(R_w1);
    mxDestroyArray(F_1w);
    mxDestroyArray(F_1b);
    mxDestroyArray(M_1bb);
    mxDestroyArray(t1);
    mxDestroyArray(R_2b);
    mxDestroyArray(v_2cb);
    mxDestroyArray(v_2c2);
    mxDestroyArray(a2);
    mxDestroyArray(R_w2);
    mxDestroyArray(F_2w);
    mxDestroyArray(F_2b);
    mxDestroyArray(M_2bb);
    mxDestroyArray(t2);
    mxDestroyArray(R_3b);
    mxDestroyArray(v_3cb);
    mxDestroyArray(v_3c3);
    mxDestroyArray(a3);
    mxDestroyArray(R_w3);
    mxDestroyArray(F_3w);
    mxDestroyArray(F_3b);
    mxDestroyArray(M_3bb);
    mxDestroyArray(t3);
    mxDestroyArray(R_4b);
    mxDestroyArray(v_4cb);
    mxDestroyArray(v_4c4);
    mxDestroyArray(a4);
    mxDestroyArray(R_w4);
    mxDestroyArray(F_4w);
    mxDestroyArray(F_4b);
    mxDestroyArray(M_4bb);
    mxDestroyArray(t4);
    mxDestroyArray(R_5b);
    mxDestroyArray(v_5cb);
    mxDestroyArray(v_5c5);
    mxDestroyArray(a5);
    mxDestroyArray(R_w5);
    mxDestroyArray(F_5w);
    mxDestroyArray(F_5b);
    mxDestroyArray(M_5bb);
    mxDestroyArray(t5);
    mxDestroyArray(R_6b);
    mxDestroyArray(v_6cb);
    mxDestroyArray(v_6c6);
    mxDestroyArray(a6);
    mxDestroyArray(R_w6);
    mxDestroyArray(F_6w);
    mxDestroyArray(F_6b);
    mxDestroyArray(M_6bb);
    mxDestroyArray(t6);
    mxDestroyArray(R_7b);
    mxDestroyArray(v_7cb);
    mxDestroyArray(v_7c7);
    mxDestroyArray(a7);
    mxDestroyArray(R_w7);
    mxDestroyArray(F_7w);
    mxDestroyArray(F_7b);
    mxDestroyArray(M_7bb);
    mxDestroyArray(t7);
    mxDestroyArray(R_8b);
    mxDestroyArray(v_8cb);
    mxDestroyArray(v_8c8);
    mxDestroyArray(a8);
    mxDestroyArray(R_w8);
    mxDestroyArray(F_8w);
    mxDestroyArray(F_8b);
    mxDestroyArray(M_8bb);
    mxDestroyArray(t8);
    mxDestroyArray(de);
    mxDestroyArray(vr);
    mxDestroyArray(veh);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return td;
}
