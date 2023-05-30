/*
 * MATLAB Compiler: 2.0.1
 * Date: Wed Oct 17 15:45:44 2001
 * Arguments: "-h" "-x" "vxdot.m" 
 */
#include "tau_damp.h"
#include "a2clcd.h"
#include "a2clcdxc.h"
#include "vp.h"

static double __Array0_r[3] = { 1.0, 0.0, 0.0 };

static double __Array1_r[1] = { 0.0 };

static double __Array2_r[3] = { 0.0, 0.0, 1.0 };

static double __Array3_r[1] = { 0.0 };

static double __Array4_r[1] = { 0.0 };

static double __Array5_r[3] = { 0.0, 1.0, 0.0 };

static double __Array6_r[1] = { 0.0 };

static double __Array7_r[1] = { 0.0 };

static double __Array8_r[1] = { 0.0 };

static double __Array9_r[3] = { 0.0, 1.0, 0.0 };

static double __Array10_r[3] = { 0.0, 1.0, 0.0 };

static double __Array11_r[1] = { 0.0 };

static double __Array12_r[3] = { 0.0, 1.0, 0.0 };

static double __Array13_r[3] = { 0.0, 1.0, 0.0 };

static double __Array14_r[1] = { 0.0 };

static double __Array15_r[3] = { 0.0, -1.0, 0.0 };

static double __Array16_r[3] = { 0.0, 1.0, 0.0 };

static double __Array17_r[1] = { 0.0 };

static double __Array18_r[3] = { 0.0, -1.0, 0.0 };

static double __Array19_r[3] = { 0.0, 1.0, 0.0 };

static double __Array20_r[1] = { 0.0 };

static double __Array21_r[3] = { 0.0, 1.0, 0.0 };

static double __Array22_r[3] = { 0.0, 1.0, 0.0 };

static double __Array23_r[1] = { 0.0 };

static double __Array24_r[3] = { 0.0, 1.0, 0.0 };

static double __Array25_r[3] = { 0.0, 1.0, 0.0 };

static double __Array26_r[1] = { 0.0 };

static double __Array27_r[3] = { 0.0, -1.0, 0.0 };

static double __Array28_r[3] = { 0.0, 1.0, 0.0 };

static double __Array29_r[1] = { 0.0 };

static double __Array30_r[3] = { 0.0, -1.0, 0.0 };

static double __Array31_r[3] = { 0.0, 1.0, 0.0 };

static double __Array32_r[1] = { 0.0 };

/*
 * The function "Mtau_damp" is the implementation version of the "tau_damp"
 * M-function from file "E:\RTW\NEW3\TAU_DAMP.M" (lines 1-307). It contains the
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
    mxArray * td = mclGetUninitializedArray();
    mxArray * F_1b = mclGetUninitializedArray();
    mxArray * F_1w = mclGetUninitializedArray();
    mxArray * F_2b = mclGetUninitializedArray();
    mxArray * F_2w = mclGetUninitializedArray();
    mxArray * F_3b = mclGetUninitializedArray();
    mxArray * F_3w = mclGetUninitializedArray();
    mxArray * F_4b = mclGetUninitializedArray();
    mxArray * F_4w = mclGetUninitializedArray();
    mxArray * F_5b = mclGetUninitializedArray();
    mxArray * F_5w = mclGetUninitializedArray();
    mxArray * F_6b = mclGetUninitializedArray();
    mxArray * F_6w = mclGetUninitializedArray();
    mxArray * F_7b = mclGetUninitializedArray();
    mxArray * F_7w = mclGetUninitializedArray();
    mxArray * F_8b = mclGetUninitializedArray();
    mxArray * F_8w = mclGetUninitializedArray();
    mxArray * F_Bb = mclGetUninitializedArray();
    mxArray * F_Bw = mclGetUninitializedArray();
    mxArray * M_1bb = mclGetUninitializedArray();
    mxArray * M_2bb = mclGetUninitializedArray();
    mxArray * M_3bb = mclGetUninitializedArray();
    mxArray * M_4bb = mclGetUninitializedArray();
    mxArray * M_5bb = mclGetUninitializedArray();
    mxArray * M_6bb = mclGetUninitializedArray();
    mxArray * M_7bb = mclGetUninitializedArray();
    mxArray * M_8bb = mclGetUninitializedArray();
    mxArray * M_Bbb = mclGetUninitializedArray();
    mxArray * Pf_b = mclGetUninitializedArray();
    mxArray * R_1b = mclGetUninitializedArray();
    mxArray * R_2b = mclGetUninitializedArray();
    mxArray * R_3b = mclGetUninitializedArray();
    mxArray * R_4b = mclGetUninitializedArray();
    mxArray * R_5b = mclGetUninitializedArray();
    mxArray * R_6b = mclGetUninitializedArray();
    mxArray * R_7b = mclGetUninitializedArray();
    mxArray * R_8b = mclGetUninitializedArray();
    mxArray * R_fb = mclGetUninitializedArray();
    mxArray * R_w1 = mclGetUninitializedArray();
    mxArray * R_w2 = mclGetUninitializedArray();
    mxArray * R_w3 = mclGetUninitializedArray();
    mxArray * R_w4 = mclGetUninitializedArray();
    mxArray * R_w5 = mclGetUninitializedArray();
    mxArray * R_w6 = mclGetUninitializedArray();
    mxArray * R_w7 = mclGetUninitializedArray();
    mxArray * R_w8 = mclGetUninitializedArray();
    mxArray * R_wf = mclGetUninitializedArray();
    mxArray * a1 = mclGetUninitializedArray();
    mxArray * a2 = mclGetUninitializedArray();
    mxArray * a3 = mclGetUninitializedArray();
    mxArray * a4 = mclGetUninitializedArray();
    mxArray * a5 = mclGetUninitializedArray();
    mxArray * a6 = mclGetUninitializedArray();
    mxArray * a7 = mclGetUninitializedArray();
    mxArray * a8 = mclGetUninitializedArray();
    mxArray * af = mclGetUninitializedArray();
    mxArray * cd = mclGetUninitializedArray();
    mxArray * cl = mclGetUninitializedArray();
    mxArray * i_fb = mclGetUninitializedArray();
    mxArray * k_fb = mclGetUninitializedArray();
    mxArray * sf = mclGetUninitializedArray();
    mxArray * t1 = mclGetUninitializedArray();
    mxArray * t2 = mclGetUninitializedArray();
    mxArray * t3 = mclGetUninitializedArray();
    mxArray * t4 = mclGetUninitializedArray();
    mxArray * t5 = mclGetUninitializedArray();
    mxArray * t6 = mclGetUninitializedArray();
    mxArray * t7 = mclGetUninitializedArray();
    mxArray * t8 = mclGetUninitializedArray();
    mxArray * tf = mclGetUninitializedArray();
    mxArray * v_1c1 = mclGetUninitializedArray();
    mxArray * v_1cb = mclGetUninitializedArray();
    mxArray * v_2c2 = mclGetUninitializedArray();
    mxArray * v_2cb = mclGetUninitializedArray();
    mxArray * v_3c3 = mclGetUninitializedArray();
    mxArray * v_3cb = mclGetUninitializedArray();
    mxArray * v_4c4 = mclGetUninitializedArray();
    mxArray * v_4cb = mclGetUninitializedArray();
    mxArray * v_5c5 = mclGetUninitializedArray();
    mxArray * v_5cb = mclGetUninitializedArray();
    mxArray * v_6c6 = mclGetUninitializedArray();
    mxArray * v_6cb = mclGetUninitializedArray();
    mxArray * v_7c7 = mclGetUninitializedArray();
    mxArray * v_7cb = mclGetUninitializedArray();
    mxArray * v_8c8 = mclGetUninitializedArray();
    mxArray * v_8cb = mclGetUninitializedArray();
    mxArray * v_Bcb = mclGetUninitializedArray();
    mxArray * v_Bcf = mclGetUninitializedArray();
    mxArray * xcp = mclGetUninitializedArray();
    mclValidateInputs("tau_damp", 3, &veh, &vr, &de);
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
      mlfMtimes(
        mlfMrdivide(mlfPi(), mlfScalar(4.0)),
        mlfFeval(
          mclValueVarargout(),
          mlxMpower,
          mlfIndexRef(veh, ".d"),
          mlfScalar(2.0),
          NULL)));
    /*
     * 
     * % relative velocity of B_b wrt c in b
     * v_Bcb=vr(1:3)+vp(vr(4:6),veh.B_b);
     */
    mlfAssign(
      &v_Bcb,
      mlfPlus(
        mlfIndexRef(vr, "(?)", mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL)),
        mlfFeval(
          mclValueVarargout(),
          mlxVp,
          mlfIndexRef(
            vr, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL)),
          mlfIndexRef(veh, ".B_b"),
          NULL)));
    /*
     * 
     * % fuselage rotation matrix
     * i_fb=[1; 0; 0];
     */
    mlfAssign(&i_fb, mlfDoubleMatrix(3, 1, __Array0_r, NULL));
    /*
     * if   norm([0; v_Bcb(2); v_Bcb(3)])<1e-12, k_fb=[0; 0; 1];
     */
    if (mlfTobool(
          mlfLt(
            mlfNorm(
              mlfVertcat(
                mlfDoubleMatrix(1, 1, __Array1_r, NULL),
                mlfHorzcat(mlfIndexRef(v_Bcb, "(?)", mlfScalar(2.0)), NULL),
                mlfHorzcat(mlfIndexRef(v_Bcb, "(?)", mlfScalar(3.0)), NULL),
                NULL),
              NULL),
            mlfScalar(1e-12)))) {
        mlfAssign(&k_fb, mlfDoubleMatrix(3, 1, __Array2_r, NULL));
    /*
     * else k_fb=[0; v_Bcb(2); v_Bcb(3)]/norm([0; v_Bcb(2); v_Bcb(3)]);end
     */
    } else {
        mlfAssign(
          &k_fb,
          mlfMrdivide(
            mlfVertcat(
              mlfDoubleMatrix(1, 1, __Array3_r, NULL),
              mlfHorzcat(mlfIndexRef(v_Bcb, "(?)", mlfScalar(2.0)), NULL),
              mlfHorzcat(mlfIndexRef(v_Bcb, "(?)", mlfScalar(3.0)), NULL),
              NULL),
            mlfNorm(
              mlfVertcat(
                mlfDoubleMatrix(1, 1, __Array4_r, NULL),
                mlfHorzcat(mlfIndexRef(v_Bcb, "(?)", mlfScalar(2.0)), NULL),
                mlfHorzcat(mlfIndexRef(v_Bcb, "(?)", mlfScalar(3.0)), NULL),
                NULL),
              NULL)));
    }
    /*
     * R_fb=[i_fb, vp(k_fb,i_fb), k_fb];
     */
    mlfAssign(&R_fb, mlfHorzcat(i_fb, mlfVp(k_fb, i_fb), k_fb, NULL));
    /*
     * 
     * % relative velocity of B_b wrt c in f
     * v_Bcf=R_fb'*v_Bcb;
     */
    mlfAssign(&v_Bcf, mlfMtimes(mlfCtranspose(R_fb), v_Bcb));
    /*
     * 
     * % attack angle
     * af=atan2(v_Bcf(3),v_Bcf(1));
     */
    mlfAssign(
      &af,
      mlfAtan2(
        mlfIndexRef(v_Bcf, "(?)", mlfScalar(3.0)),
        mlfIndexRef(v_Bcf, "(?)", mlfScalar(1.0))));
    /*
     * 
     * % wind frame rotation matrix
     * R_wf=[cos(af) 0 -sin(af); 0 1 0; sin(af) 0 cos(af)];
     */
    mlfAssign(
      &R_wf,
      mlfVertcat(
        mlfHorzcat(mlfCos(af), mlfScalar(0.0), mlfUminus(mlfSin(af)), NULL),
        mlfDoubleMatrix(1, 3, __Array5_r, NULL),
        mlfHorzcat(mlfSin(af), mlfScalar(0.0), mlfCos(af), NULL),
        NULL));
    /*
     * 
     * % cl cd xcp computation
     * [cl,cd,xcp]=a2clcdxc(af);
     */
    mlfAssign(&cl, mlfA2clcdxc(&cd, &xcp, af));
    /*
     * 
     * % damping forces on B wrt w 
     * F_Bw=-0.5*veh.rho*sf*v_Bcf'*v_Bcf*[cd; 0; cl];
     */
    mlfAssign(
      &F_Bw,
      mlfMtimes(
        mlfMtimes(
          mlfMtimes(
            mlfMtimes(
              mlfFeval(
                mclValueVarargout(),
                mlxMtimes,
                mlfScalar(-0.5),
                mlfIndexRef(veh, ".rho"),
                NULL),
              sf),
            mlfCtranspose(v_Bcf)),
          v_Bcf),
        mlfVertcat(
          mlfHorzcat(cd, NULL),
          mlfDoubleMatrix(1, 1, __Array6_r, NULL),
          mlfHorzcat(cl, NULL),
          NULL)));
    /*
     * 
     * % damping forces on B wrt b 
     * F_Bb=R_fb*R_wf*F_Bw;
     */
    mlfAssign(&F_Bb, mlfMtimes(mlfMtimes(R_fb, R_wf), F_Bw));
    /*
     * 
     * % force application point
     * Pf_b=[-xcp*veh.l; 0; 0];
     */
    mlfAssign(
      &Pf_b,
      mlfVertcat(
        mlfHorzcat(
          mlfFeval(
            mclValueVarargout(),
            mlxMtimes,
            mlfUminus(xcp),
            mlfIndexRef(veh, ".l"),
            NULL),
          NULL),
        mlfDoubleMatrix(1, 1, __Array7_r, NULL),
        mlfDoubleMatrix(1, 1, __Array8_r, NULL),
        NULL));
    /*
     * 
     * % moments on B with pole in b wrt b
     * M_Bbb=vp(Pf_b,F_Bb);
     */
    mlfAssign(&M_Bbb, mlfVp(Pf_b, F_Bb));
    /*
     * 
     * tf=[F_Bb;M_Bbb];
     */
    mlfAssign(
      &tf, mlfVertcat(mlfHorzcat(F_Bb, NULL), mlfHorzcat(M_Bbb, NULL), NULL));
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
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(1.0))),
          mlfScalar(0.0),
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(1.0))),
          NULL),
        mlfDoubleMatrix(1, 3, __Array9_r, NULL),
        mlfHorzcat(
          mlfUminus(mlfSin(mlfIndexRef(de, "(?)", mlfScalar(1.0)))),
          mlfScalar(0.0),
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(1.0))),
          NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin1 middle point wrt c in b
     * v_1cb=vr(1:3)+vp(vr(4:6),veh.P1_b);
     */
    mlfAssign(
      &v_1cb,
      mlfPlus(
        mlfIndexRef(vr, "(?)", mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL)),
        mlfFeval(
          mclValueVarargout(),
          mlxVp,
          mlfIndexRef(
            vr, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL)),
          mlfIndexRef(veh, ".P1_b"),
          NULL)));
    /*
     * 
     * % relative velocity of fin1 middle point wrt c in 1
     * v_1c1=R_1b'*v_1cb;
     */
    mlfAssign(&v_1c1, mlfMtimes(mlfCtranspose(R_1b), v_1cb));
    /*
     * 
     * % attack angle
     * a1=atan2(v_1c1(3),v_1c1(1));
     */
    mlfAssign(
      &a1,
      mlfAtan2(
        mlfIndexRef(v_1c1, "(?)", mlfScalar(3.0)),
        mlfIndexRef(v_1c1, "(?)", mlfScalar(1.0))));
    /*
     * 
     * % fin1 wind frame rotation matrix
     * R_w1=[cos(a1) 0 -sin(a1); 0 1 0; sin(a1) 0 cos(a1)];
     */
    mlfAssign(
      &R_w1,
      mlfVertcat(
        mlfHorzcat(mlfCos(a1), mlfScalar(0.0), mlfUminus(mlfSin(a1)), NULL),
        mlfDoubleMatrix(1, 3, __Array10_r, NULL),
        mlfHorzcat(mlfSin(a1), mlfScalar(0.0), mlfCos(a1), NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a1);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, a1));
    /*
     * 
     * % damping forces on 1 wrt w 
     * F_1w=-0.5*veh.rho*veh.sw*(v_1c1(1)^2+v_1c1(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_1w,
      mlfMtimes(
        mlfMtimes(
          mlfFeval(
            mclValueVarargout(),
            mlxMtimes,
            mlfFeval(
              mclValueVarargout(),
              mlxMtimes,
              mlfScalar(-0.5),
              mlfIndexRef(veh, ".rho"),
              NULL),
            mlfIndexRef(veh, ".sw"),
            NULL),
          mlfPlus(
            mlfMpower(
              mlfIndexRef(v_1c1, "(?)", mlfScalar(1.0)), mlfScalar(2.0)),
            mlfMpower(
              mlfIndexRef(v_1c1, "(?)", mlfScalar(3.0)), mlfScalar(2.0)))),
        mlfVertcat(
          mlfHorzcat(cd, NULL),
          mlfDoubleMatrix(1, 1, __Array11_r, NULL),
          mlfHorzcat(cl, NULL),
          NULL)));
    /*
     * 
     * % damping forces on 1 wrt b 
     * F_1b=R_1b*R_w1*F_1w;
     */
    mlfAssign(&F_1b, mlfMtimes(mlfMtimes(R_1b, R_w1), F_1w));
    /*
     * 
     * % moments on 1 with pole in b wrt b
     * M_1bb=vp(veh.P1_b,F_1b);
     */
    mlfAssign(
      &M_1bb,
      mlfFeval(
        mclValueVarargout(), mlxVp, mlfIndexRef(veh, ".P1_b"), F_1b, NULL));
    /*
     * 
     * t1=[F_1b;M_1bb];
     */
    mlfAssign(
      &t1, mlfVertcat(mlfHorzcat(F_1b, NULL), mlfHorzcat(M_1bb, NULL), NULL));
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
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(2.0))),
          mlfScalar(0.0),
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(2.0))),
          NULL),
        mlfHorzcat(
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(2.0))),
          mlfScalar(0.0),
          mlfUminus(mlfCos(mlfIndexRef(de, "(?)", mlfScalar(2.0)))),
          NULL),
        mlfDoubleMatrix(1, 3, __Array12_r, NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin2 middle point wrt c in b
     * v_2cb=vr(1:3)+vp(vr(4:6),veh.P2_b);
     */
    mlfAssign(
      &v_2cb,
      mlfPlus(
        mlfIndexRef(vr, "(?)", mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL)),
        mlfFeval(
          mclValueVarargout(),
          mlxVp,
          mlfIndexRef(
            vr, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL)),
          mlfIndexRef(veh, ".P2_b"),
          NULL)));
    /*
     * 
     * % relative velocity of fin2 middle point wrt c in 2
     * v_2c2=R_2b'*v_2cb;
     */
    mlfAssign(&v_2c2, mlfMtimes(mlfCtranspose(R_2b), v_2cb));
    /*
     * 
     * % attack angle
     * a2=atan2(v_2c2(3),v_2c2(1));
     */
    mlfAssign(
      &a2,
      mlfAtan2(
        mlfIndexRef(v_2c2, "(?)", mlfScalar(3.0)),
        mlfIndexRef(v_2c2, "(?)", mlfScalar(1.0))));
    /*
     * 
     * % fin2 wind frame rotation matrix
     * R_w2=[cos(a2) 0 -sin(a2); 0 1 0; sin(a2) 0 cos(a2)];
     */
    mlfAssign(
      &R_w2,
      mlfVertcat(
        mlfHorzcat(mlfCos(a2), mlfScalar(0.0), mlfUminus(mlfSin(a2)), NULL),
        mlfDoubleMatrix(1, 3, __Array13_r, NULL),
        mlfHorzcat(mlfSin(a2), mlfScalar(0.0), mlfCos(a2), NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a2);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, a2));
    /*
     * 
     * % damping forces on 2 wrt w 
     * F_2w=-0.5*veh.rho*veh.sw*(v_2c2(1)^2+v_2c2(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_2w,
      mlfMtimes(
        mlfMtimes(
          mlfFeval(
            mclValueVarargout(),
            mlxMtimes,
            mlfFeval(
              mclValueVarargout(),
              mlxMtimes,
              mlfScalar(-0.5),
              mlfIndexRef(veh, ".rho"),
              NULL),
            mlfIndexRef(veh, ".sw"),
            NULL),
          mlfPlus(
            mlfMpower(
              mlfIndexRef(v_2c2, "(?)", mlfScalar(1.0)), mlfScalar(2.0)),
            mlfMpower(
              mlfIndexRef(v_2c2, "(?)", mlfScalar(3.0)), mlfScalar(2.0)))),
        mlfVertcat(
          mlfHorzcat(cd, NULL),
          mlfDoubleMatrix(1, 1, __Array14_r, NULL),
          mlfHorzcat(cl, NULL),
          NULL)));
    /*
     * 
     * % damping forces on 2 wrt b 
     * F_2b=R_2b*R_w2*F_2w;
     */
    mlfAssign(&F_2b, mlfMtimes(mlfMtimes(R_2b, R_w2), F_2w));
    /*
     * 
     * % moments on 2 with pole in b wrt b
     * M_2bb=vp(veh.P2_b,F_2b);
     */
    mlfAssign(
      &M_2bb,
      mlfFeval(
        mclValueVarargout(), mlxVp, mlfIndexRef(veh, ".P2_b"), F_2b, NULL));
    /*
     * 
     * t2=[F_2b;M_2bb];
     */
    mlfAssign(
      &t2, mlfVertcat(mlfHorzcat(F_2b, NULL), mlfHorzcat(M_2bb, NULL), NULL));
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
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(3.0))),
          mlfScalar(0.0),
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(3.0))),
          NULL),
        mlfDoubleMatrix(1, 3, __Array15_r, NULL),
        mlfHorzcat(
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(3.0))),
          mlfScalar(0.0),
          mlfUminus(mlfCos(mlfIndexRef(de, "(?)", mlfScalar(3.0)))),
          NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin3 middle point wrt c in b
     * v_3cb=vr(1:3)+vp(vr(4:6),veh.P3_b);
     */
    mlfAssign(
      &v_3cb,
      mlfPlus(
        mlfIndexRef(vr, "(?)", mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL)),
        mlfFeval(
          mclValueVarargout(),
          mlxVp,
          mlfIndexRef(
            vr, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL)),
          mlfIndexRef(veh, ".P3_b"),
          NULL)));
    /*
     * 
     * % relative velocity of fin3 middle point wrt c in 3
     * v_3c3=R_3b'*v_3cb;
     */
    mlfAssign(&v_3c3, mlfMtimes(mlfCtranspose(R_3b), v_3cb));
    /*
     * 
     * % attack angle
     * a3=atan2(v_3c3(3),v_3c3(1));
     */
    mlfAssign(
      &a3,
      mlfAtan2(
        mlfIndexRef(v_3c3, "(?)", mlfScalar(3.0)),
        mlfIndexRef(v_3c3, "(?)", mlfScalar(1.0))));
    /*
     * 
     * % fin3 wind frame rotation matrix
     * R_w3=[cos(a3) 0 -sin(a3); 0 1 0; sin(a3) 0 cos(a3)];
     */
    mlfAssign(
      &R_w3,
      mlfVertcat(
        mlfHorzcat(mlfCos(a3), mlfScalar(0.0), mlfUminus(mlfSin(a3)), NULL),
        mlfDoubleMatrix(1, 3, __Array16_r, NULL),
        mlfHorzcat(mlfSin(a3), mlfScalar(0.0), mlfCos(a3), NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a3);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, a3));
    /*
     * 
     * % damping forces on 3 wrt w 
     * F_3w=-0.5*veh.rho*veh.sw*(v_3c3(1)^2+v_3c3(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_3w,
      mlfMtimes(
        mlfMtimes(
          mlfFeval(
            mclValueVarargout(),
            mlxMtimes,
            mlfFeval(
              mclValueVarargout(),
              mlxMtimes,
              mlfScalar(-0.5),
              mlfIndexRef(veh, ".rho"),
              NULL),
            mlfIndexRef(veh, ".sw"),
            NULL),
          mlfPlus(
            mlfMpower(
              mlfIndexRef(v_3c3, "(?)", mlfScalar(1.0)), mlfScalar(2.0)),
            mlfMpower(
              mlfIndexRef(v_3c3, "(?)", mlfScalar(3.0)), mlfScalar(2.0)))),
        mlfVertcat(
          mlfHorzcat(cd, NULL),
          mlfDoubleMatrix(1, 1, __Array17_r, NULL),
          mlfHorzcat(cl, NULL),
          NULL)));
    /*
     * 
     * % damping forces on 3 wrt b 
     * F_3b=R_3b*R_w3*F_3w;
     */
    mlfAssign(&F_3b, mlfMtimes(mlfMtimes(R_3b, R_w3), F_3w));
    /*
     * 
     * % moments on 3 with pole in b wrt b
     * M_3bb=vp(veh.P3_b,F_3b);
     */
    mlfAssign(
      &M_3bb,
      mlfFeval(
        mclValueVarargout(), mlxVp, mlfIndexRef(veh, ".P3_b"), F_3b, NULL));
    /*
     * 
     * t3=[F_3b;M_3bb];
     */
    mlfAssign(
      &t3, mlfVertcat(mlfHorzcat(F_3b, NULL), mlfHorzcat(M_3bb, NULL), NULL));
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
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(4.0))),
          mlfScalar(0.0),
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(4.0))),
          NULL),
        mlfHorzcat(
          mlfUminus(mlfSin(mlfIndexRef(de, "(?)", mlfScalar(4.0)))),
          mlfScalar(0.0),
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(4.0))),
          NULL),
        mlfDoubleMatrix(1, 3, __Array18_r, NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin4 middle point wrt c in b
     * v_4cb=vr(1:3)+vp(vr(4:6),veh.P4_b);
     */
    mlfAssign(
      &v_4cb,
      mlfPlus(
        mlfIndexRef(vr, "(?)", mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL)),
        mlfFeval(
          mclValueVarargout(),
          mlxVp,
          mlfIndexRef(
            vr, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL)),
          mlfIndexRef(veh, ".P4_b"),
          NULL)));
    /*
     * 
     * % relative velocity of fin4 middle point wrt c in 4
     * v_4c4=R_4b'*v_4cb;
     */
    mlfAssign(&v_4c4, mlfMtimes(mlfCtranspose(R_4b), v_4cb));
    /*
     * 
     * % attack angle
     * a4=atan2(v_4c4(3),v_4c4(1));
     */
    mlfAssign(
      &a4,
      mlfAtan2(
        mlfIndexRef(v_4c4, "(?)", mlfScalar(3.0)),
        mlfIndexRef(v_4c4, "(?)", mlfScalar(1.0))));
    /*
     * 
     * % fin4 wind frame rotation matrix
     * R_w4=[cos(a4) 0 -sin(a4); 0 1 0; sin(a4) 0 cos(a4)];
     */
    mlfAssign(
      &R_w4,
      mlfVertcat(
        mlfHorzcat(mlfCos(a4), mlfScalar(0.0), mlfUminus(mlfSin(a4)), NULL),
        mlfDoubleMatrix(1, 3, __Array19_r, NULL),
        mlfHorzcat(mlfSin(a4), mlfScalar(0.0), mlfCos(a4), NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a4);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, a4));
    /*
     * 
     * % damping forces on 4 wrt w 
     * F_4w=-0.5*veh.rho*veh.sw*(v_4c4(1)^2+v_4c4(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_4w,
      mlfMtimes(
        mlfMtimes(
          mlfFeval(
            mclValueVarargout(),
            mlxMtimes,
            mlfFeval(
              mclValueVarargout(),
              mlxMtimes,
              mlfScalar(-0.5),
              mlfIndexRef(veh, ".rho"),
              NULL),
            mlfIndexRef(veh, ".sw"),
            NULL),
          mlfPlus(
            mlfMpower(
              mlfIndexRef(v_4c4, "(?)", mlfScalar(1.0)), mlfScalar(2.0)),
            mlfMpower(
              mlfIndexRef(v_4c4, "(?)", mlfScalar(3.0)), mlfScalar(2.0)))),
        mlfVertcat(
          mlfHorzcat(cd, NULL),
          mlfDoubleMatrix(1, 1, __Array20_r, NULL),
          mlfHorzcat(cl, NULL),
          NULL)));
    /*
     * 
     * % damping forces on 4 wrt b 
     * F_4b=R_4b*R_w4*F_4w;
     */
    mlfAssign(&F_4b, mlfMtimes(mlfMtimes(R_4b, R_w4), F_4w));
    /*
     * 
     * % moments on 4 with pole in b wrt b
     * M_4bb=vp(veh.P4_b,F_4b);
     */
    mlfAssign(
      &M_4bb,
      mlfFeval(
        mclValueVarargout(), mlxVp, mlfIndexRef(veh, ".P4_b"), F_4b, NULL));
    /*
     * 
     * t4=[F_4b;M_4bb];
     */
    mlfAssign(
      &t4, mlfVertcat(mlfHorzcat(F_4b, NULL), mlfHorzcat(M_4bb, NULL), NULL));
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
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(5.0))),
          mlfScalar(0.0),
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(5.0))),
          NULL),
        mlfDoubleMatrix(1, 3, __Array21_r, NULL),
        mlfHorzcat(
          mlfUminus(mlfSin(mlfIndexRef(de, "(?)", mlfScalar(5.0)))),
          mlfScalar(0.0),
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(5.0))),
          NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin5 middle point wrt c in b
     * v_5cb=vr(1:3)+vp(vr(4:6),veh.P5_b);
     */
    mlfAssign(
      &v_5cb,
      mlfPlus(
        mlfIndexRef(vr, "(?)", mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL)),
        mlfFeval(
          mclValueVarargout(),
          mlxVp,
          mlfIndexRef(
            vr, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL)),
          mlfIndexRef(veh, ".P5_b"),
          NULL)));
    /*
     * 
     * % relative velocity of fin5 middle point wrt c in 5
     * v_5c5=R_5b'*v_5cb;
     */
    mlfAssign(&v_5c5, mlfMtimes(mlfCtranspose(R_5b), v_5cb));
    /*
     * 
     * % attack angle
     * a5=atan2(v_5c5(3),v_5c5(1));
     */
    mlfAssign(
      &a5,
      mlfAtan2(
        mlfIndexRef(v_5c5, "(?)", mlfScalar(3.0)),
        mlfIndexRef(v_5c5, "(?)", mlfScalar(1.0))));
    /*
     * 
     * % fin5 wind frame rotation matrix
     * R_w5=[cos(a5) 0 -sin(a5); 0 1 0; sin(a5) 0 cos(a5)];
     */
    mlfAssign(
      &R_w5,
      mlfVertcat(
        mlfHorzcat(mlfCos(a5), mlfScalar(0.0), mlfUminus(mlfSin(a5)), NULL),
        mlfDoubleMatrix(1, 3, __Array22_r, NULL),
        mlfHorzcat(mlfSin(a5), mlfScalar(0.0), mlfCos(a5), NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a5);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, a5));
    /*
     * 
     * % damping forces on 5 wrt w 
     * F_5w=-0.5*veh.rho*veh.st*(v_5c5(1)^2+v_5c5(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_5w,
      mlfMtimes(
        mlfMtimes(
          mlfFeval(
            mclValueVarargout(),
            mlxMtimes,
            mlfFeval(
              mclValueVarargout(),
              mlxMtimes,
              mlfScalar(-0.5),
              mlfIndexRef(veh, ".rho"),
              NULL),
            mlfIndexRef(veh, ".st"),
            NULL),
          mlfPlus(
            mlfMpower(
              mlfIndexRef(v_5c5, "(?)", mlfScalar(1.0)), mlfScalar(2.0)),
            mlfMpower(
              mlfIndexRef(v_5c5, "(?)", mlfScalar(3.0)), mlfScalar(2.0)))),
        mlfVertcat(
          mlfHorzcat(cd, NULL),
          mlfDoubleMatrix(1, 1, __Array23_r, NULL),
          mlfHorzcat(cl, NULL),
          NULL)));
    /*
     * 
     * % damping forces on 5 wrt b 
     * F_5b=R_5b*R_w5*F_5w;
     */
    mlfAssign(&F_5b, mlfMtimes(mlfMtimes(R_5b, R_w5), F_5w));
    /*
     * 
     * % moments on 5 with pole in b wrt b
     * M_5bb=vp(veh.P5_b,F_5b);
     */
    mlfAssign(
      &M_5bb,
      mlfFeval(
        mclValueVarargout(), mlxVp, mlfIndexRef(veh, ".P5_b"), F_5b, NULL));
    /*
     * 
     * t5=[F_5b;M_5bb];
     */
    mlfAssign(
      &t5, mlfVertcat(mlfHorzcat(F_5b, NULL), mlfHorzcat(M_5bb, NULL), NULL));
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
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(6.0))),
          mlfScalar(0.0),
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(6.0))),
          NULL),
        mlfHorzcat(
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(6.0))),
          mlfScalar(0.0),
          mlfUminus(mlfCos(mlfIndexRef(de, "(?)", mlfScalar(6.0)))),
          NULL),
        mlfDoubleMatrix(1, 3, __Array24_r, NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin6 middle point wrt c in b
     * v_6cb=vr(1:3)+vp(vr(4:6),veh.P6_b);
     */
    mlfAssign(
      &v_6cb,
      mlfPlus(
        mlfIndexRef(vr, "(?)", mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL)),
        mlfFeval(
          mclValueVarargout(),
          mlxVp,
          mlfIndexRef(
            vr, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL)),
          mlfIndexRef(veh, ".P6_b"),
          NULL)));
    /*
     * 
     * % relative velocity of fin6 middle point wrt c in 6
     * v_6c6=R_6b'*v_6cb;
     */
    mlfAssign(&v_6c6, mlfMtimes(mlfCtranspose(R_6b), v_6cb));
    /*
     * 
     * % attack angle
     * a6=atan2(v_6c6(3),v_6c6(1));
     */
    mlfAssign(
      &a6,
      mlfAtan2(
        mlfIndexRef(v_6c6, "(?)", mlfScalar(3.0)),
        mlfIndexRef(v_6c6, "(?)", mlfScalar(1.0))));
    /*
     * 
     * % fin6 wind frame rotation matrix
     * R_w6=[cos(a6) 0 -sin(a6); 0 1 0; sin(a6) 0 cos(a6)];
     */
    mlfAssign(
      &R_w6,
      mlfVertcat(
        mlfHorzcat(mlfCos(a6), mlfScalar(0.0), mlfUminus(mlfSin(a6)), NULL),
        mlfDoubleMatrix(1, 3, __Array25_r, NULL),
        mlfHorzcat(mlfSin(a6), mlfScalar(0.0), mlfCos(a6), NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a6);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, a6));
    /*
     * 
     * % damping forces on 6 wrt w 
     * F_6w=-0.5*veh.rho*veh.st*(v_6c6(1)^2+v_6c6(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_6w,
      mlfMtimes(
        mlfMtimes(
          mlfFeval(
            mclValueVarargout(),
            mlxMtimes,
            mlfFeval(
              mclValueVarargout(),
              mlxMtimes,
              mlfScalar(-0.5),
              mlfIndexRef(veh, ".rho"),
              NULL),
            mlfIndexRef(veh, ".st"),
            NULL),
          mlfPlus(
            mlfMpower(
              mlfIndexRef(v_6c6, "(?)", mlfScalar(1.0)), mlfScalar(2.0)),
            mlfMpower(
              mlfIndexRef(v_6c6, "(?)", mlfScalar(3.0)), mlfScalar(2.0)))),
        mlfVertcat(
          mlfHorzcat(cd, NULL),
          mlfDoubleMatrix(1, 1, __Array26_r, NULL),
          mlfHorzcat(cl, NULL),
          NULL)));
    /*
     * 
     * % damping forces on 6 wrt b 
     * F_6b=R_6b*R_w6*F_6w;
     */
    mlfAssign(&F_6b, mlfMtimes(mlfMtimes(R_6b, R_w6), F_6w));
    /*
     * 
     * % moments on 6 with pole in b wrt b
     * M_6bb=vp(veh.P6_b,F_6b);
     */
    mlfAssign(
      &M_6bb,
      mlfFeval(
        mclValueVarargout(), mlxVp, mlfIndexRef(veh, ".P6_b"), F_6b, NULL));
    /*
     * 
     * t6=[F_6b;M_6bb];
     */
    mlfAssign(
      &t6, mlfVertcat(mlfHorzcat(F_6b, NULL), mlfHorzcat(M_6bb, NULL), NULL));
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
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(7.0))),
          mlfScalar(0.0),
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(7.0))),
          NULL),
        mlfDoubleMatrix(1, 3, __Array27_r, NULL),
        mlfHorzcat(
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(7.0))),
          mlfScalar(0.0),
          mlfUminus(mlfCos(mlfIndexRef(de, "(?)", mlfScalar(7.0)))),
          NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin7 middle point wrt c in b
     * v_7cb=vr(1:3)+vp(vr(4:6),veh.P7_b);
     */
    mlfAssign(
      &v_7cb,
      mlfPlus(
        mlfIndexRef(vr, "(?)", mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL)),
        mlfFeval(
          mclValueVarargout(),
          mlxVp,
          mlfIndexRef(
            vr, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL)),
          mlfIndexRef(veh, ".P7_b"),
          NULL)));
    /*
     * 
     * % relative velocity of fin7 middle point wrt c in 7
     * v_7c7=R_7b'*v_7cb;
     */
    mlfAssign(&v_7c7, mlfMtimes(mlfCtranspose(R_7b), v_7cb));
    /*
     * 
     * % attack angle
     * a7=atan2(v_7c7(3),v_7c7(1));
     */
    mlfAssign(
      &a7,
      mlfAtan2(
        mlfIndexRef(v_7c7, "(?)", mlfScalar(3.0)),
        mlfIndexRef(v_7c7, "(?)", mlfScalar(1.0))));
    /*
     * 
     * % fin7 wind frame rotation matrix
     * R_w7=[cos(a7) 0 -sin(a7); 0 1 0; sin(a7) 0 cos(a7)];
     */
    mlfAssign(
      &R_w7,
      mlfVertcat(
        mlfHorzcat(mlfCos(a7), mlfScalar(0.0), mlfUminus(mlfSin(a7)), NULL),
        mlfDoubleMatrix(1, 3, __Array28_r, NULL),
        mlfHorzcat(mlfSin(a7), mlfScalar(0.0), mlfCos(a7), NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a7);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, a7));
    /*
     * 
     * % damping forces on 7 wrt w 
     * F_7w=-0.5*veh.rho*veh.st*(v_7c7(1)^2+v_7c7(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_7w,
      mlfMtimes(
        mlfMtimes(
          mlfFeval(
            mclValueVarargout(),
            mlxMtimes,
            mlfFeval(
              mclValueVarargout(),
              mlxMtimes,
              mlfScalar(-0.5),
              mlfIndexRef(veh, ".rho"),
              NULL),
            mlfIndexRef(veh, ".st"),
            NULL),
          mlfPlus(
            mlfMpower(
              mlfIndexRef(v_7c7, "(?)", mlfScalar(1.0)), mlfScalar(2.0)),
            mlfMpower(
              mlfIndexRef(v_7c7, "(?)", mlfScalar(3.0)), mlfScalar(2.0)))),
        mlfVertcat(
          mlfHorzcat(cd, NULL),
          mlfDoubleMatrix(1, 1, __Array29_r, NULL),
          mlfHorzcat(cl, NULL),
          NULL)));
    /*
     * 
     * % damping forces on 7 wrt b 
     * F_7b=R_7b*R_w7*F_7w;
     */
    mlfAssign(&F_7b, mlfMtimes(mlfMtimes(R_7b, R_w7), F_7w));
    /*
     * 
     * % moments on 7 with pole in b wrt b
     * M_7bb=vp(veh.P7_b,F_7b);
     */
    mlfAssign(
      &M_7bb,
      mlfFeval(
        mclValueVarargout(), mlxVp, mlfIndexRef(veh, ".P7_b"), F_7b, NULL));
    /*
     * 
     * t7=[F_7b;M_7bb];
     */
    mlfAssign(
      &t7, mlfVertcat(mlfHorzcat(F_7b, NULL), mlfHorzcat(M_7bb, NULL), NULL));
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
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(8.0))),
          mlfScalar(0.0),
          mlfSin(mlfIndexRef(de, "(?)", mlfScalar(8.0))),
          NULL),
        mlfHorzcat(
          mlfUminus(mlfSin(mlfIndexRef(de, "(?)", mlfScalar(8.0)))),
          mlfScalar(0.0),
          mlfCos(mlfIndexRef(de, "(?)", mlfScalar(8.0))),
          NULL),
        mlfDoubleMatrix(1, 3, __Array30_r, NULL),
        NULL));
    /*
     * 
     * % relative velocity of fin8 middle point wrt c in b
     * v_8cb=vr(1:3)+vp(vr(4:6),veh.P8_b);
     */
    mlfAssign(
      &v_8cb,
      mlfPlus(
        mlfIndexRef(vr, "(?)", mlfColon(mlfScalar(1.0), mlfScalar(3.0), NULL)),
        mlfFeval(
          mclValueVarargout(),
          mlxVp,
          mlfIndexRef(
            vr, "(?)", mlfColon(mlfScalar(4.0), mlfScalar(6.0), NULL)),
          mlfIndexRef(veh, ".P8_b"),
          NULL)));
    /*
     * 
     * % relative velocity of fin8 middle point wrt c in 8
     * v_8c8=R_8b'*v_8cb;
     */
    mlfAssign(&v_8c8, mlfMtimes(mlfCtranspose(R_8b), v_8cb));
    /*
     * 
     * % attack angle
     * a8=atan2(v_8c8(3),v_8c8(1));
     */
    mlfAssign(
      &a8,
      mlfAtan2(
        mlfIndexRef(v_8c8, "(?)", mlfScalar(3.0)),
        mlfIndexRef(v_8c8, "(?)", mlfScalar(1.0))));
    /*
     * 
     * % fin8 wind frame rotation matrix
     * R_w8=[cos(a8) 0 -sin(a8); 0 1 0; sin(a8) 0 cos(a8)];
     */
    mlfAssign(
      &R_w8,
      mlfVertcat(
        mlfHorzcat(mlfCos(a8), mlfScalar(0.0), mlfUminus(mlfSin(a8)), NULL),
        mlfDoubleMatrix(1, 3, __Array31_r, NULL),
        mlfHorzcat(mlfSin(a8), mlfScalar(0.0), mlfCos(a8), NULL),
        NULL));
    /*
     * 
     * % cl and cd computation
     * [cl,cd]=a2clcd(a8);
     */
    mlfAssign(&cl, mlfA2clcd(&cd, a8));
    /*
     * 
     * % damping forces on 8 wrt w 
     * F_8w=-0.5*veh.rho*veh.st*(v_8c8(1)^2+v_8c8(3)^2)*[cd; 0; cl];
     */
    mlfAssign(
      &F_8w,
      mlfMtimes(
        mlfMtimes(
          mlfFeval(
            mclValueVarargout(),
            mlxMtimes,
            mlfFeval(
              mclValueVarargout(),
              mlxMtimes,
              mlfScalar(-0.5),
              mlfIndexRef(veh, ".rho"),
              NULL),
            mlfIndexRef(veh, ".st"),
            NULL),
          mlfPlus(
            mlfMpower(
              mlfIndexRef(v_8c8, "(?)", mlfScalar(1.0)), mlfScalar(2.0)),
            mlfMpower(
              mlfIndexRef(v_8c8, "(?)", mlfScalar(3.0)), mlfScalar(2.0)))),
        mlfVertcat(
          mlfHorzcat(cd, NULL),
          mlfDoubleMatrix(1, 1, __Array32_r, NULL),
          mlfHorzcat(cl, NULL),
          NULL)));
    /*
     * 
     * % damping forces on 8 wrt b 
     * F_8b=R_8b*R_w8*F_8w;
     */
    mlfAssign(&F_8b, mlfMtimes(mlfMtimes(R_8b, R_w8), F_8w));
    /*
     * 
     * % moments on 8 with pole in b wrt b
     * M_8bb=vp(veh.P8_b,F_8b);
     */
    mlfAssign(
      &M_8bb,
      mlfFeval(
        mclValueVarargout(), mlxVp, mlfIndexRef(veh, ".P8_b"), F_8b, NULL));
    /*
     * 
     * t8=[F_8b;M_8bb];
     */
    mlfAssign(
      &t8, mlfVertcat(mlfHorzcat(F_8b, NULL), mlfHorzcat(M_8bb, NULL), NULL));
    /*
     * 
     * % --------------------------------------------------------------
     * % resulting hydrodynamic force and moment with pole in b wrt b
     * 
     * td=tf+t1+t2+t3+t4+t5+t6+t7+t8; 
     */
    mlfAssign(
      &td,
      mlfPlus(
        mlfPlus(
          mlfPlus(
            mlfPlus(mlfPlus(mlfPlus(mlfPlus(mlfPlus(tf, t1), t2), t3), t4), t5),
            t6),
          t7),
        t8));
    mclValidateOutputs("tau_damp", 1, nargout_, &td);
    mxDestroyArray(F_1b);
    mxDestroyArray(F_1w);
    mxDestroyArray(F_2b);
    mxDestroyArray(F_2w);
    mxDestroyArray(F_3b);
    mxDestroyArray(F_3w);
    mxDestroyArray(F_4b);
    mxDestroyArray(F_4w);
    mxDestroyArray(F_5b);
    mxDestroyArray(F_5w);
    mxDestroyArray(F_6b);
    mxDestroyArray(F_6w);
    mxDestroyArray(F_7b);
    mxDestroyArray(F_7w);
    mxDestroyArray(F_8b);
    mxDestroyArray(F_8w);
    mxDestroyArray(F_Bb);
    mxDestroyArray(F_Bw);
    mxDestroyArray(M_1bb);
    mxDestroyArray(M_2bb);
    mxDestroyArray(M_3bb);
    mxDestroyArray(M_4bb);
    mxDestroyArray(M_5bb);
    mxDestroyArray(M_6bb);
    mxDestroyArray(M_7bb);
    mxDestroyArray(M_8bb);
    mxDestroyArray(M_Bbb);
    mxDestroyArray(Pf_b);
    mxDestroyArray(R_1b);
    mxDestroyArray(R_2b);
    mxDestroyArray(R_3b);
    mxDestroyArray(R_4b);
    mxDestroyArray(R_5b);
    mxDestroyArray(R_6b);
    mxDestroyArray(R_7b);
    mxDestroyArray(R_8b);
    mxDestroyArray(R_fb);
    mxDestroyArray(R_w1);
    mxDestroyArray(R_w2);
    mxDestroyArray(R_w3);
    mxDestroyArray(R_w4);
    mxDestroyArray(R_w5);
    mxDestroyArray(R_w6);
    mxDestroyArray(R_w7);
    mxDestroyArray(R_w8);
    mxDestroyArray(R_wf);
    mxDestroyArray(a1);
    mxDestroyArray(a2);
    mxDestroyArray(a3);
    mxDestroyArray(a4);
    mxDestroyArray(a5);
    mxDestroyArray(a6);
    mxDestroyArray(a7);
    mxDestroyArray(a8);
    mxDestroyArray(af);
    mxDestroyArray(cd);
    mxDestroyArray(cl);
    mxDestroyArray(i_fb);
    mxDestroyArray(k_fb);
    mxDestroyArray(sf);
    mxDestroyArray(t1);
    mxDestroyArray(t2);
    mxDestroyArray(t3);
    mxDestroyArray(t4);
    mxDestroyArray(t5);
    mxDestroyArray(t6);
    mxDestroyArray(t7);
    mxDestroyArray(t8);
    mxDestroyArray(tf);
    mxDestroyArray(v_1c1);
    mxDestroyArray(v_1cb);
    mxDestroyArray(v_2c2);
    mxDestroyArray(v_2cb);
    mxDestroyArray(v_3c3);
    mxDestroyArray(v_3cb);
    mxDestroyArray(v_4c4);
    mxDestroyArray(v_4cb);
    mxDestroyArray(v_5c5);
    mxDestroyArray(v_5cb);
    mxDestroyArray(v_6c6);
    mxDestroyArray(v_6cb);
    mxDestroyArray(v_7c7);
    mxDestroyArray(v_7cb);
    mxDestroyArray(v_8c8);
    mxDestroyArray(v_8cb);
    mxDestroyArray(v_Bcb);
    mxDestroyArray(v_Bcf);
    mxDestroyArray(xcp);
    return td;
}

/*
 * The function "mlfTau_damp" contains the normal interface for the "tau_damp"
 * M-function from file "E:\RTW\NEW3\TAU_DAMP.M" (lines 1-307). This function
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
 * M-function from file "E:\RTW\NEW3\TAU_DAMP.M" (lines 1-307). The feval
 * function calls the implementation version of tau_damp through this function.
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxTau_damp(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[3];
    mxArray * mplhs[1];
    int i;
    if (nlhs > 1) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: tau_damp Line: 1 Column:"
            " 0 The function \"tau_damp\" was called with m"
            "ore than the declared number of outputs (1)"));
    }
    if (nrhs > 3) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: tau_damp Line: 1 Column:"
            " 0 The function \"tau_damp\" was called with m"
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
    mplhs[0] = Mtau_damp(nlhs, mprhs[0], mprhs[1], mprhs[2]);
    mlfRestorePreviousContext(0, 3, mprhs[0], mprhs[1], mprhs[2]);
    plhs[0] = mplhs[0];
}
