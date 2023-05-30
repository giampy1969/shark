/*
 * MATLAB Compiler: 2.1
 * Date: Wed Oct 17 16:15:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-h" "-x" "-W" "mex"
 * "-L" "C" "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "vxdot.m" 
 */
#include "rpy2J.h"
#include "libmatlbm.h"
#include "rpy2R_eb.h"

static mxChar _array1_[128] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'r', 'p', 'y', '2', 'J',
                                ' ', 'L', 'i', 'n', 'e', ':', ' ', '1', ' ',
                                'C', 'o', 'l', 'u', 'm', 'n', ':', ' ', '1',
                                ' ', 'T', 'h', 'e', ' ', 'f', 'u', 'n', 'c',
                                't', 'i', 'o', 'n', ' ', '"', 'r', 'p', 'y',
                                '2', 'J', '"', ' ', 'w', 'a', 's', ' ', 'c',
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
                                'l', 'e', ':', ' ', 'r', 'p', 'y', '2', 'J',
                                ' ', 'L', 'i', 'n', 'e', ':', ' ', '1', ' ',
                                'C', 'o', 'l', 'u', 'm', 'n', ':', ' ', '1',
                                ' ', 'T', 'h', 'e', ' ', 'f', 'u', 'n', 'c',
                                't', 'i', 'o', 'n', ' ', '"', 'r', 'p', 'y',
                                '2', 'J', '"', ' ', 'w', 'a', 's', ' ', 'c',
                                'a', 'l', 'l', 'e', 'd', ' ', 'w', 'i', 't',
                                'h', ' ', 'm', 'o', 'r', 'e', ' ', 't', 'h',
                                'a', 'n', ' ', 't', 'h', 'e', ' ', 'd', 'e',
                                'c', 'l', 'a', 'r', 'e', 'd', ' ', 'n', 'u',
                                'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ', 'i',
                                'n', 'p', 'u', 't', 's', ' ', '(', '1', ')',
                                '.' };
static mxArray * _mxarray2_;
static mxArray * _mxarray4_;
static mxArray * _mxarray5_;
static mxArray * _mxarray6_;

void InitializeModule_rpy2J(void) {
    _mxarray0_ = mclInitializeString(128, _array1_);
    _mxarray2_ = mclInitializeString(127, _array3_);
    _mxarray4_ = mclInitializeDouble(3.0);
    _mxarray5_ = mclInitializeDouble(1.0);
    _mxarray6_ = mclInitializeDouble(0.0);
}

void TerminateModule_rpy2J(void) {
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray5_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mrpy2J(int nargout_, mxArray * rpy);

_mexLocalFunctionTable _local_function_table_rpy2J
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfRpy2J" contains the normal interface for the "rpy2J"
 * M-function from file "D:\RTW\new3\rpy2J.m" (lines 1-14). This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
mxArray * mlfRpy2J(mxArray * rpy) {
    int nargout = 1;
    mxArray * J = mclGetUninitializedArray();
    mlfEnterNewContext(0, 1, rpy);
    J = Mrpy2J(nargout, rpy);
    mlfRestorePreviousContext(0, 1, rpy);
    return mlfReturnValue(J);
}

/*
 * The function "mlxRpy2J" contains the feval interface for the "rpy2J"
 * M-function from file "D:\RTW\new3\rpy2J.m" (lines 1-14). The feval function
 * calls the implementation version of rpy2J through this function. This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
void mlxRpy2J(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
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
    mplhs[0] = Mrpy2J(nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mrpy2J" is the implementation version of the "rpy2J"
 * M-function from file "D:\RTW\new3\rpy2J.m" (lines 1-14). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function J=rpy2J(rpy)
 */
static mxArray * Mrpy2J(int nargout_, mxArray * rpy) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_rpy2J);
    mxArray * J = mclGetUninitializedArray();
    mxArray * ct = mclGetUninitializedArray();
    mxArray * tt = mclGetUninitializedArray();
    mxArray * cf = mclGetUninitializedArray();
    mxArray * sf = mclGetUninitializedArray();
    mclCopyArray(&rpy);
    /*
     * 
     * % J=rpy2J(rpy); computes generalised Jacobian matrix which 
     * % transforms ni into eta derivatives, given in input roll pitch
     * % and yaw angles.
     * 
     * sf = sin(rpy(1));
     */
    mlfAssign(&sf, mlfSin(mclVe(mclIntArrayRef1(mclVsa(rpy, "rpy"), 1))));
    /*
     * cf = cos(rpy(1));
     */
    mlfAssign(&cf, mlfCos(mclVe(mclIntArrayRef1(mclVsa(rpy, "rpy"), 1))));
    /*
     * tt = tan(rpy(2));
     */
    mlfAssign(&tt, mlfTan(mclVe(mclIntArrayRef1(mclVsa(rpy, "rpy"), 2))));
    /*
     * ct = cos(rpy(2));
     */
    mlfAssign(&ct, mlfCos(mclVe(mclIntArrayRef1(mclVsa(rpy, "rpy"), 2))));
    /*
     * 
     * J = [ rpy2R_eb(rpy)' zeros(3,3)
     */
    mlfAssign(
      &J,
      mlfVertcat(
        mlfHorzcat(
          mlfCtranspose(mclVe(mlfRpy2R_eb(mclVa(rpy, "rpy")))),
          mclVe(mlfZeros(_mxarray4_, _mxarray4_, NULL)),
          NULL),
        mlfHorzcat(
          mclVe(mlfZeros(_mxarray4_, _mxarray4_, NULL)),
          mlfVertcat(
            mlfHorzcat(
              _mxarray5_,
              mclMtimes(mclVv(sf, "sf"), mclVv(tt, "tt")),
              mclMtimes(mclVv(cf, "cf"), mclVv(tt, "tt")),
              NULL),
            mlfHorzcat(
              _mxarray6_, mclVv(cf, "cf"), mclUminus(mclVv(sf, "sf")), NULL),
            mlfHorzcat(
              _mxarray6_,
              mclMrdivide(mclVv(sf, "sf"), mclVv(ct, "ct")),
              mclMrdivide(mclVv(cf, "cf"), mclVv(ct, "ct")),
              NULL),
            NULL),
          NULL),
        NULL));
    mclValidateOutput(J, 1, nargout_, "J", "rpy2J");
    mxDestroyArray(sf);
    mxDestroyArray(cf);
    mxDestroyArray(tt);
    mxDestroyArray(ct);
    mxDestroyArray(rpy);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return J;
    /*
     * zeros(3,3)     [1 sf*tt cf*tt; 0 cf -sf; 0 sf/ct cf/ct]  ];
     */
}
