/*
 * MATLAB Compiler: 2.1
 * Date: Wed Oct 17 16:15:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-h" "-x" "-W" "mex"
 * "-L" "C" "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "vxdot.m" 
 */
#include "vp.h"
#include "libmatlbm.h"

static mxChar _array1_[122] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'v', 'p', ' ', 'L', 'i',
                                'n', 'e', ':', ' ', '1', ' ', 'C', 'o', 'l',
                                'u', 'm', 'n', ':', ' ', '1', ' ', 'T', 'h',
                                'e', ' ', 'f', 'u', 'n', 'c', 't', 'i', 'o',
                                'n', ' ', '"', 'v', 'p', '"', ' ', 'w', 'a',
                                's', ' ', 'c', 'a', 'l', 'l', 'e', 'd', ' ',
                                'w', 'i', 't', 'h', ' ', 'm', 'o', 'r', 'e',
                                ' ', 't', 'h', 'a', 'n', ' ', 't', 'h', 'e',
                                ' ', 'd', 'e', 'c', 'l', 'a', 'r', 'e', 'd',
                                ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o',
                                'f', ' ', 'o', 'u', 't', 'p', 'u', 't', 's',
                                ' ', '(', '1', ')', '.' };
static mxArray * _mxarray0_;

static mxChar _array3_[121] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'v', 'p', ' ', 'L', 'i',
                                'n', 'e', ':', ' ', '1', ' ', 'C', 'o', 'l',
                                'u', 'm', 'n', ':', ' ', '1', ' ', 'T', 'h',
                                'e', ' ', 'f', 'u', 'n', 'c', 't', 'i', 'o',
                                'n', ' ', '"', 'v', 'p', '"', ' ', 'w', 'a',
                                's', ' ', 'c', 'a', 'l', 'l', 'e', 'd', ' ',
                                'w', 'i', 't', 'h', ' ', 'm', 'o', 'r', 'e',
                                ' ', 't', 'h', 'a', 'n', ' ', 't', 'h', 'e',
                                ' ', 'd', 'e', 'c', 'l', 'a', 'r', 'e', 'd',
                                ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o',
                                'f', ' ', 'i', 'n', 'p', 'u', 't', 's', ' ',
                                '(', '2', ')', '.' };
static mxArray * _mxarray2_;
static mxArray * _mxarray4_;

void InitializeModule_vp(void) {
    _mxarray0_ = mclInitializeString(122, _array1_);
    _mxarray2_ = mclInitializeString(121, _array3_);
    _mxarray4_ = mclInitializeDouble(0.0);
}

void TerminateModule_vp(void) {
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mvp(int nargout_, mxArray * x, mxArray * y);

_mexLocalFunctionTable _local_function_table_vp
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfVp" contains the normal interface for the "vp" M-function
 * from file "D:\RTW\new3\vp.m" (lines 1-13). This function processes any input
 * arguments and passes them to the implementation version of the function,
 * appearing above.
 */
mxArray * mlfVp(mxArray * x, mxArray * y) {
    int nargout = 1;
    mxArray * z = mclGetUninitializedArray();
    mlfEnterNewContext(0, 2, x, y);
    z = Mvp(nargout, x, y);
    mlfRestorePreviousContext(0, 2, x, y);
    return mlfReturnValue(z);
}

/*
 * The function "mlxVp" contains the feval interface for the "vp" M-function
 * from file "D:\RTW\new3\vp.m" (lines 1-13). The feval function calls the
 * implementation version of vp through this function. This function processes
 * any input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
void mlxVp(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
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
    mplhs[0] = Mvp(nlhs, mprhs[0], mprhs[1]);
    mlfRestorePreviousContext(0, 2, mprhs[0], mprhs[1]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mvp" is the implementation version of the "vp" M-function from
 * file "D:\RTW\new3\vp.m" (lines 1-13). It contains the actual compiled code
 * for that M-function. It is a static function and must only be called from
 * one of the interface functions, appearing below.
 */
/*
 * function z=vp(x,y)
 */
static mxArray * Mvp(int nargout_, mxArray * x, mxArray * y) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_vp);
    int nargin_ = mclNargin(2, x, y, NULL);
    mxArray * z = mclGetUninitializedArray();
    mclCopyArray(&x);
    mclCopyArray(&y);
    /*
     * 
     * % z=vp(x,y); z = 3d cross product of x and y
     * % vp(x) is the 3d cross product matrix : vp(x)*y=vp(x,y).
     * %
     * % by Giampiero Campa.  
     * 
     * z=[  0    -x(3)   x(2);
     */
    mlfAssign(
      &z,
      mlfVertcat(
        mlfHorzcat(
          _mxarray4_,
          mclUminus(mclVe(mclIntArrayRef1(mclVsa(x, "x"), 3))),
          mclVe(mclIntArrayRef1(mclVsa(x, "x"), 2)),
          NULL),
        mlfHorzcat(
          mclVe(mclIntArrayRef1(mclVsa(x, "x"), 3)),
          _mxarray4_,
          mclUminus(mclVe(mclIntArrayRef1(mclVsa(x, "x"), 1))),
          NULL),
        mlfHorzcat(
          mclUminus(mclVe(mclIntArrayRef1(mclVsa(x, "x"), 2))),
          mclVe(mclIntArrayRef1(mclVsa(x, "x"), 1)),
          _mxarray4_,
          NULL),
        NULL));
    /*
     * x(3)    0    -x(1);
     * -x(2)   x(1)    0   ];
     * 
     * if nargin>1, z=z*y; end
     */
    if (nargin_ > 1) {
        mlfAssign(&z, mclMtimes(mclVv(z, "z"), mclVa(y, "y")));
    }
    mclValidateOutput(z, 1, nargout_, "z", "vp");
    mxDestroyArray(y);
    mxDestroyArray(x);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return z;
}
