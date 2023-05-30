/*
 * MATLAB Compiler: 2.1
 * Date: Wed Oct 17 16:15:00 2001
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-h" "-x" "-W" "mex"
 * "-L" "C" "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "vxdot.m" 
 */
#include "rpy2R_eb.h"
#include "libmatlbm.h"

static mxChar _array1_[134] = { 'R', 'u', 'n', '-', 't', 'i', 'm', 'e', ' ',
                                'E', 'r', 'r', 'o', 'r', ':', ' ', 'F', 'i',
                                'l', 'e', ':', ' ', 'r', 'p', 'y', '2', 'R',
                                '_', 'e', 'b', ' ', 'L', 'i', 'n', 'e', ':',
                                ' ', '1', ' ', 'C', 'o', 'l', 'u', 'm', 'n',
                                ':', ' ', '1', ' ', 'T', 'h', 'e', ' ', 'f',
                                'u', 'n', 'c', 't', 'i', 'o', 'n', ' ', '"',
                                'r', 'p', 'y', '2', 'R', '_', 'e', 'b', '"',
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
                                'l', 'e', ':', ' ', 'r', 'p', 'y', '2', 'R',
                                '_', 'e', 'b', ' ', 'L', 'i', 'n', 'e', ':',
                                ' ', '1', ' ', 'C', 'o', 'l', 'u', 'm', 'n',
                                ':', ' ', '1', ' ', 'T', 'h', 'e', ' ', 'f',
                                'u', 'n', 'c', 't', 'i', 'o', 'n', ' ', '"',
                                'r', 'p', 'y', '2', 'R', '_', 'e', 'b', '"',
                                ' ', 'w', 'a', 's', ' ', 'c', 'a', 'l', 'l',
                                'e', 'd', ' ', 'w', 'i', 't', 'h', ' ', 'm',
                                'o', 'r', 'e', ' ', 't', 'h', 'a', 'n', ' ',
                                't', 'h', 'e', ' ', 'd', 'e', 'c', 'l', 'a',
                                'r', 'e', 'd', ' ', 'n', 'u', 'm', 'b', 'e',
                                'r', ' ', 'o', 'f', ' ', 'i', 'n', 'p', 'u',
                                't', 's', ' ', '(', '1', ')', '.' };
static mxArray * _mxarray2_;

void InitializeModule_rpy2R_eb(void) {
    _mxarray0_ = mclInitializeString(134, _array1_);
    _mxarray2_ = mclInitializeString(133, _array3_);
}

void TerminateModule_rpy2R_eb(void) {
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static mxArray * Mrpy2R_eb(int nargout_, mxArray * rpy);

_mexLocalFunctionTable _local_function_table_rpy2R_eb
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfRpy2R_eb" contains the normal interface for the "rpy2R_eb"
 * M-function from file "D:\RTW\new3\rpy2R_eb.m" (lines 1-16). This function
 * processes any input arguments and passes them to the implementation version
 * of the function, appearing above.
 */
mxArray * mlfRpy2R_eb(mxArray * rpy) {
    int nargout = 1;
    mxArray * R_eb = mclGetUninitializedArray();
    mlfEnterNewContext(0, 1, rpy);
    R_eb = Mrpy2R_eb(nargout, rpy);
    mlfRestorePreviousContext(0, 1, rpy);
    return mlfReturnValue(R_eb);
}

/*
 * The function "mlxRpy2R_eb" contains the feval interface for the "rpy2R_eb"
 * M-function from file "D:\RTW\new3\rpy2R_eb.m" (lines 1-16). The feval
 * function calls the implementation version of rpy2R_eb through this function.
 * This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxRpy2R_eb(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
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
    mplhs[0] = Mrpy2R_eb(nlhs, mprhs[0]);
    mlfRestorePreviousContext(0, 1, mprhs[0]);
    plhs[0] = mplhs[0];
}

/*
 * The function "Mrpy2R_eb" is the implementation version of the "rpy2R_eb"
 * M-function from file "D:\RTW\new3\rpy2R_eb.m" (lines 1-16). It contains the
 * actual compiled code for that M-function. It is a static function and must
 * only be called from one of the interface functions, appearing below.
 */
/*
 * function R_eb=rpy2R_eb(rpy)
 */
static mxArray * Mrpy2R_eb(int nargout_, mxArray * rpy) {
    mexLocalFunctionTable save_local_function_table_ = mclSetCurrentLocalFunctionTable(
                                                         &_local_function_table_rpy2R_eb);
    mxArray * R_eb = mclGetUninitializedArray();
    mxArray * cp = mclGetUninitializedArray();
    mxArray * sp = mclGetUninitializedArray();
    mxArray * ct = mclGetUninitializedArray();
    mxArray * st = mclGetUninitializedArray();
    mxArray * cf = mclGetUninitializedArray();
    mxArray * sf = mclGetUninitializedArray();
    mclCopyArray(&rpy);
    /*
     * 
     * % R_eb=rpy2R_eb(rpy), computes the rotation matrix of e-frame wrt b-frame,
     * % given in input a vector containing roll, pitch and yaw angles.
     * 
     * sf = sin(rpy(1));
     */
    mlfAssign(&sf, mlfSin(mclVe(mclIntArrayRef1(mclVsa(rpy, "rpy"), 1))));
    /*
     * cf = cos(rpy(1));
     */
    mlfAssign(&cf, mlfCos(mclVe(mclIntArrayRef1(mclVsa(rpy, "rpy"), 1))));
    /*
     * st = sin(rpy(2));
     */
    mlfAssign(&st, mlfSin(mclVe(mclIntArrayRef1(mclVsa(rpy, "rpy"), 2))));
    /*
     * ct = cos(rpy(2));
     */
    mlfAssign(&ct, mlfCos(mclVe(mclIntArrayRef1(mclVsa(rpy, "rpy"), 2))));
    /*
     * sp = sin(rpy(3));
     */
    mlfAssign(&sp, mlfSin(mclVe(mclIntArrayRef1(mclVsa(rpy, "rpy"), 3))));
    /*
     * cp = cos(rpy(3));
     */
    mlfAssign(&cp, mlfCos(mclVe(mclIntArrayRef1(mclVsa(rpy, "rpy"), 3))));
    /*
     * 
     * R_eb = [         +ct*cp               +ct*sp         -st
     */
    mlfAssign(
      &R_eb,
      mlfVertcat(
        mlfHorzcat(
          mclMtimes(mclUplus(mclVv(ct, "ct")), mclVv(cp, "cp")),
          mclMtimes(mclUplus(mclVv(ct, "ct")), mclVv(sp, "sp")),
          mclUminus(mclVv(st, "st")),
          NULL),
        mlfHorzcat(
          mclMinus(
            mclMtimes(
              mclMtimes(mclUplus(mclVv(sf, "sf")), mclVv(st, "st")),
              mclVv(cp, "cp")),
            mclMtimes(mclVv(cf, "cf"), mclVv(sp, "sp"))),
          mclPlus(
            mclMtimes(
              mclMtimes(mclUplus(mclVv(sf, "sf")), mclVv(st, "st")),
              mclVv(sp, "sp")),
            mclMtimes(mclVv(cf, "cf"), mclVv(cp, "cp"))),
          mclMtimes(mclUplus(mclVv(sf, "sf")), mclVv(ct, "ct")),
          NULL),
        mlfHorzcat(
          mclPlus(
            mclMtimes(
              mclMtimes(mclUplus(mclVv(cf, "cf")), mclVv(st, "st")),
              mclVv(cp, "cp")),
            mclMtimes(mclVv(sf, "sf"), mclVv(sp, "sp"))),
          mclMinus(
            mclMtimes(
              mclMtimes(mclUplus(mclVv(cf, "cf")), mclVv(st, "st")),
              mclVv(sp, "sp")),
            mclMtimes(mclVv(sf, "sf"), mclVv(cp, "cp"))),
          mclMtimes(mclUplus(mclVv(cf, "cf")), mclVv(ct, "ct")),
          NULL),
        NULL));
    mclValidateOutput(R_eb, 1, nargout_, "R_eb", "rpy2R_eb");
    mxDestroyArray(sf);
    mxDestroyArray(cf);
    mxDestroyArray(st);
    mxDestroyArray(ct);
    mxDestroyArray(sp);
    mxDestroyArray(cp);
    mxDestroyArray(rpy);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
    return R_eb;
    /*
     * +sf*st*cp-cf*sp      +sf*st*sp+cf*cp      +sf*ct
     * +cf*st*cp+sf*sp      +cf*st*sp-sf*cp      +cf*ct ];
     */
}
