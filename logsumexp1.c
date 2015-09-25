#include "mex.h"
#include "matrix.h"
#include <math.h>

/*   Copyright (c) 2015, Lloyd T. Elliott.   */

int validate(int nlhs,
        mxArray *plhs[],
        int nrhs,
        const mxArray *prhs[]) {
    
    if (nrhs != 1) {
        mexErrMsgTxt("Error: logsumexp1 requires 1 argument");
        return 1;
    } else if (nlhs == 0) {
        return 1;
    } else if (nlhs > 1) {
        mexErrMsgTxt("Error: logsumexp1 requires at most 1 output argument");
        return 1;
    } else if (!mxIsDouble(prhs[0]) ||
            !mxIsNumeric(prhs[0]) ||
            mxIsSparse(prhs[0]) ||
            mxIsCell(prhs[0]) ||
            mxIsStruct(prhs[0])) {
        
        mexErrMsgTxt("Error: logsumexp1 argument must be double array or scalar");
        return 1;
    } else if (mxIsEmpty(prhs[0])) {
        plhs[0] = mxCreateDoubleScalar(NAN);
        return 1;
    }
    return 0;
}

void logsumexp(int nlhs,
        mxArray *plhs[],
        int nrhs,
        const mxArray *prhs[]) {

    size_t n = mxGetNumberOfElements(prhs[0]);
    double *x = mxGetPr(prhs[0]);
    double max = x[0];
    double y = 0.0;
    size_t i;
    for (i = 1; i < n; i++) {
        if (x[i] > max) {
            max = x[i];
        }
    }
    for (i = 0; i < n; i++) {
        y += exp(x[i] - max);
    }
    plhs[0] = mxCreateDoubleScalar(log(y) + max);
}

void mexFunction(int nlhs,
        mxArray *plhs[],
        int nrhs,
        const mxArray *prhs[]) {
    
    /*if (validate(nlhs, plhs, nrhs, prhs) != 0) {
        return;
    }*/
    logsumexp(nlhs, plhs, nrhs, prhs);
}