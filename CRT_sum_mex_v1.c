
#include <mex.h>
/*L = CRT_sum_mex(x,r);*/


void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    mwSize Lenx;
    mwIndex i, j, ij;
    double *x, *RND, *Lsum, *prob;
    double maxx, r;
    Lenx = mxGetM(prhs[0])*mxGetN(prhs[0]);
    x = mxGetPr(prhs[0]);
    r = mxGetScalar(prhs[1]);
    
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    Lsum = mxGetPr(plhs[0]);
    
    
    for(i=0, maxx=0;i<Lenx;i++)
        if (maxx<x[i]) maxx = x[i];
    
    prob = (double *) mxCalloc(maxx, sizeof(double));
    
    for(i=0;i<maxx;i++)
        prob[i] = r/(r+i);
    
    for(ij=0, i=0, Lsum[0]=0;i<Lenx;i++)
        for(j=0;j<x[i];j++) {
            /*if ( ((double) randomMT() / (double) 4294967296.0) <prob[j])     Lsum[0]++;*/
            /*if  ((double) randomMT() <= prob[j]*RAND_MAX_32)     Lsum[0]++;*/
            if ( (double) rand() <= prob[j]*RAND_MAX)     Lsum[0]++;
        }
}
