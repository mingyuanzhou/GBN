
#include <mex.h>
#define MAX(a,b) ((a) > (b) ? a : b)
/* L = CRT_sum_mex_matrix(X,r);   
 X is a K*N sparse matrix, r is a 1*N vector, L is a 1*N vector */
        

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    mwSize Lenx;
    mwIndex i, j, ij, k;
    double *ZSDS, *RND, *Lsum, *prob;
    double maxx, *r_i;    
    
    double  *pr;
    mwIndex *ir, *jc;
    mwIndex Vsize, Nsize, Ksize;
    mwIndex starting_row_index, stopping_row_index, current_row_index;
    
    
    pr = mxGetPr(prhs[0]);
    ir = mxGetIr(prhs[0]);
    jc = mxGetJc(prhs[0]); 
    Ksize = mxGetM(prhs[0]);
    Nsize = mxGetN(prhs[0]);
    
    
    
    r_i = mxGetPr(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(1, Nsize, mxREAL);
    Lsum = mxGetPr(plhs[0]);
    
    for (j=0;j<Nsize;j++) {
        starting_row_index = jc[j];
        stopping_row_index = jc[j+1];
        if (starting_row_index == stopping_row_index)
            continue;
        else {
            for (current_row_index =  starting_row_index; current_row_index<stopping_row_index; current_row_index++) {
                maxx = MAX(maxx,pr[current_row_index]);
            }
            prob = (double *) mxCalloc(maxx, sizeof(double));
            for(i=0;i<maxx;i++)
                prob[i] = r_i[j]/(r_i[j]+i);
            
            for (Lsum[j]=0, current_row_index =  starting_row_index; current_row_index<stopping_row_index; current_row_index++)
                /*k = ir[current_row_index];*/   
                for(i=0;i<pr[current_row_index];i++) {
                    if  ((double) random() <= prob[i]*RAND_MAX)     
                        Lsum[j]++;
                }
        }
    }
}
