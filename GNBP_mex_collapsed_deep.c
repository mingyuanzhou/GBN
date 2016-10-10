#include "mex.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
/*///#include "cokus.c"
 * //#define RAND_MAX_32 4294967295.0*/



mwIndex BinarySearch(double probrnd, double *prob_cumsum, mwSize Ksize) {
    mwIndex k, kstart, kend;
    if (probrnd <=prob_cumsum[0])
        return(0);
    else {
        for (kstart=1, kend=Ksize-1; ; ) {
            if (kstart >= kend) {
                /*//k = kend;*/
                return(kend);
            }
            else {
                k = kstart+ (kend-kstart)/2;
                if (prob_cumsum[k-1]>probrnd && prob_cumsum[k]>probrnd)
                    kend = k-1;
                else if (prob_cumsum[k-1]<probrnd && prob_cumsum[k]<probrnd)
                    kstart = k+1;
                else
                    return(k);
            }
        }
    }
    return(k);
}



void GNBP_mex_collapsed_deep(double *ZSDS, double *ZSWS, double *n_dot_k,
        double *ZS, double *WS, double *DS, mwSize Vsize, mwSize Nsize, mwSize Ksize, mwSize WordNum,
        double *shape, double *eta, double *prob_cumsum)
{
    double cum_sum, probrnd, sumeta;
    mwIndex i, k, j, v, Kvalid;
    /*    for (sumr=0,j=0;j<Nsize;j++){
        sumr += r_k[j];
    }
    for (sumeta=0,v=0;v<Vsize;v++){
        sumeta += eta[v];
    }
    //Kvalid = Ksize;*/
    for (i=0;i<WordNum;i++){
        /*mexPrintf("%d ",i);*/
        v = (mwIndex) WS[i] -1;
        j = (mwIndex) DS[i] -1;
        k = (mwIndex) ZS[i] -1;
        if(ZS[i]>0){
            ZSDS[k+Ksize*j]--;
            ZSWS[k+Ksize*v]--;
            n_dot_k[k]--;
            /*//             if((mwIndex)n_dot_k[k]==0)
//                 Kvalid = Kvalid-1;*/
        }
        
        for (cum_sum=0, k=0; k<Ksize; k++) {
            cum_sum += (eta[0]+ ZSWS[k+Ksize*v])/((double)Vsize *eta[0] + n_dot_k[k])*(ZSDS[k+Ksize*j] + shape[k+Ksize*j]);
            prob_cumsum[k] = cum_sum;
        }
        
        /*
        cum_sum = cum_sum + r_star/((double)Vsize);
        prob_cumsum[Ksize] = cum_sum;
         */
        probrnd = (double)rand()/(double)RAND_MAX*cum_sum;
        /* //probrnd = (double) randomMT()/RAND_MAX_32*cum_sum;*/
        
        k = BinarySearch(probrnd, prob_cumsum, Ksize);
        ZS[i] = k+1;
        ZSDS[k+Ksize*j]++;
        ZSWS[k+Ksize*v]++;
        n_dot_k[k]++;
        /*
        if (k==Ksize){
            return(i+1);
        }
        else{
            ZSDS[k+Ksize*j]++;
            ZSWS[k+Ksize*v]++;
            n_dot_k[k]++;
        }*/
    }
   
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *end_index;
    double *ZSDS, *ZSWS, *n_dot_k, *ZS;
    double *WS, *DS, *shape, *eta, *prob_cumsum;
    mwSize Vsize, Nsize, Ksize, WordNum, Kcurrent;
    /*//mxArray * ZSDSworkspace, *ZSWSworkspace, *sumZSDSworkspace, *sumZSWSworkspace;*/
    mxArray *prob_cumsum_Array;
    
   
    
    WS = mxGetPr(prhs[4]);
    DS = mxGetPr(prhs[5]);
    
    shape = mxGetPr(prhs[6]);
    eta = mxGetPr(prhs[7]);
    
    Ksize = mxGetM(prhs[0]);
    Nsize = mxGetN(prhs[0]);
    Vsize = mxGetN(prhs[1]);
    
    WordNum = mxGetM(prhs[4])*mxGetN(prhs[4]);
    
    plhs[0] = mxDuplicateArray(prhs[0]);
    plhs[1] = mxDuplicateArray(prhs[1]);
    plhs[2] = mxDuplicateArray(prhs[2]);
    plhs[3] = mxDuplicateArray(prhs[3]);
    plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    ZSDS = mxGetPr(plhs[0]);
    ZSWS = mxGetPr(plhs[1]);
    n_dot_k = mxGetPr(plhs[2]);
    ZS = mxGetPr(plhs[3]);
    
    
    prob_cumsum_Array = mxCreateDoubleMatrix(Ksize,1,mxREAL);
    prob_cumsum = mxGetPr(prob_cumsum_Array);
    
    
    GNBP_mex_collapsed_deep(ZSDS, ZSWS, n_dot_k, ZS, WS, DS, Vsize, Nsize, Ksize, WordNum, shape, eta, prob_cumsum);
    
    mxDestroyArray(prob_cumsum_Array);
    
}
