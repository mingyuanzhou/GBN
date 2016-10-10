
#include <mex.h>
#include <math.h>
/*BesselCount  = Rand_Truncated_PMF_bessel(PMF,Mode,alpha)*/


void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    mwSize Msize, Nsize;
    mwIndex i, j;
    double *PMF, *Mode, *alpha;
    double *Count, prob;
    mwIndex u;
    
    Msize = mxGetM(prhs[0]);
    Nsize = mxGetN(prhs[0]);
    PMF = mxGetPr(prhs[0]);
    Mode = mxGetPr(prhs[1]);
    alpha = mxGetPr(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(Msize,Nsize,mxREAL);
    
    Count = mxGetPr(plhs[0]);
    
    for(i=0;i<Msize*Nsize;i++){
       
        for (j=0;PMF[i]>0;j++){
            u = (mwIndex) Mode[i] + j;
            prob = exp(u*(2*log(alpha[i])-log(4))-lgamma(u)-lgamma(u+1) -alpha[i]);
            if (prob/PMF[i] >= (double) rand() / RAND_MAX){
                Count[i]=u;
                break;
            }
            else{
                PMF[i] -= prob;
            }
            if(j>0 & j<(mwIndex) Mode[i]){
                u = (mwIndex) Mode[i] - j;
                prob = exp(u*(2*log(alpha[i])-log(4))-lgamma(u)-lgamma(u+1) -alpha[i]);
                if (prob/PMF[i] >= (double) rand() / RAND_MAX){
                    Count[i]=u;
                    break;
                }
                else{
                    PMF[i] -= prob;
                }
            }
        }
        if (PMF[i]<=0){
            mexWarnMsgTxt("warning, use mode");
            Count[i]=(mwIndex) Mode[i];
        }
            
    }
      
}
