/*==========================================================
 * Mult_Sparse.c - 
 *
 *
 * The calling syntax is:
 *
 *		Xout = Mult_Sparse(Xmask,Phi,Theta);
 *      Xmask is a sparse matrix
 *      Phi and Theta are full matrices
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2015 Mingyuan Zhou
 *
 *========================================================*/
/* $Revision: 0.1 $ */

#include "mex.h"

void Mult_Matrix(double *Phi, double *Theta, mwIndex *ir, mwIndex *jc, double *pr, mwSize Vsize, mwSize Nsize, mwSize Ksize) 
{    
  
    mwIndex k, j, v, token,  total=0; 
    mwIndex starting_row_index, stopping_row_index, current_row_index;
    double cum_sum;
    
    for (j=0;j<Nsize;j++) {
        starting_row_index = jc[j];
        stopping_row_index = jc[j+1];
        if (starting_row_index == stopping_row_index)
            continue;
        else {
            for (current_row_index =  starting_row_index; current_row_index<stopping_row_index; current_row_index++) {
                v = ir[current_row_index];                
                for (cum_sum=0,k=0; k<Ksize; k++) {
                    /*//prob[k] = Phi(v+ k*Vsize)*Theta(k + Ksize*i);*/
                    cum_sum += Phi[v+ k*Vsize]*Theta[k + Ksize*j];                     
                }
                pr[total++] = cum_sum;
            }
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *Phi, *Theta;
    double  *pr;
    mwIndex *ir, *jc;
    mwIndex Vsize, Nsize, Ksize;

      
     Vsize = mxGetM(prhs[0]);
     Nsize = mxGetN(prhs[0]);
     Ksize = mxGetN(prhs[1]);
     Phi = mxGetPr(prhs[1]);
     Theta = mxGetPr(prhs[2]);

    
     plhs[0] = mxDuplicateArray(prhs[0]);
     pr = mxGetPr(plhs[0]);
     ir = mxGetIr(plhs[0]);
     jc = mxGetJc(plhs[0]);      
 
   
   Mult_Matrix(Phi, Theta, ir, jc, pr,  Vsize, Nsize, Ksize); 
}