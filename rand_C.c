#include "mex.h" 
void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {
    /* Now generate 5 pseudo-random numbers */
    int i;
    for (i=0; i<1000; i++)
    {
      mexPrintf ("x(end+1)= %f;",
        i, (double)rand()/RAND_MAX);
    }
    return;
  }