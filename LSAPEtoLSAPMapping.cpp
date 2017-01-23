#include <mex.h>
#include <iostream>
#include <cassert>
#include <cstring>

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray*prhs[]) {
  
  int n = (int)mxGetDimensions(prhs[0])[0];
  int m = (int)mxGetDimensions(prhs[1])[0];
  int * G1toG2 = (int*)mxGetData(prhs[0]); //Take care with the given mapping : munkres give a double * and hungarianLSAP an int *
  int * G2toG1 = (int*)mxGetData(prhs[1]); //Take care with the given mapping : munkres give a double * and hungarianLSAP an int *
  plhs[0]=mxCreateDoubleMatrix(n+m,1,mxREAL);
  double * mapping = mxGetPr(plhs[0]);
  int index_eps_G2 = m+1;
  for(int i=0;i<n;i++)
    if(G1toG2[i] >= m){
      mapping[i] = index_eps_G2;
      index_eps_G2 = index_eps_G2 +1;
    }
    else
      mapping[i] = G1toG2[i]+1;
  
  
  int index_eps_G1 = n;
  for (int j=0;j<m;j++)
    if(G2toG1[j] >= n){
      mapping[index_eps_G1] = j + 1;
      index_eps_G1 = index_eps_G1 + 1;
    }
  for(int k = index_eps_G1;k<n+m;k++)
    mapping[k] = k+1;
}
