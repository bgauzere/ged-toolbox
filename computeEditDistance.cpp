#include <mex.h>
#include <iostream>
#include <cassert>
#include <cstring>
#include "graph.h"
#include "GraphEditDistance.h"

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray*prhs[]) {
  
  /*Graphs loading
    prhs[0] -> G1*/
  int n1 = (int)mxGetDimensions(prhs[0])[0]; 
  double * tmp = mxGetPr(prhs[0]);
  int * am1 = new int[n1*n1];
  for (int i=0;i<n1*n1;i++)
    am1[i] = (int)tmp[i];
  Graph * g1 = new Graph(am1,n1);
  
#if DEBUG
  mexPrintf("Noeud de G1 \n");
  for (int i = 0; i<n1;i++)
    mexPrintf("%d ,", (*g1)[i]->attr);
  mexPrintf("\n");
#endif

  /*prhs[1] -> G2*/
  int n2 = (int)mxGetDimensions(prhs[1])[0];
  tmp = mxGetPr(prhs[1]);
  int * am2 = new int[n2*n2];
  for (int i=0;i<n2*n2;i++)
    am2[i] = (int)(tmp[i]);
  Graph * g2 = new Graph(am2,n2);

#if DEBUG
  mexPrintf("Noeud de G2 \n");
  for (int i = 0; i<n2;i++)
    mexPrintf("%d ,", (*g2)[i]->attr);
  mexPrintf("\n");
#endif
  
  /*mapping  loading
    prhs[2] -> A(i) −> mapping of node  i of g1 onto g2 */
  int dim = (int)mxGetDimensions(prhs[2])[0];
  int * mapping = (int*)mxGetData(prhs[2]); //Take care with the given mapping : munkres give a double * and hungarianLSAP an int *
#if DEBUG
  mexPrintf("Appariemment (%d) \n", dim);
  for(int j=0;j<dim;j++)
    mexPrintf("%d -> %d \n", j,mapping[j]-1);
#endif// int dim = (int)mxGetDimensions(prhs[2])[0];
//   double * mapping = mxGetPr(prhs[2]);
//   for(int j=0;j<dim;j++){
//     mexPrintf("data lue : %d \n", (int)mapping[j]);
//   }
// #if DEBUG
//   mexPrintf("MF :Appariemment (%d) \n", dim);
//   for(int j=0;j<dim;j++)
//     mexPrintf("%d -> %d \n", j,(int)mapping[j]-1);
// #endif

  /*Cost loading*/
  int cns = mxGetScalar(prhs[3]);
  int cni = mxGetScalar(prhs[4]);
  int ces = mxGetScalar(prhs[5]);
  int cei = mxGetScalar(prhs[6]);
  EditDistanceCost * cf = new EditDistanceCost(cns,cni,ces,cei);
  
  
  plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  double cost = GraphEditDistance::GedFromMapping(g1,g2,cf, mapping, dim);
  *mxGetPr(plhs[0])=  cost;
  delete cf;
  delete g2; delete[] am2; delete[] am1; delete g1;
  
  
}
