/*
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * % 
 * % quadratic term computes X^T * D
 * %
 * % Compute Linearization Matrix: compute the A, a matrix 2-dimensional
 * % linearization of quadratic cost, based on parameters:
 * %           A1       first matrix of adjacency
 * %           A2       second matrix of adjacency
 * %           cei      cost of edge insertion
 * %           ced      cost of edge deletion
 * %           ces      cost of edge substitution
 * %           M        matrix of permutation of the form [ P Q ]
 * %                                                      [ S T ]
 * %           Indices  matrix 2-col correspond to i -> phi(i) mappings
 * %   Returns
 * %           A :  X'D matrix
 * %
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */
#include <mex.h>
#include <cassert>
#include <cstring>
#define A1(a,b) a1[(a-1)+(b-1)*N]
#define A2(a,b) a2[(a-1)+(b-1)*M]
#define M(a,b) m[(a-1)+(b-1)*(N+M)]
#define I(a,b) indices[(a)+(b*dimIndices)]
// #define PHI_I(i) indices[(i-1)+(N+M)] 
#define A(a,b) Ajl[(a-1)+(b-1)*(N+M)]
//TODO : macros for S and R 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) {
    
  //get matrix A1
  int N = (int)mxGetDimensions(prhs[0])[0];
  double * a1 = mxGetPr(prhs[0]);
  // mexPrintf("QuadraticTerm.cpp : N=%d\n", N);
    
  //get matrix A2
  int M = (int)mxGetDimensions(prhs[1])[0];
  double * a2 = mxGetPr(prhs[1]);
  // mexPrintf("QuadraticTerm.cpp : M=%d\n", M);
    
  //get values of cei, ced, ces
  int cei = mxGetScalar(prhs[2]);
  int ced = mxGetScalar(prhs[3]);
  int ces = mxGetScalar(prhs[4]);
  //mexPrintf("cei=%d\n", cei);
  //mexPrintf("ced=%d\n", ced);
  //mexPrintf("ces=%d\n", ces);
    
  //get matrix P
  int dimMap = (int)mxGetDimensions(prhs[5])[0]; //n+m
  // mexPrintf("QuadraticTerm.cpp : dimMap=%d\n", dimMap);
  assert(dimMap == N+M);
  double * m = mxGetPr(prhs[5]); 
  
  //get mapping information
  int dimIndices = (int)mxGetM(prhs[6]); 
  double * indices = mxGetPr(prhs[6]);
  // mexPrintf("dimIndices=%d\n", dimIndices);
  
  
  plhs[0] = mxCreateDoubleMatrix(dimMap,dimMap,mxREAL);
  double * Ajl = mxGetPr(plhs[0]);
  memset((void*) Ajl,0,sizeof(double)* dimMap*dimMap);
  for(int j = 1; j < dimMap+1; j++){
    for(int l = 1; l < dimMap+1;l++){
	 double sum = 0;
	   // mexPrintf("------- \n");
	 for(int count = 0; count < dimIndices; count++){
	   //Traversal of all couples of j->l
	   int i = I(count,0); // j != count
	   int k = I(count,1);
	   // mexPrintf("(i,k) = (%d,%d)\n",i,k);
	   // mexPrintf("D(%d,%d,%d,%d)\n",i+1,j+1,k+1,l+1);
	   // mexPrintf("%d ->  %d \n",i,k);
	   bool eps_i,eps_j,eps_k,eps_l;
	   eps_i = (i > N);
	   eps_j = (j > N);
	   eps_k = (k > M);
	   eps_l = (l > M);
	
	   int delta_e1 = 0; //delta_e1
	   int delta_e2 = 0; //delta_e2
                
	   //Test delta E1 (i,j)
	   if (i <= N && j <= N) // check for not exceeding matrix
	     if (A1(i,j) != 0)
	       delta_e1 = 1; //(i,j) exists in G1
	   //mexPrintf("eij=%d\n", eij);
                
	   //Test delta E2 (k,l)
	   if (k <= M && l <= M)
	     if (A2(k,l) != 0)
	       delta_e2 = 1;
	   //mexPrintf("ekl=%d\n", ekl);
	
	   float cost = 0;

	   //Computation of D i j k l => too factorize
	   //If (i,j) and (k,l) are both same nodes, no edges between them, so delta_e1 and delta_e2 are both 0, and so the cost
	   if((i != j) && (k != l)){
	     //Substitution case
	  
	     if ( (!eps_i) && (!eps_j) && (!eps_k) && (!eps_l)){
	       if (delta_e1 && delta_e2) // sub
		 cost = ces*(A2(k,l) != A1(i,j));
	       else if (delta_e1 && !delta_e2) //deletion
		 cost = ced;
	       else if (! delta_e1 && delta_e2)
		 cost = cei;
	     }	  
	     //l is epsilon, (k,l) do not exist => deletion of (i,j)
	     else if (i <= N && j <= N && k <= M && l > M){
	       cost = ced*delta_e1;
	     }
	     //k is epsilon => (k,l) do not exist so delete (i,j) if exists
	     else if (i <= N && j <= N && k > M && l <= M){
	       cost = ced*delta_e1;
	     }
	     //(k,l) do not exists, del (i,j) if exists (factorizable with previous case)
	     else if (i <= N && j <= N && k > M && l > M){
	       cost = ced*delta_e1;
	     }
	     //i is epsilon => add (k,l) if exists
	     else if (i > N && j <= N && k <= M && l <= M){
	       cost = cei*delta_e2;
	     }
	     //j is epsilon => add (k,l) if exists
	     else if (i <= N && j > N && k <= M && l <= M){
	       cost = cei*delta_e2;
	     }
	     //i,j are epsilon, (i,j) do not exists, add (k,l) if it exists
	     else if (i > N && j > N && k <= M && l <= M){
	       cost = cei*delta_e2;
	     }
	     // // no (i,j) edges => insert (k,l)
	     // else if (i > N && j > N && k <= M && l <= M){
	     //   cost = cei*delta_e2;
	     // }
// #ifdef DEBUG
	     // if (M(i,k))
	     // 	   mexPrintf("** cost %d,%d %d %d = %f\n",i,j,k,l,cost);
	     // mexPrintf("M(%d,%d)=%f\n",i,k,M(i,k));
// #endif
	   sum += cost * M(i,k);
	   // if((j==6) && (l ==2))
	   //   if (cost * M(i,k))
	   //     mexPrintf("** cost %d,%d %d %d = %f\n",i,j,k,l,cost);
	   // mexPrintf("(%d,%d) -> (%d,%d) =  %f\n",j,l,i,k, cost * M(i,k));   
	 }
	 
	 
      }
	 A(j,l) = sum;
    }
  }
}
