/** -----------------------------------------------------------
 * Matlab interface to the:
 * Hungarian algorithm for solving the Linear Sum Assignment
 * Problem with Edition (LSAPE)
 * authors: Sebastien Bougleux and Luc Brun
 * institution: Normandie Universite
 *              GREYC UMR 6072
 *              CNRS - Universite Caen Normandie - ENSICAEN
 * reference: S. Bougleux and L. Brun
 *            Linear Sum Assignment Problem with Edition
 *            Research Report, Normandie Universite, GREYC UMR 6072, March 2016
 *
 *      [rho,varrho,u,v]  = hungarianLSAPE(C);
 *
 * Compute an assignment with edition between two sets
 * U and V of size n and m respectively, provided 
 * a (n+1)x(m+1) edit cost matrix C, and such that the
 * total cost sum of the assignment is miniminal.
 *
 * The resulting assignment is provided as two mappings
 * rho:U->V\cup{epsilon} and varrho:V->U\cup{epsilon}.
 *
 * execute matlab file 'compile_mex' to compile this function
   -----------------------------------------------------------
   
   This file is part of LSAPE.
   
   LSAPE is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   LSAPE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.
   
   You should have received a copy of the GNU General Public License,
   and a copy of the GNU Lesser General Public License, along with 
   LSAPE. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mex.h>
#include <cstdio>
#include <string>
#include <hungarian-lsape.hh>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) 
{
  if (nrhs != 1 || nlhs < 3 || nlhs > 4) mexErrMsgTxt("USAGE: [rho,varrho,u,v] = hungarianLSAPE(C)\nUSAGE: [X,u,v] = hungarianLSAPE(C)");
  //------------------------------------------------------------------
  // 1st input : edit cost matrix
  const double *C = mxGetPr(prhs[0]);
  // get dimensions
  int nrows = (int)mxGetDimensions(prhs[0])[0];
  int ncols = (int)mxGetDimensions(prhs[0])[1];
  int nr1 = nrows-1, nc1 = ncols-1;
  //------------------------------------------------------------------
  if (nlhs == 4)
  {
    int dimr[2] = { nr1, 1 }, dimc[2] = { 1, nc1 };
    plhs[0] = mxCreateNumericArray(2, dimr, mxINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, dimc, mxINT32_CLASS, mxREAL);
    int *rho = (int*)mxGetData(plhs[0]);
    int *varrho = (int*)mxGetData(plhs[1]);
    //------------------------------------------------------------------
    plhs[2] = mxCreateNumericArray(2, dimr, mxDOUBLE_CLASS, mxREAL);
    double *u = mxGetPr(plhs[2]);
    plhs[3] = mxCreateNumericArray(2, dimc, mxDOUBLE_CLASS, mxREAL);
    double *v = mxGetPr(plhs[3]);
    hungarianLSAPE<int,double>(C,nrows,ncols,rho,varrho,u,v);
  }
  else if (nlhs == 3)
  {
    int dimx[2] = { nrows, ncols }, dimr[2] = { nr1, 1 }, dimc[2] = { 1, nc1 };
    plhs[0] = mxCreateNumericArray(2, dimx, mxINT32_CLASS, mxREAL);
    int *X = (int*)mxGetData(plhs[0]);
    int *rho = new int[nrows-1], *varrho = new int[ncols-1];
    //------------------------------------------------------------------
    plhs[1] = mxCreateNumericArray(2, dimr, mxDOUBLE_CLASS, mxREAL);
    double *u = mxGetPr(plhs[1]);
    plhs[2] = mxCreateNumericArray(2, dimc, mxDOUBLE_CLASS, mxREAL);
    double *v = mxGetPr(plhs[2]);
    hungarianLSAPE<int,double>(C,nrows,ncols,rho,varrho,u,v);
    reconstructMtx<int,int>(rho,nr1,varrho,nc1,X);
    delete[] rho; delete[] varrho;
  }
}
