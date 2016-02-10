// ======================================================
// ======================================================
#include <mex.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) 
{
  if (nrhs != 2 || nlhs < 1 || nlhs > 2) mexErrMsgTxt("USAGE: [Wx,Idxab] = labeledKron(Wa, Wb)");
  //------------------------------------------------------------------
  // 1rst input : Wa : adjaceny matrix with (Wa)_{ii}=label of the node and (Wa)_{i,j}=label of the edge
  double *Wa = (double*)mxGetPr(prhs[0]);
  // get dimensions
  int anrows = (int)mxGetDimensions(prhs[0])[0];
  int ancols = (int)mxGetDimensions(prhs[0])[1];
  //------------------------------------------------------------------
  // 2nd input : source points
  // 1rst input : Wb : adjaceny matrix with (Wb)_{ii}=label of the node and (Wb)_{i,j}=label of the edge
  double *Wb = (double*)mxGetPr(prhs[1]);
  // get dimensions
  int bnrows = (int)mxGetDimensions(prhs[1])[0];
  int bncols = (int)mxGetDimensions(prhs[1])[1];
  //------------------------------------------------------------------
  int i, j, ia, ib, k, l, ja, jb;
  std::vector<std::pair<int,int> > mab;
  //------------------------------------------------------------------
  // 1st output : labeled Kronecker product
  int dimr[2];
  if (nlhs == 2)
  {
    for (i = 0; i < anrows; i++) // for each node of Wa
    {
      ia = anrows*i + i;
      for (j = 0; j < bnrows; j++) // for each node of Wb
      {
          ib = bnrows*j+j;
          if (Wa[ia] == Wb[ib]) // same label
          {    mab.push_back(std::make_pair(i,j)); }
      }
    }
    dimr[0] = mab.size(); dimr[1] = mab.size();
  }
  else { dimr[0] = anrows*bnrows; dimr[1] = ancols*bncols; }
  plhs[0] = mxCreateNumericArray(2, dimr, mxDOUBLE_CLASS, mxREAL);
  double *Wx = (double*)mxGetPr(plhs[0]);
  //------------------------------------------------------------------
  // 2nd output : indicies correspondances between original nodes and product node
  if (nlhs == 2)
  {
    int dimi[2];
    dimi[0] = 2; dimi[1] = mab.size();
    plhs[1] = mxCreateNumericArray(2, dimi, mxINT32_CLASS, mxREAL);
    int *Ixab = (int*)mxGetPr(plhs[1]);
    //------------------------------------------------------------------
    for (i = 0; i < (int)mab.size(); i++) // for each node of Wx
    {
      Ixab[2*i] = mab[i].first+1;
      Ixab[2*i+1] = mab[i].second+1;
      Wx[mab.size()*i + i] = 0;
      for (j = i+1; j < (int)mab.size(); j++) // for each node of Wx
      {
          ia = anrows*mab[i].first + mab[j].first;
          ib = bnrows*mab[i].second + mab[j].second;
	  if (Wa[ia] == Wb[ib]) // same edge label
              Wx[mab.size()*j + i] = Wx[mab.size()*i + j] = Wa[ia];
          else Wx[mab.size()*j + i] = Wx[mab.size()*i + j] = 0;
      }
    }
  }
  else // nlhs == 1
  {
    for (i = 0; i < anrows; i++) // for each node of Wa
    {
      ia = anrows*i + i;
      for (j = 0; j < ancols; j++)
      {
	ja = anrows*j + j;
	for (k = 0; k < bnrows; k++) // for each node of Wb
	{
          ib = bnrows*k+k;
	  for (l = 0; l < bncols; l++)
	  {
	    jb = bnrows*l+l;
	    if (Wa[ia] == Wb[ib] && Wa[ja] == Wb[jb]) // same node label
	    {
	      if (Wa[anrows*j+i] == Wb[bnrows*l+k] && Wa[anrows*j+i] > 0 && Wb[bnrows*l+k] > 0) // same edge label
		Wx[anrows*bnrows*(j*bncols+l)+(i*bnrows+k)] = Wa[anrows*j+i];
	      else Wx[anrows*bnrows*(j*bncols+l)+(i*bnrows+k)] = 0.0;
	    }
	    else Wx[anrows*bnrows*(j*bncols+l)+(i*bnrows+k)] = 0.0;
	  }
	}
      }
    }
  }
}
