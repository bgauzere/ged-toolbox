// ======================================================
// ======================================================
#include <mex.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) 
{
  if (nrhs != 3 || nlhs < 1) mexErrMsgTxt("USAGE: [Lkmin,...,Lk] = bagOfLabelSequences(W,n,k)");
  //------------------------------------------------------------------
  // 1rst input : W : adjaceny matrix with (W)_{ii}=label of the node and (W)_{i,j}=label of the edge
  double *W = (double*)mxGetPr(prhs[0]);
  // get dimensions
  int nrows = (int)mxGetDimensions(prhs[0])[0];
  //------------------------------------------------------------------
  // 2nd input : n = index of the node, starting paths
  int n = (int)mxGetScalar(prhs[1])-1;
  //------------------------------------------------------------------
  // 3rd input : k = length of paths from node n
  short k = (short)mxGetScalar(prhs[2]);
  //------------------------------------------------------------------
  int sk = k+1, j, lab;
  std::vector<std::vector<int>* > bos;
  bos.push_back(new std::vector<int>);
  bos[0]->push_back(n);
  std::stack<int> st;
  st.push(0);
  const short skk = sk+1;
  int nbP[skk]; // nbP[i] is the number of paths with i nodes
  std::memset(&nbP[0],0,skk*sizeof(int));
  nbP[1] = 1;
  //------------------------------------------------------------------
  while (!st.empty())
  {
    std::vector<int> *seq = bos[st.top()];
    st.pop();
    n = (*seq)[seq->size()-1];
    if ((int) seq->size()  == sk) continue; // do not update if the depth is reached
    for (j = 0; j < nrows; j++)
    {
      if (j == n) continue; // avoid diagonal
      lab = (int)W[nrows*n+j];
      if (lab == 0) continue; // avoid no edge
      if (std::find(seq->begin(),seq->end(),j) != seq->end()) continue; // avoid passing through a node in the path
      std::vector<int> *vn = new std::vector<int>(*seq);
      bos.push_back(vn);
      vn->push_back(j);
      nbP[vn->size()]++;
      st.push(bos.size()-1);
    }
  }
  //------------------------------------------------------------------
  // outputs
  int dimr[2], kc, beg = std::max(sk-nlhs+1,1);
  int nbPC[skk];
  std::memset(&nbPC[0],0,skk*sizeof(int));
  for (j = beg; j < skk; j++)
  {
    dimr[0] = 2*j-1; dimr[1] = nbP[j];
    plhs[j-beg] = mxCreateNumericArray(2, dimr, mxDOUBLE_CLASS, mxREAL);
  }
  for (j = 0; j < (int)bos.size(); j++) // for each path
  {
    if ((int)bos[j]->size() < beg) // avoid paths having less than 2 nodes
    {
      delete bos[j]; bos[j] = NULL;
      continue;
    }
    double *Pj = (double*)mxGetPr(plhs[bos[j]->size()-beg]);
    kc = nbPC[bos[j]->size()];
    nbPC[bos[j]->size()]++;
    for (n = 0; n < (int)bos[j]->size()-1; n++) // for each node
    {
      Pj[(2*bos[j]->size()-1)*kc+2*n] = W[(*bos[j])[n]*nrows+(*bos[j])[n]];
      Pj[(2*bos[j]->size()-1)*kc+2*n+1] = W[(*bos[j])[n+1]*nrows+(*bos[j])[n]];
    }
    if ((int)bos[j]->size() > 0) Pj[(2*bos[j]->size()-1)*kc+2*n] = W[(*bos[j])[n]*nrows+(*bos[j])[n]];
    delete bos[j]; bos[j] = NULL;
  }
  //------------------------------------------------------------------
}
