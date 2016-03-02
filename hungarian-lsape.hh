// =========================================================================
/* Hungarian algorithm for solving the Linear Sum Assignment
   Problem with Edition (LSAPE)
   
   authors: Sebastien Bougleux and Luc Brun
   institution: GREYC UMR 6072
                CNRS - Universite Caen Normandie - ENSICAEN
 
   reference: S. Bougleux and L. Brun
              Linear Sum Assignment Problem with Edition
              GREYC Technical Report, Oct. 2015
 
   hungarianLSAPE(C,n+1,m+1,rho,varrho)

   Compute an assignment with edition between two sets
   U and V of size n and m respectively, provided 
   a (n+1)x(m+1) edit cost matrix C, and such that the
   total cost sum of the assignment is miniminal.
   The resulting assignment is provided as two mappings
   rho:U->V\cup{epsilon} and varrho:V->U\cup{epsilon}.
   - Time complexity 0(n^2m+nm^2)
   - Memory space 0(nm)
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
// =========================================================================

#ifndef _HUNGARIAN_LSAPE_
#define _HUNGARIAN_LSAPE_

#include <limits>

// -----------------------------------------------------------
// Compute a initial partial assignment (rho,varrho) and associated dual variables (u,v)
// according to min on rows and then to min on reduced columns
// nass and mass are the number assigned elements in U and V respectively
// -----------------------------------------------------------
template<typename IT, class DT>
void basic_preproc(const DT *C, const IT &nrows, const IT &ncols,
		   IT *rho, IT *varrho, DT *u, DT *v, IT &nass, IT &mass)
{
  const IT n = nrows-1, m = ncols-1;
  IT i = 0, j;
  bool *U = new bool[n], *V = new bool[m];
  DT mn, val;
  nass = mass = 0;
  
  // find the min of each row
  for (; i < n; i++)
  {
    mn = std::numeric_limits<DT>::max();
    for (j = 0; j < ncols; j++)
    {
      const DT &c = C[j*nrows+i];
      if (c < mn) mn = c;
    }
    u[i] = mn;
  }
  
  // find de min of each column
  for (j = 0; j < m; j++)
  {
    mn = std::numeric_limits<DT>::max();
    for (i = 0; i < n; i++)
    {
      val = C[j*nrows+i] - u[i];
      if (val < mn) mn = val;
    }
    // last row
    const DT &c = C[j*nrows+n];
    if (c < mn) mn = c;
    v[j] = mn;
    V[j] = false;
    varrho[j] = -1;
  }
  
  // assign
  for (i = 0; i < n; i++)
  {
    U[i] = false;
    rho[i] = -1;
    for (j = 0; j < m; j++)
    {
      val = C[j*nrows+i] - u[i] - v[j];
      if (U[i] == false && V[j] == false && val == 0)
	{ rho[i] = j; varrho[j] = i; U[i] = true; V[j] = true; nass++; mass++; break; }
    }
    // last column
    if (U[i] == false && C[m*nrows+i]-u[i] == 0) { rho[i] = m; U[i] = true; nass++; }
  }
  
  // last row
  for (j = 0; j < m; j++)
    if (V[j] == false && C[j*nrows+n]-v[j] == 0) { varrho[j] = n; V[j] = true; mass++; }

  delete[] U; delete[] V;
}

static int nb_aug, nb_dual;

// -----------------------------------------------------------
// Construct an alternating tree rooted in an element of V
// and perform dual updates until an augmenting path is found
// -----------------------------------------------------------
template<typename IT, class DT>
void augmentCol(const IT &k, const DT *C, const IT &nrows, const IT &m, const IT *rho, const IT *varrho, 
		IT *U, IT *SV, IT *pred, DT *u, DT *v, DT *pi, IT &zi, IT &zj)
{
  const IT n = nrows-1;
  IT i = 0, j = k, r = 0, *SVptr = SV, *ulutop = U, *uluptr = NULL;
  DT delta = 0, mx = std::numeric_limits<DT>::max(), cred = 0;
  const IT *svptr = NULL, *luptr = NULL, *uend = U+n, *lusutop = U;
  bool lstrw = false;
  *SV = -1;
  zj = zi = -1;
  
  for (i = 0; i < n; i++) { pi[i] = mx; U[i] = i; }
  
  while (true)
  {
    *SVptr = j; *(++SVptr) = -1;
    
    // last row: null element epsilon, it is a sink of an alternating path
    if (varrho[j] < n && C[j*nrows+n] == v[j]) { zi = n; zj = j; return; }
    
    for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
    {
      r = *uluptr;
      cred = C[j*nrows+r] - u[r] - v[j];
      if (cred < pi[r])
      {
	pred[r] = j;
	pi[r] = cred;
	if (cred == 0)
	{
	  if (rho[r] == -1 || rho[r] == m) { zi = r; zj = -1; return; }
	  i = *ulutop; *ulutop = r; *uluptr = i; ++ulutop;
	}
      }
    }
    
    if (lusutop == ulutop) // dual update
    {
      nb_dual++;
      delta = mx; lstrw = false;
      for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
	if (pi[*uluptr] < delta) delta = pi[*uluptr];
      for (svptr = SV; *svptr != -1; ++svptr) // last row
      {
	cred = C[*svptr*nrows+n] - v[*svptr];
	if (cred <= delta) { delta = cred; lstrw = true; zj = *svptr; }
      }
      for (svptr = SV; *svptr != -1; ++svptr) v[*svptr] += delta;
      for (luptr = U; luptr != ulutop; ++luptr) u[*luptr] -= delta;
      if (lstrw) { zi = n; return; }
      for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
      {
	pi[*uluptr] -= delta;
	if (pi[*uluptr] == 0) 
	{
	  if (rho[*uluptr] == -1 || rho[*uluptr] == m) { zi = *uluptr; zj = -1; return; }
	  r = *ulutop; *ulutop = *uluptr; *uluptr = r; ++ulutop;
	}
      }
    } // end dual update
    i = *lusutop; ++lusutop; // i is now in SU
    //if (rho[i] == -1 || rho[i] == m) { zi = i; zj = -1; break; }
    j = rho[i];
  }
}

// -----------------------------------------------------------
// Construct an alternating tree rooted in an element of U
// and perform dual updates until an augmenting path is found
// -----------------------------------------------------------
template<typename IT, class DT>
void augmentRow(const IT &k, const DT *C, const IT &nrows, const IT &m, const IT *rho, const IT *varrho, 
		IT *V, IT *SU, IT *pred, DT *u, DT *v, DT *pi, IT &zi, IT &zj)
{
  const IT n = nrows-1, *suptr = NULL, *lvptr = NULL, *vend = V+m, *lvsvtop = V;
  IT i = k, j = 0, c = 0, *SUptr = SU, *vlvtop = V, *vlvptr = NULL;
  DT delta = 0, mx = std::numeric_limits<DT>::max(), cred = 0;
  bool lstcl = false;
  *SU = -1;
  zj = zi = -1;
  
  for (j = 0; j < m; j++) { pi[j] = mx; V[j] = j; }
  
  while (true)
  {
    *SUptr = i; *(++SUptr) = -1;
    
    // last column: null element epsilon, it is a sink of an alternating path
    if (rho[i] < m && C[m*nrows+i] == u[i]) { zi = i; zj = m; return; }
    
    for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // U\LU
    {
      c = *vlvptr;
      cred = C[c*nrows+i] - u[i] - v[c];
      if (cred < pi[c])
      {
	pred[c] = i;
	pi[c] = cred;
	if (cred == 0)
	{
	  if (varrho[c] == -1 || varrho[c] == n) { zi = -1; zj = c; return; }
	  j = *vlvtop; *vlvtop = c; *vlvptr = j; ++vlvtop;
	}
      }
    }
        
    if (lvsvtop == vlvtop) // dual update
    {
      nb_dual++;
      delta = mx; lstcl = false;
      for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	if (pi[*vlvptr] < delta) delta = pi[*vlvptr];
      for (suptr = SU; *suptr != -1; ++suptr) // last column
      {
	cred = C[m*nrows+*suptr] - u[*suptr];
	if (cred <= delta) { delta = cred; lstcl = true; zi = *suptr; }
      }
      for (suptr = SU; *suptr != -1; ++suptr) u[*suptr] += delta;
      for (lvptr = V; lvptr != vlvtop; ++lvptr) v[*lvptr] -= delta;
      if (lstcl) { zj = m; return; }
      for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	{
	  pi[*vlvptr] -= delta;
	  if (pi[*vlvptr] == 0)
	  {
	    if (varrho[*vlvptr] == -1 || varrho[*vlvptr] == n) { zi = -1; zj = *vlvptr; return; }
	    c = *vlvtop; *vlvtop = *vlvptr; *vlvptr = c; ++vlvtop;
	  }
	}
    } // end dual update
    j = *lvsvtop; ++lvsvtop; // j is now in SV
    //if (varrho[j] == -1 || varrho[j] == n) { zi = -1; zj = j; break; }
    i = varrho[j];
  }
}

// -----------------------------------------------------------
// Main Hungarian algorithm
// -----------------------------------------------------------
template<typename IT, class DT>
void hungarianLSAPE(const DT *C, const IT &nrows, const IT &ncols, IT *rho, IT *varrho, DT *u, DT *v)
{
  const IT n = nrows-1, m = ncols-1;
  IT nass = 0, mass = 0, i = -1, j = -1, r = -1, c = -1, k = -1, *S = NULL, *U = NULL, *V = NULL, *pred = NULL;
  DT *pi = NULL;
  nb_aug = 0; nb_dual = 0;
  
  // initialization -----------------------------------------------
  basic_preproc(C,nrows,ncols,rho,varrho,u,v,nass,mass);

  // augmentation of columns --------------------------------------
  if (mass < m)
  {
    U = new IT[nrows]; S = new IT[ncols];  pi = new DT[n]; pred = new IT[n];
    for (k = 0; k < m; k++)
      if (varrho[k] == -1)
      {
	nb_aug++;
	augmentCol(k,C,nrows,m,rho,varrho,U,S,pred,u,v,pi,i,j);  // augment always finds an augmenting path
	if (i == n) { r = varrho[j]; varrho[j] = i; i = r; }
	else j = -1;
	for (; j != k;)  // update primal solution = new partial assignment
	{
	  j = pred[i]; rho[i] = j;
	  r = varrho[j]; varrho[j] = i; i = r;
	}
      }
    delete[] U; delete[] pred; delete[] pi; delete[] S;
  }
  // augmentation of rows --------------------------------------
  if (nass < n)
  {
    V = new IT[ncols]; S = new IT[nrows];  pi = new DT[m]; pred = new IT[m];
    for (k = 0; k < n; k++)
      if (rho[k] == -1)
      {
	nb_aug++;
	augmentRow(k,C,nrows,m,rho,varrho,V,S,pred,u,v,pi,i,j);  // augment always finds an augmenting path
	if (j == m) { c = rho[i]; rho[i] = j; j = c; }
	else i = -1;
	for (; i != k;)  // update primal solution = new partial assignment
	{
	  i = pred[j]; varrho[j] = i;
	  c = rho[i]; rho[i] = j; j = c;
	}
      }
    delete[] V; delete[] pred; delete[] pi; delete[] S;
  }
}

// -----------------------------------------------------------
// Reconstruct the list of insertions epsilon->V
// return false if no insertion is found
// -----------------------------------------------------------
template <typename IT>
IT reconstructInsertions(const IT *varrho, const IT &n, const IT &m, IT **rhoeps)
{
  IT nb = 0;
  for (IT j = 0; j < m; j++) if (varrho[j] == n) nb++;
  if (nb > 0) *rhoeps = new IT[nb];
  else return nb;
  for (IT j = 0, k = 0; j < m; j++) if (varrho[j] == n) { (*rhoeps)[k] = j; k++; }
  return nb;
}

// -----------------------------------------------------------
#endif
