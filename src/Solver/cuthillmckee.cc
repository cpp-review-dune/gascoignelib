/**
 *
 * Copyright (C) 2004, 2008, 2009, 2011 by the Gascoigne 3D authors
 *
 * This file is part of Gascoigne 3D
 *
 * Gascoigne 3D is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version.
 *
 * Gascoigne 3D is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * Please refer to the file LICENSE.TXT for further information
 * on this license.
 *
 **/

#include "cuthillmckee.h"
#include "ilupermutate.h"
#include "metis.h"

using namespace std;

/* --------------------------------------------------------------- */

namespace Gascoigne {
// extern "C" void METIS_NodeND
// (int*, int*, int*, int*, int*, int*, int*);

/* --------------------------------------------------------------- */

CuthillMcKee::CuthillMcKee(const StencilInterface* s)
{
  S = 0;
  DS = 0;
  CS = 0;
  CS = dynamic_cast<const ColumnStencil*>(s);
  DS = dynamic_cast<const DynamicStencil*>(s);
  S = s;
  assert(CS || DS);
  assert(S);
}

/* --------------------------------------------------------------- */

CuthillMcKee::CuthillMcKee()
{
  //   M=0;
  DS = 0;
  CS = 0;
  //   dimension=0;
}

// /* --------------------------------------------------------------- */

// void CuthillMcKee::Permutate    (IndexVector &perm, const Vertex2d v)
// {
//   dimension=2;
// //   dir2d=v;
//   Permutate(perm);
// }

// /* --------------------------------------------------------------- */

// void CuthillMcKee::Permutate    (IndexVector &perm, const Vertex3d v)
// {
//   dimension=3;
// //   dir3d=v;
//   Permutate(perm);
// }

/* --------------------------------------------------------------- */

void
CuthillMcKee::Permutate(IndexVector& perm)
{
  // mit metis graph aufbauen
  //
  // diese methode wird in StdSolver::PermutateIlu aufgerufen
  // unmittelbar davor *wird* perm schon richtig gesized und mit iota gefuellt
  // dieser code:
  //  IndexType n = S->n();
  //  perm.resize(n);
  // sollte hier nicht aufgerufen werden, sondern ueberpruefen ob perm von der
  // groesse her stimmt
  idx_t n = S->n();
  assert(n == perm.size());

  vector<idx_t> adj(n + 1, 0);
  vector<idx_t> adjncy;

  IndexType count = 0;
  IndexType c;
  adj[0] = 0;
  assert(CS || DS);
  if (CS)
    for (IndexType r = 0; r < n; ++r) {
      for (IndexType p = CS->start(r); p < CS->stop(r); ++p) {
        c = CS->col(p);
        if (r == c)
          continue;
        ++count;
        adjncy.push_back(c);
      }
      adj[r + 1] = count;
    }
  else if (DS)
    for (IndexType r = 0; r < n; ++r) {
      for (list<IndexType>::const_iterator it = DS->cstart(r);
           it != DS->cstop(r);
           ++it) {
        c = *it;
        if (r == c)
          continue;
        ++count;
        adjncy.push_back(c);
      }
      adj[r + 1] = count;
    }

  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 0;

  vector<idx_t> iperm(n);
  vector<idx_t> Mperm(n);

  for (IndexType i = 0; i < n; ++i)
    Mperm[i] = perm[i];
  METIS_NodeND(
    &n, &adj[0], &adjncy[0], NULL, &options[0], &Mperm[0], &iperm[0]);
  for (IndexType i = 0; i < n; ++i)
    perm[i] = Mperm[i];

  //    assert((dimension==0)||(M->dimension()==dimension));
  //      				   // Liste mit Nachbarn aufbauen
  //    typedef SparseStructure::const_iterator SetIterator;
  //    IndexType n = S->n();
  //    neighbors.resize(n);
  //    for (unsigned IndexType i=0;i<n;++i)
  //      {
  //        neighbors[i]=S->rowsize(i);
  //      }

  //  				   // Startvektor finden
  //  //  stable_sort(perm.begin(),perm.end(),*this);
  //    IndexType starting_node=0;
  //    for (IndexType i=1;i<perm.size();++i)
  //      if (operator()(perm[i],perm[starting_node])) starting_node=i;

  //  				   // jetzt werden die Knoten so sortiert
  //  				   // dass immer die Nachbarn eines Knotens
  //  				   // kommen, sortiert nach neighbors
  //    nvector<bool> versorgt(n,false);
  //    IndexType index = 0;
  //    perm[index]=starting_node; versorgt[perm[index]]=true;
  //    IndexType st=0;
  //    IndexType en=0;
  //    IndexType st1=0;
  //    while (en<n-1)
  //      {
  //        st1=index+1;
  //        if (en<st) abort();
  //        for (IndexType i=st;i<=en;++i)
  //  	{
  //    	  for (IndexType p=S->start(perm[i]);
  //  	       p<S->stop(perm[i]);++p)
  //    	    if (!versorgt[S->col(p)])
  //    	      {
  //    		++index; perm[index]=S->col(p);
  //  		versorgt[perm[index]]=true;
  //    	      }
  //    	  if (index>st1)
  //  	    stable_sort(perm.begin()+st1,perm.begin()+index+1,*this);
  //  	}
  //        st=en+1;
  //        en=index;
  //      }

  //    vector<IndexType> iperm(n);
  //    for (IndexType i=0;i<n;++i) iperm[perm[i]]=i;

  //    char s[10];
  //    sprintf(s,"l_%d",n);
  //    ofstream aus(s);
  //    for (IndexType r=0;r<n;++r)
  //      {
  //        IndexType row = perm[r];
  //        for (IndexType p = S->start(row);p!=S->stop(row);++p)
  //  	aus << r << "\t" << iperm[S->col(p)] << endl;
  //      }

  //    aus.close();
}

// /* --------------------------------------------------------------- */

// bool CuthillMcKee::operator()(IndexType i,IndexType j) const
// {
//   if (neighbors[i]<neighbors[j]) return 1;
//   if (neighbors[i]>neighbors[j]) return 0;
//   if (dimension==2) return (((M->vertex2d(j)-M->vertex2d(i))*dir2d)>0);
//   if (dimension==3) return (((M->vertex3d(j)-M->vertex3d(i))*dir3d)>0);
//   return 0;
// }

#ifdef __WITH_THREADS__
void
CuthillMcKee::Permutate(
  IndexVector& perm,
  const IndexVector& nodes_in_domain,
  const std::vector<std::vector<std::pair<IndexType, IndexType>>>& node2domain,
  IndexType domain)
{
  // mit metis graph aufbauen
  //
  // diese methode wird in StdSolver::PermutateIlu aufgerufen
  // unmittelbar davor *wird* perm schon richtig gesized und mit iota gefuellt
  // dieser code:
  //  IndexType n = S->n();
  //  perm.resize(n);
  // sollte hier nicht aufgerufen werden, sondern ueberpruefen ob perm von der
  // groesse her stimmt
  IndexType n = nodes_in_domain.size();
  assert(n == perm.size());

  vector<IndexType> adj(n + 1, 0);
  vector<IndexType> adjncy;

  IndexType count = 0;
  IndexType c;
  IndexType globalnodeid;

  IndexType ignoredcouplings = 0;

  adj[0] = 0;
  if (CS)
    for (IndexType r = 0; r < n; ++r) {
      globalnodeid = nodes_in_domain[r];
      for (IndexType p = CS->start(globalnodeid); p < CS->stop(globalnodeid);
           ++p) {
        c = CS->col(p);
        if (globalnodeid == c)
          continue;
        // Falls  c nicht in der Subdomain liegt auch ignorieren
        IndexType localposition = -1;
        for (IndexType l = 0; l < node2domain[c].size(); l++) {
          if (node2domain[c][l].first == domain) {
            localposition = node2domain[c][l].second;
            break;
          }
        }
        if (localposition == -1) {
          ignoredcouplings++;
          continue;
        }
        // Ende der Pruefung ob c in der Subdomain liegt

        ++count;
        //	      adjncy.push_back(c);
        // use local indices to avoid metis errors, when nodes have larger
        // indices then n
        c = localposition;
        adjncy.push_back(c);
      }
      adj[r + 1] = count;
    }
  //  else if (DS)
  //      for (IndexType r=0;r<n;++r)
  //      {
  //	  for (list<IndexType>::const_iterator it =
  // DS->cstart(r);it!=DS->cstop(r);++it)
  //	  {
  //	      c = *it;
  //	      if (r==c) continue;
  //	      ++count;
  //	      adjncy.push_back(c);
  //	  }
  //	  adj[r+1]=count;
  //      }
  else
    abort();

  //  std::cout<<"Ignored Couplings: "<<ignoredcouplings<<" Used  Couplings:
  //  "<<adjncy.size()<<std::endl;

  IndexType numflag = 0;
  IndexType options[8];
  options[0] = 1;
  options[1] = 3;
  options[2] = 1;
  options[3] = 1;
  options[4] = 0;
  options[5] = 1;
  // options[5]=2;
  options[6] = 0;
  options[7] = 1;
  vector<IndexType> iperm(n);

  METIS_NodeND(
    &n, &adj[0], &adjncy[0], &numflag, &options[0], &perm[0], &iperm[0]);
  // perm ist nun der Vector bei dem perm[i] sagt, das der Knoten
  // nodes_of_domain[i] auf nodes_of_domain_[perm[i]] abgebilded  wird.
}
// End of Thread-Version
#endif
} // namespace Gascoigne
