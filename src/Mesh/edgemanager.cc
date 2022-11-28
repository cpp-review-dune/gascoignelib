/**
 *
 * Copyright (C) 2004, 2008, 2011 by the Gascoigne 3D authors
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

#include "edgemanager.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "../Common/vecalgo.h"

#include "hangsort.h"

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne {
EdgeManager::EdgeManager(vector<Edge>& e,
                         vector<Quad>& q,
                         const IndexVector& con,
                         IndexVector& eon)
  : edges(e)
  , quads(q)
  , co2n(con)
  , eo2n(eon)
  , QuadLaO(q)
{}

/*---------------------------------------------------*/

void
EdgeManager::BSETest() const
{
  cout << "BSE Tester:\n";
  IndexVector x(edges.size());
  for (IndexType i = 0; i < quads.size(); i++) {
    for (IndexType e = 0; e < 4; e++) {
      IndexType edge = quads[i].edge(e);
      x[edge]++;
    }
  }
  for (IndexType i = 0; i < x.size(); i++) {
    if (x[i] > 2) {
      cout << "problem 1 in edge " << i << " " << x[i] << endl;
    }
  }
  for (IndexType i = 0; i < quads.size(); i++) {
    for (IndexType e = 0; e < 4; e++) {
      IndexType edge = quads[i].edge(e);

      if (edge < 0)
        continue;
      const Edge& E = edges[edge];
      IndexType m = E.master();
      IndexType s = E.slave();
      IndexType ml = E.LocalMasterIndex();
      IndexType sl = E.LocalSlaveIndex();
      IndexType nachbar = -10;
      IndexType f = -1;
      if (m == i) {
        nachbar = s;
        f = sl;
      } else if (s == i) {
        nachbar = m;
        f = ml;
      }

      if (nachbar == -10) {
        cout << "problem 2 in edge " << i << " " << e << " edgenumber " << edge
             << " f= " << f << endl;
      }
    }
  }
}

/*---------------------------------------------------*/

void
EdgeManager::LoadEdgeElimination(IndexVector& edel,
                                 const IndexSet& CellCoarseList,
                                 const HangContainer2d& hangset) const
{
  edel.resize(4 * CellCoarseList.size() + 2 * hangset.NToBeDeleted());

  IndexType n = 0;

  IndexSet::const_iterator cp;

  for (cp = CellCoarseList.begin(); cp != CellCoarseList.end(); cp++) {
    if (*cp < 0)
      continue;

    for (IndexType e = 0; e < 4; e++) {
      edel[n++] = QuadLaO.GlobalInnerEdge(*cp, e);
    }
  }
  HangList<2>::const_iterator p = hangset.Deleting().begin();

  for (; p != hangset.Deleting().end(); p++) {
    for (IndexType i = 0; i < 2; i++) {
      edel[n++] = QuadLaO.GlobalChildEdge(p->first, p->second.rneighbour(), i);
    }
  }
  edel.resize(n);
}

/*---------------------------------------------------*/

void
EdgeManager::Build(const IndexSet& CellRefList, HangContainer2d& hangset)
{
  //    IndexVector SwappedEdge;

  Update();
  InnerEdges(CellRefList);
  // BSETest();
  OuterEdges(hangset);
  OldHangings(hangset, CellRefList);
  SwappedEdges();
  NeighbourTester();
  SortHangings();
}

/*---------------------------------------------------*/

void
EdgeManager::Update()
{
  for (IndexType i = 0; i < edges.size(); i++) {
    IndexType m = edges[i].master();
    IndexType s = edges[i].slave();
    if (m < 0) {
      throw std::runtime_error("EdgeManager::update(): Master negativ (i,m)" +
                               std::to_string(i) + " " + std::to_string(m));
    }
    IndexType nm = co2n[m];
    if (nm >= 0) {
      edges[i].master() = nm;
      if (s >= 0) {
        edges[i].slave() = co2n[s];
      }
    } else {
      if (s < 0) {
        edges[i].master() = -1;
        SwappedEdge.push_back(i);
        continue;
      }
      edges[i].swapping(co2n[s]);
    }
  }
}

/*---------------------------------------------------*/

void
EdgeManager::InnerEdges(const IndexSet& CellRefList)
{
  IndexType n = edges.size();
  IndexType nv1 = CellRefList.size();

  edges.reserve(n + 4 * nv1);
  edges.resize(n + 4 * nv1);

  IndexSet::const_iterator cp;

  for (cp = CellRefList.begin(); cp != CellRefList.end(); cp++) {
    IndexType f = co2n[*cp];

    for (IndexType e = 0; e < 4; e++) {
      IndexType ic = quad(f).child(e);
      IndexType ie = QuadLaO.InnerEdgeOfChild(e, 0);
      quad(ic).edge(ie) = n;

      IndexType oe = QuadLaO.OuterEdgeOfChild(e, 0);
      quads[ic].edge(oe) = -1;
      oe = QuadLaO.OuterEdgeOfChild(e, 1);
      quad(ic).edge(oe) = -1;

      if (ic < 0) {
        cout << "master = " << ic << endl;
      }

      Edge E(ic, ie);

      ic = quad(f).child((e + 1) % 4);
      ie = QuadLaO.InnerEdgeOfChild((e + 1) % 4, 1);

      quad(ic).edge(ie) = n;

      E.slave() = ic;
      E.LocalSlaveIndex() = ie;

      edges[n] = E;
      n++;
    }
  }
}

/*---------------------------------------------------*/

void
EdgeManager::OuterEdges(const HangContainer2d& hangset)
{
  IndexType n = edges.size();
  IndexType nn = n + 2 * hangset.NToBeCreated();

  edges.reserve(nn);
  edges.resize(nn);

  HangList<2>::const_iterator p = hangset.Creating().begin();

  for (; p != hangset.Creating().end(); p++) {
    for (IndexType i = 0; i < 2; i++) {
      EdgeVector edge;
      IndexType hang = p->second.hanging();
      IndexType rneigh = p->second.rneighbour();
      pair<IndexType, IndexType> cp =
        QuadLaO.GetChildEdges(edge, p->first, hang, rneigh, i);

      IndexType cellindex = cp.first;
      IndexType edgeindex = cp.second;

      Edge E(cellindex, edgeindex);

      quads[cellindex].edge(edgeindex) = n;

      edges[n] = E;

      IndexType bigslave = p->second.cneighbour();

      if (bigslave >= 0) {
        if (quads[bigslave].sleep()) {
          IndexType slave = 0;
          for (IndexType j = 0; j < 4; j++) {
            slave = quads[bigslave].child(j);
            edgeindex = QuadLaO.local_edge_index(slave, edge);
            if (edgeindex >= 0)
              break;
          }
          if (edgeindex < 0) {
            throw std::runtime_error(std::to_string(slave) + " ###smallind 2 " +
                                     std::to_string(edgeindex));
          }
          quads[slave].edge(edgeindex) = n;
          edges[n].slave() = slave;
          edges[n].LocalSlaveIndex() = edgeindex;
        }
      }
      n++;
    }
  }
}

/*---------------------------------------------------*/

void
EdgeManager::OldHangings(HangContainer2d& hangset, const IndexSet& CellRefList)
{
  HangList<2>::iterator p = hangset.NotAnyMore().begin();

  for (; p != hangset.NotAnyMore().end(); p++) {
    if (p->second.hanging() < 0)
      continue;

    int bigslave = p->second.cneighbour();
    int bigmaster = p->second.rneighbour();
    int hang = p->second.hanging();

    if (hang < 0)
      continue;
    if (bigslave < 0)
      continue;
    if (bigmaster < 0)
      continue;

    if (!quads[bigslave].sleep()) {
      throw std::runtime_error("BIGSLAVE sleeps");
    }
    for (IndexType i = 0; i < 2; i++) {
      EdgeVector edge;

      IndexType rneigh = p->second.rneighbour();
      pair<IndexType, IndexType> cp =
        QuadLaO.GetChildEdges(edge, p->first, hang, rneigh, i);

      IndexType gedge = quads[cp.first].edge(cp.second);
      if (gedge < 0) {
        swap(bigmaster, bigslave);
        swap(p->second.rneighbour(), p->second.cneighbour());

        rneigh = p->second.rneighbour();
        cp = QuadLaO.GetChildEdges(edge, p->first, hang, rneigh, i);

        gedge = quads[cp.first].edge(cp.second);
        if (gedge < 0) {
          cout << "###gedge negativ:" << bigmaster << " " << bigslave << endl;
        }
      }
      IndexType ledge = -1;
      IndexType slave = -1;
      IndexType master = cp.first;

      for (IndexType j = 0; j < 4; j++) {
        slave = quads[bigslave].child(j);
        ledge = QuadLaO.local_edge_index(slave, edge);
        if (ledge >= 0 && ledge != -1)
          break;
      }
      if (ledge < 0) {
        throw std::runtime_error(std::to_string(slave) + " ###ledge " +
                                 std::to_string(ledge));
      }
      Edge& E = edges[gedge];

      if (E.master() != master) {
        if (E.slave() == master) {
          E.master() = master;
          E.LocalMasterIndex() = E.LocalSlaveIndex();
        } else {
          throw std::runtime_error("bad master" + std::to_string(master) +
                                   ": " + std::to_string(E.master()) + "-" +
                                   std::to_string(E.slave()));
        }
      }
      quads[slave].edge(ledge) = gedge;
      E.slave() = slave;
      E.LocalSlaveIndex() = ledge;
    }
  }
}

/*---------------------------------------------------*/

void
EdgeManager::SwappedEdges()
{
  IndexType n = 0;
  IndexType m = 0;
  for (IndexType i = 0; i < quads.size(); i++) {
    const Quad& q = quads[i];
    if (q.sleep())
      continue;
    for (IndexType e = 0; e < 4; e++) {
      if (q.edge(e) < 0)
        m++;
    }
  }
  if (m != SwappedEdge.size()) {
    throw std::out_of_range(
      "ERROR: Inconsistency in SwappedEdges2d: " + std::to_string(m) + " " +
      std::to_string(SwappedEdge.size()));
  }
  for (IndexType i = 0; i < quads.size(); i++) {
    Quad& q = quads[i];
    if (q.sleep())
      continue;
    for (IndexType e = 0; e < 4; e++) {
      if (q.edge(e) < 0) {
        IndexType ei = SwappedEdge[n++];
        q.edge(e) = ei;

        edges[ei].setmaster(i, e);
      }
    }
  }
}

/*---------------------------------------------------*/

std::array<IndexType, 2>
EdgeManager::ChildrenOfEdge(IndexType e) const
{
  IndexType s = edges[e].slave();
  IndexType is = edges[e].LocalSlaveIndex();

  if (s < 0) {
    throw std::runtime_error("EdgeManager::ChildrenOfEdge(): no slave");
  }
  std::array<IndexType, 2> f;
  for (IndexType ii = 0; ii < 2; ii++) {
    IndexType ic = quad(s).child(QuadLaO.ChildsOfEdge(is, ii));
    f[ii] = quad(ic).edge(QuadLaO.ChildEdge(is));
  }
  return f;
}

/*---------------------------------------------------*/

void
EdgeManager::DeleteEdges()
{
  compress(edges, eo2n);

  for (IndexType i = 0; i < quads.size(); i++) {
    Quad& Q = quads[i];
    if (co2n[i] < 0) {
      Q.edges().fill(-1);
    } else {
      for (IndexType e = 0; e < 4; e++) {
        IndexType ne = eo2n[Q.edge(e)];
        if (ne < 0) {
          throw std::runtime_error("eo2n: " + std::to_string(ne));
        }
        Q.edge(e) = ne;
      }
    }
  }
}

/*---------------------------------------------------*/

bool
EdgeManager::EdgeIsHanging(IndexType e) const
{
  IndexType m = edges[e].master();
  IndexType s = edges[e].slave();
  if (s < 0)
    return 0;
  if (quad(m).sleep() && !quad(s).sleep())
    return 1;
  if (quad(s).sleep() && !quad(m).sleep())
    return 1;
  return 0;
}

/*---------------------------------------------------*/

bool
EdgeManager::EdgeIsHanging(const Edge& e) const
{
  IndexType m = e.master();
  IndexType s = e.slave();
  if (s < 0)
    return 0;
  if (quad(m).sleep() && !quad(s).sleep())
    return 1;
  if (quad(s).sleep() && !quad(m).sleep())
    return 1;
  return 0;
}

/*---------------------------------------------------*/

void
EdgeManager::SortHangings()
{
  // edges with hanging nodes swapped to the end of list

  vector<IndexType> perm(edges.size());
  iota(perm.begin(), perm.end(), 0);
  stable_sort(perm.begin(), perm.end(), HangEdgeSort(*this));
  stable_sort(edges.begin(), edges.end(), HangEdgeSort2(*this));

  IndexType i = edges.size();
  while (EdgeIsHanging(i - 1))
    i--;

  vector<IndexType> permi(perm.size());
  for (IndexType i = 0; i < perm.size(); i++)
    permi[perm[i]] = i;

  for (IndexType i = 0; i < quads.size(); i++) {
    for (IndexType ii = 0; ii < 4; ii++)
      quads[i].edge(ii) = permi[quad(i).edge(ii)];
  }
  for (IndexType i = 0; i < eo2n.size(); i++) {
    if (eo2n[i] >= 0)
      eo2n[i] = permi[eo2n[i]];
  }

  // master of each edge is allways the coarser quad

  for (IndexType i = 0; i < edges.size(); i++) {
    IndexType m = edges[i].master();
    IndexType s = edges[i].slave();
    if (s < 0)
      continue;
    if (quad(m).sleep() && !quad(s).sleep()) {
      swap(edges[i].master(), edges[i].slave());
      swap(edges[i].LocalMasterIndex(), edges[i].LocalSlaveIndex());
    }
  }
}

/*---------------------------------------------------*/

void
EdgeManager::InitEdges()
{
  std::unordered_map<EdgeVector, IndexType, EdgeHash<2>> H;

  EdgeVector e;
  for (IndexType i = 0; i < quads.size(); i++) {
    for (IndexType j = 0; j < 4; j++) {
      e[0] = quads[i].vertex(j);
      e[1] = quads[i].vertex((j + 1) % 4);

      if (e[1] < e[0])
        swap(e[0], e[1]);

      auto yes = H.find(e);

      if (yes != H.end()) {
        IndexType k = yes->second;
        edges[k].slave() = i;
        edges[k].LocalSlaveIndex() = j;
        quads[i].edge(j) = k;
        H.erase(yes);
      } else {
        Edge E(i, j);
        IndexType n = edges.size();
        edges.push_back(E);
        H.insert(make_pair(e, n));
        quads[i].edge(j) = n;
      }
    }
  }
}

/*---------------------------------------------------*/

void
EdgeManager::NeighbourTester() const
{
  IndexType n = quads.size();
  vector<std::array<IndexType, 4>> vecino(n);

  for (IndexType q = 0; q < quads.size(); q++) {
    for (IndexType e = 0; e < 4; e++) {
      IndexType ind = quads[q].edge(e);
      if (ind < 0) {
        cout << "* Q=" << q << " edge=" << e << " " << ind << endl;
      }
    }
  }
  for (IndexType q = 0; q < quads.size(); q++) {
    for (IndexType e = 0; e < 4; e++) {
      IndexType ind = quads[q].edge(e);
      IndexType m = edges[ind].master();
      IndexType s = edges[ind].slave();
      if (m == q) {
        vecino[q][e] = s;
      } else if (s == q) {
        vecino[q][e] = m;
      } else {
        vecino[q][e] = -2;
        cout << "*q=" << q << " el=" << e << " eg=" << ind << ": m=" << m
             << " s=" << s << endl;
      }
    }
  }
}
} // namespace Gascoigne

/*------------------------------------------------------*/
