/**
 *
 * Copyright (C) 2004, 2006, 2007, 2008 by the Gascoigne 3D authors
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

#include "levelmesh2d.h"

#include <algorithm>

#include "../Common/set2vec.h"
#include "../Interface/gascoignehash.h"

#include "leveljumper.h"
#include "levelsorter.h"

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne {
LevelMesh2d::LevelMesh2d(const HierarchicalMesh* hmp)
  : Index()
{
  HMP = dynamic_cast<const HierarchicalMesh2d*>(hmp);
}

/*---------------------------------------------------*/

LevelMesh2d::~LevelMesh2d() {}

/*---------------------------------------------------*/

bool
LevelMesh2d::EdgeIsHangingGlobalIndex(IndexType i) const
{
  IndexType igm = HMP->edge(i).master();
  IndexType igs = HMP->edge(i).slave();

  // rand oder kleien hedges
  if (igs < 0)
    return 0;

  bool m_in_lmesh = (Quadg2l().find(igm) != Quadg2l().end());
  bool s_in_lmesh = (Quadg2l().find(igs) != Quadg2l().end());

  if (m_in_lmesh && s_in_lmesh)
    return 0;

  IndexType ivg = HMP->NodeOnEdge(i);
  if (Vertexg2l().find(ivg) == Vertexg2l().end())
    return 0;

  return 1;
}

/*---------------------------------------------------*/

void
LevelMesh2d::BasicInit(const IndexSet& newq, const IndexSet& oldq)
{
  // doch sortiert

  IndexType n = newq.size() + oldq.size();
  Index::Quadl2g().memory(n);
  IndexVector::const_iterator p = set_union(newq.begin(),
                                            newq.end(),
                                            oldq.begin(),
                                            oldq.end(),
                                            Index::Quadl2g().begin());
  n = p - Index::Quadl2g().begin();

  InitCells(n);
  InitNodes(n);
  //  InitEdges(n);
}

/*-----------------------------------------*/

void
LevelMesh2d::InitCells(IndexType n)
{
  Index::Quadl2g().memory(n);

  sort(Index::Quadl2g().begin(), Index::Quadl2g().end(), LevelSorter2d(*HMP));

  Index::InitQuads();
}

/*-----------------------------------------*/

void
LevelMesh2d::InitNodes(IndexType n)
{
  IndexSet nodes;
  for (IndexType i = 0; i < n; i++) {
    IndexType ind = Index::Quadl2g()[i];
    for (IndexType ii = 0; ii < 4; ii++) {
      nodes.insert(HMP->vertex_of_cell(ind, ii));
    }
  }
  Index::InitNodes(nodes);
}

/*-----------------------------------------*/

void
LevelMesh2d::InitEdges(IndexType n)
{
  // edges
  IndexSet edges;
  for (IndexType i = 0; i < n; i++) {
    const Quad& Q = HMP->quad(Index::Quadl2g()[i]);
    for (IndexType ii = 0; ii < 4; ii++) {
      edges.insert(Q.edge(ii));
    }
  }
  Index::InitEdges(edges);

  // sortiere haengende edges nach hinten

  stable_sort(Edgel2g().begin(), Edgel2g().end(), HangEdgeSort3(*this));

  Edgeg2l().clear();
  for (IndexType i = 0; i < Edgel2g().size(); i++) {
    Edgeg2l().insert(make_pair(Edgel2g()[i], i));
  }
}

/*-----------------------------------------*/

bool
LevelMesh2d::BuildFathers(IndexSet& Vaeter) const
{
  for (IndexType i = 0; i < ncells(); i++) {
    const Quad& q = quad(i);
    IndexType findex = q.father();
    if (findex == -1) {
      return 0;
    }

    const Quad& qf = HMP->quad(findex);
    for (IndexType ii = 0; ii < qf.nchilds(); ii++) {
      IndexType cindex = qf.child(ii);
      if (Quadg2lCheck(cindex) == -2) {
        return 0;
      }
    }
    Vaeter.insert(findex);
  }
  return 1;
}

/*-----------------------------------------*/

bool
LevelMesh2d::ConstructCellIndOfPatch(IndexVector& dst) const
{
  IndexSet Vaeter;
  BuildFathers(Vaeter);

  IndexType nq = ncells() / 4;
  dst.reservesize(nq);

  IndexType count = 0;
  IndexSet::const_iterator pf = Vaeter.begin();
  while (pf != Vaeter.end()) {
    IndexType findex = *pf;
    dst[count] = findex;

    count++;
    pf++;
  }
  return 1;
}

/*-----------------------------------------*/

void
LevelMesh2d::ConstructIndOfPatch(nvector<IndexVector>& dst) const
{
  IndexSet Vaeter;
  BuildFathers(Vaeter);

  IndexType nq = ncells() / 4;
  dst.reserve(nq);
  dst.resize(nq, IndexVector(9));

  IndexType count = 0;
  IndexSet::const_iterator pf = Vaeter.begin();
  while (pf != Vaeter.end()) {
    IndexType findex = *pf;
    const Quad& qf = HMP->quad(findex);

    std::array<IndexType, 4> FineQuads;
    for (IndexType ii = 0; ii < qf.nchilds(); ii++) {
      FineQuads[ii] = qf.child(ii);
    }
    dst[count][0] = Vertexg2l(HMP->quad(FineQuads[0]).vertex(0));
    dst[count][1] = Vertexg2l(HMP->quad(FineQuads[0]).vertex(1));
    dst[count][2] = Vertexg2l(HMP->quad(FineQuads[1]).vertex(1));
    dst[count][3] = Vertexg2l(HMP->quad(FineQuads[0]).vertex(3));
    dst[count][4] = Vertexg2l(HMP->quad(FineQuads[0]).vertex(2));
    dst[count][5] = Vertexg2l(HMP->quad(FineQuads[1]).vertex(2));
    dst[count][6] = Vertexg2l(HMP->quad(FineQuads[3]).vertex(3));
    dst[count][7] = Vertexg2l(HMP->quad(FineQuads[3]).vertex(2));
    dst[count][8] = Vertexg2l(HMP->quad(FineQuads[2]).vertex(2));

    count++;
    pf++;
  }
}

/*---------------------------------------------------*/

void
LevelMesh2d::ConstructHangingStructureQuadratic(
  QuadraticHNStructure3& hnq2) const
{
  hnq2.clear();
  IndexSet habschon;
  for (IndexType i = 0; i < ncells(); i++) {
    const Quad& q = HMP->quad(Quadl2g(i));
    IndexType father = q.father();
    if (father < 0)
      continue;
    if (habschon.find(father) != habschon.end())
      continue;
    habschon.insert(father);
    for (IndexType in = 0; in < 4; in++) {
      IndexType neighbour = HMP->neighbour(father, in);
      if (neighbour < 0)
        continue;
      const Quad& qfn = HMP->quad(neighbour);
      if (qfn.nchilds() == 0)
        continue;

      std::array<IndexType, 2> childs;
      IndexType ne = in;

      {
        IndexType start = 0;
        IndexType neighbourneighbour = HMP->neighbour(neighbour, ne);
        while ((neighbourneighbour != father) && (start < 4)) {
          start++;
          ne = (ne + 1) % 4;
          neighbourneighbour = HMP->neighbour(neighbour, ne);
        }

        assert(neighbourneighbour == father);
      }
      HMP->QuadLawOrder().childs_of_edge(childs, qfn, ne);
      IndexType child = childs[0];
      if (Quadg2lCheck(child) >= 0)
        continue;
      const Quad& qfc = HMP->quad(child);

      if (qfc.nchilds() == 0)
        continue;

      IndexType enkel = qfc.child(0);
      if (Quadg2lCheck(enkel) < 0)
        continue;

      // jetzt haengt er
      IndexType hn = Vertexg2l(HMP->QuadLawOrder().edge_vertex(qfc, ne));
      std::array<IndexType, 3> F;
      F[0] = qfn[ne];
      F[1] = HMP->QuadLawOrder().edge_vertex(qfn, ne);
      F[2] = qfn[(ne + 1) % 4];

      assert((qfc[0] == F[0]) || (qfc[1] == F[0]) || (qfc[2] == F[0]) ||
             (qfc[3] == F[0]));

      for (IndexType k = 0; k < 3; k++)
        F[k] = Vertexg2l(F[k]);

      hnq2.insert(make_pair(hn, F));

      const Quad& qfc2 = HMP->quad(childs[1]);
      hn = Vertexg2l(HMP->QuadLawOrder().edge_vertex(qfc2, ne));

      swap(F[0], F[2]);

      hnq2.insert(make_pair(hn, F));
    }
  }
}

/*---------------------------------------------------*/

void
LevelMesh2d::ConstructHangingStructureQuartic(QuarticHNStructure5& hnq4) const
{
  hnq4.clear();
  assert(HMP->patchdepth() >= 2);
  IndexType count = 0;
  IndexSet habschon;
  for (IndexType i = 0; i < ncells(); ++i) {
    const Quad& q_c = HMP->quad(Quadl2g(i));
    IndexType father_c = q_c.father();
    assert(father_c >= 0);
    const Quad& f_c = HMP->quad(father_c);
    IndexType opa_c = f_c.father();
    assert(opa_c >= 0);

    // wurde der opa schon behandelt?
    if (habschon.find(opa_c) != habschon.end())
      continue;
    habschon.insert(opa_c);

    // alle Nachbarn
    for (IndexType in = 0; in < 4; ++in) {
      IndexType opa_r = HMP->neighbour(opa_c, in);
      // Nachbar existiert nicht
      if (opa_r < 0)
        continue;
      const Quad& o_r = HMP->quad(opa_r);
      // muss einmal verfeinert sein
      assert(o_r.nchilds() != 0);

      // nachbar-nachbar finden
      IndexType ne = in;
      {
        IndexType start = 0;
        IndexType neighbourneighbour = HMP->neighbour(opa_r, ne);
        while ((neighbourneighbour != opa_c) && (start < 4)) {
          start++;
          ne = (ne + 1) % 4;
          neighbourneighbour = HMP->neighbour(opa_r, ne);
        }
        assert(neighbourneighbour == opa_c);
      }
      // beide Vaeter entlang der gemeinsamen Kante
      std::array<IndexType, 2> fathers_r;
      HMP->QuadLawOrder().childs_of_edge(fathers_r, o_r, ne);

      // wenn gleiches Level, dann haengt nix.
      const Quad& f0_r = HMP->quad(fathers_r[0]);
      const Quad& f1_r = HMP->quad(fathers_r[1]);
      // groeber?
      if (f0_r.nchilds() == 0)
        continue;
      // gleiches level?
      const Quad& q0_r = HMP->quad(f0_r.child(0));
      if (q0_r.nchilds() == 0)
        continue;

      // der Enkel muss aktiv sein.
      assert(Quadg2lCheck(q0_r.child(0)) >= 0);

      // Jetzt haengt ne Menge

      // die vier enkel verfeinert
      vector<std::array<IndexType, 2>> enkels_r(2);
      HMP->QuadLawOrder().childs_of_edge(enkels_r[0], f0_r, ne);
      HMP->QuadLawOrder().childs_of_edge(enkels_r[1], f1_r, ne);
      assert(enkels_r[0][0] >= 0);
      assert(enkels_r[0][1] >= 0);
      assert(enkels_r[1][0] >= 0);
      assert(enkels_r[1][1] >= 0);

      // die 5 Knoten
      std::array<IndexType, 6> tmp;
      tmp[0] = HMP->quad(enkels_r[0][0])[ne];
      tmp[1] = HMP->quad(enkels_r[0][1])[ne];
      tmp[2] = HMP->quad(enkels_r[1][0])[ne];
      tmp[3] = HMP->quad(enkels_r[1][1])[ne];
      tmp[4] = HMP->quad(enkels_r[1][1])[(ne + 1) % 4];

      // die haengen
      std::array<IndexType, 4> hn;
      hn[0] = HMP->QuadLawOrder().edge_vertex(HMP->quad(enkels_r[0][0]), ne);
      hn[1] = HMP->QuadLawOrder().edge_vertex(HMP->quad(enkels_r[0][1]), ne);
      hn[2] = HMP->QuadLawOrder().edge_vertex(HMP->quad(enkels_r[1][0]), ne);
      hn[3] = HMP->QuadLawOrder().edge_vertex(HMP->quad(enkels_r[1][1]), ne);

      tmp[5] = 0;
      hnq4.insert(std::make_pair(hn[0], tmp));
      tmp[5] = 1;
      hnq4.insert(std::make_pair(hn[1], tmp));
      std::swap(tmp[0], tmp[4]);
      std::swap(tmp[1], tmp[3]);
      tmp[5] = 1;
      hnq4.insert(std::make_pair(hn[2], tmp));
      tmp[5] = 0;
      hnq4.insert(std::make_pair(hn[3], tmp));

      ++count;
    }
  }
}

/*---------------------------------------------------*/

void
LevelMesh2d::check_leveljump() const
{
  LevelJumper Phi;
  for (IndexType c = 0; c < ncells(); c++) {
    Phi.update(quad(c));
  }
  assert(!Phi.check());
  // if(Phi.check()) cerr << "LevelMesh2d::check_leveljump() aenderb\n";
}

/*---------------------------------------------------*/

void
LevelMesh2d::fill_opis(IndexSet& dst, IndexSet& oldquads) const
{
  for (IndexType i = 0; i < ncells(); i++) {
    const Quad& Q = quad(i);

    IndexType f = Q.father();
    assert(f >= 0);

    IndexType opi = HMP->quad(f).father();

    if (opi < 0) {
      IndexType j = Quadl2g(i);
      oldquads.insert(j);
    } else {
      dst.insert(opi);
    }
  }
}

/*---------------------------------------------------*/

void
LevelMesh2d::fill_childs(IndexSet& dst, const Quad& Q) const
{
  for (IndexType ii = 0; ii < Q.nchilds(); ii++) {
    IndexType qccindex = Q.child(ii);
    dst.insert(qccindex);
  }
}

/*---------------------------------------------------*/

void
LevelMesh2d::fill_enkel(IndexSet& oldquads, const Quad& Q) const
{
  for (IndexType ii = 0; ii < Q.nchilds(); ii++) {
    IndexType qccindex = Q.child(ii);
    const Quad& qcc = HMP->quad(qccindex);
    for (IndexType iii = 0; iii < qcc.nchilds(); iii++) {
      IndexType qcindex = qcc.child(iii);
      if (Quadg2lCheck(qcindex) >= 0) {
        oldquads.insert(qcindex);
      }
    }
  }
}

/*---------------------------------------------------*/

bool
LevelMesh2d::EnkelUniform(const Quad& Q) const
{
  for (IndexType ii = 0; ii < Q.nchilds(); ii++) {
    IndexType qccindex = Q.child(ii);
    const Quad& qcc = HMP->quad(qccindex);
    for (IndexType iii = 0; iii < qcc.nchilds(); iii++) {
      IndexType qcindex = qcc.child(iii);
      if (Quadg2lCheck(qcindex) == -2) {
        return 0;
      }
    }
  }
  return 1;
}

/*---------------------------------------------------*/

void
LevelMesh2d::construct_lists(IndexSet& newquads, IndexSet& oldquads) const
{
  newquads.clear();
  oldquads.clear();

  check_leveljump();

  IndexSet Opis;
  fill_opis(Opis, oldquads);
  for (IndexSet::const_iterator p = Opis.begin(); p != Opis.end(); p++) {
    const Quad& Q = HMP->quad(*p);

    if (EnkelUniform(Q)) {
      fill_childs(newquads, Q);
    } else {
      fill_enkel(oldquads, Q);
    }
  }

  // Iteration zum Regulaer machen (LevelJump)
  while (1) {
    LevelJumper Phi;
    IndexSet::const_iterator p;
    for (p = newquads.begin(); p != newquads.end(); p++) {
      Phi.update(HMP->quad(*p));
    }
    for (p = oldquads.begin(); p != oldquads.end(); p++) {
      Phi.update(HMP->quad(*p));
    }
    if (!Phi.check())
      break;

    // IndexType rep=0;
    IndexSet help(newquads);
    for (p = newquads.begin(); p != newquads.end(); p++) {
      const Quad& q = HMP->quad(*p);
      if (!Phi.VertexOK(q)) {
        // rep++;
        const Quad& qf = HMP->quad(q.father());
        for (IndexType ii = 0; ii < 4; ii++) {
          help.erase(qf.child(ii));
        }
        fill_enkel(oldquads, qf);
      }
    }
    newquads = help;
    // cerr << "\t Regular Iteration\t" << count++ << " " << rep << endl;
  }
}

/*---------------------------------------------------*/

void
LevelMesh2d::InitBoundaryHandler(BoundaryIndexHandler& BI,
                                 const PatchIndexHandler& PIH) const
{
  IndexSet blines;
  for (IndexType i = 0; i < HMP->nblines(); i++) {
    IndexType q = HMP->bline(i).of_quad();
    if (Quadg2l().find(q) != Quadg2l().end()) {
      blines.insert(i);
    }
  }

  // which colors are there ?
  BI.clear();
  for (IndexSet::const_iterator p = blines.begin(); p != blines.end(); p++) {
    const BoundaryLine& bl = HMP->bline(*p);
    IndexType col = bl.material();

    BI.GetColors().insert(col);
  }
  IndexVector colorvec;
  Set2Vec(colorvec, BI.GetColors());

  // compute inverse positions
  map<IndexType, IndexType> inv;

  for (IndexType i = 0; i < colorvec.size(); i++) {
    inv.insert(make_pair(colorvec[i], i));
  }

  IndexType nc = colorvec.size();
  vector<IndexSet> H1(nc); // for verteces
  // for cells and local indices
  vector<set<std::array<IndexType, 2>>> H2(nc);

  for (IndexSet::const_iterator q = blines.begin(); q != blines.end(); q++) {
    const BoundaryLine& bl = HMP->bline(*q);
    IndexType col = bl.material();

    map<IndexType, IndexType>::const_iterator p = inv.find(col);
    if (p == inv.end()) {
      cout << "LevelMesh2d::BuildVertexOnBoundary()" << endl;
      abort();
    }
    IndexType pos = p->second;

    for (IndexType ii = 0; ii < 2; ii++) {
      IndexType vindex = Vertexg2l(bl.vertex(ii));
      H1[pos].insert(vindex);
    }
    std::array<IndexType, 2> ind;
    ind[0] = Quadg2l(bl.of_quad());
    ind[1] = bl.edge_in_quad();
    H2[pos].insert(ind);
  }
  BI.CopySetToVector(H1, colorvec, BI.GetVertex());

  for (IndexType i = 0; i < H2.size(); i++) {
    IndexVector v1(H2[i].size());
    IndexVector v2(H2[i].size());
    IndexType j = 0;

    set<std::array<IndexType, 2>>::const_iterator p;
    for (p = H2[i].begin(); p != H2[i].end(); p++) {
      v1[j] = (*p)[0];
      v2[j] = (*p)[1];
      j++;
    }
    IndexType color = colorvec[i];

    BI.GetCell().insert(make_pair(color, v1));
    BI.GetLocal().insert(make_pair(color, v2));
  }

  const nvector<IndexVector>& patch2cell = PIH.GetAllPatch2Cell();

  nvector<IndexType> cell2patch(PIH.npatches() << 2);
  for (IndexType p = 0; p < patch2cell.size(); ++p)
    for (IndexType i = 0; i < patch2cell[p].size(); ++i)
      cell2patch[patch2cell[p][i]] = p;

  for (IndexSet::const_iterator c = BI.GetColors().begin();
       c != BI.GetColors().end();
       c++) {
    IndexType col = *c;
    const IndexVector& cells = BI.Cells(col);
    const IndexVector& locals = BI.Localind(col);
    std::unordered_set<IndexType> habschon;

    IndexVector p1;
    IndexVector p2;

    for (IndexType i = 0; i < cells.size(); i++) {
      IndexType iq = cells[i];
      IndexType ip = cell2patch[iq];
      IndexType ile = locals[i];

      // gabs den patch schon
      if (habschon.find((ip << 2) + ile) != habschon.end())
        continue;
      habschon.insert((ip << 2) + ile);
      p1.push_back(ip);
      p2.push_back(ile);
    }
    BI.GetPatch().insert(make_pair(col, p1));
    BI.GetLocalPatch().insert(make_pair(col, p2));
  }
}
} // namespace Gascoigne
#undef HASHSET
/*--------------------------------------------------------------*/
