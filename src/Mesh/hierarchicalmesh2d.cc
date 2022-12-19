/**
 *
 * Copyright (C) 2004, 2005, 2006, 2011 by the Gascoigne 3D authors
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

#include "hierarchicalmesh2d.h"

#include <fstream>
#include <iterator>
#include <numeric>
#include <set>
#include <stdio.h>

#include "../Common/set2vec.h"
#include "../Common/stlio.h"
#include "../Common/vecalgo.h"

#include "coarsehierarchicalmesh2d.h"
#include "deletecells.h"
#include "edgemanager.h"
#include "levelcomparer2d.h"
#include "regular_update.h"

using namespace std;

/*------------------------------------------------------*/

namespace Gascoigne {
HierarchicalMesh2d::HierarchicalMesh2d()
  : HierarchicalMesh()
  , QuadLaO(quads)
{}

/*------------------------------------------------------*/

HierarchicalMesh2d::HierarchicalMesh2d(const ParamFile& paramfile)
  : HierarchicalMesh()
  , QuadLaO(quads)
{
  BasicInit(paramfile);
}

/*------------------------------------------------------*/

HierarchicalMesh2d::HierarchicalMesh2d(const HierarchicalMesh2d& H)
  : HierarchicalMesh()
  , QuadLaO(quads)
{
  *this = H;
}

/*------------------------------------------------------*/

HierarchicalMesh2d&
HierarchicalMesh2d::operator=(const HierarchicalMesh2d& H)
{
  HierarchicalMesh::operator=(H);
  // copy all data
  vertexs2d = H.vertex2d();
  quads = H.quad();
  Blines = H.bline();
  LineHang = H.linehang();
  // line shapes werden nicht kopiert sonst geht patchrefine in die britze
  quadofcurved = H.GetQuadOfCurved();

  return *this;
}

/*------------------------------------------------------*/

pair<IndexType, IndexType>
HierarchicalMesh2d::GetBoundaryInformation(IndexType i) const
{
  int material = -1;
  int le = -1;
  IndexType ib = GetBoundaryCellOfCurved(i);
  if (ib >= 0) {
    material = bline(ib).material();
    le = bline(ib).edge_in_quad();
  }
  return make_pair(material, le);
}

/*------------------------------------------------------*/

IndexType
HierarchicalMesh2d::FindPatchDepth() const
{
  // simple version, sucht nur p=1, p=0
  for (IndexType i = 0; i < ncells(); i++) {
    const Quad& q = quad(i);
    if (q.sleep())
      continue;
    int father = q.father();
    if (father == -1)
      return 0;
    const Quad& qf = quad(father);
    for (IndexType ii = 0; ii < 4; ii++) {
      if (quad(qf.child(ii)).sleep())
        return 0;
    }
  }
  return 1;
}

/*------------------------------------------------------*/

const BoundaryFunction<2>*
HierarchicalMesh2d::line_shape(IndexType i) const
{
  if (GetCurvedShapes().empty())
    return NULL;
  if (GetCurvedShapes().Curved(i))
    return &(GetCurvedShapes().GetShape(i));
  return NULL;
}

/*------------------------------------------------------*/

set<IndexType>
HierarchicalMesh2d::GetColors() const
{
  set<IndexType> coleur;
  for (IndexType i = 0; i < nblines(); i++) {
    coleur.insert(bline(i).material());
  }
  return coleur;
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::InitQuadOfCurved()
{
  quadofcurved.clear();

  if (GetCurvedShapes().empty())
    return;

  for (IndexType il = 0; il < nblines(); ++il) {
    const BoundaryLine& B = bline(il);

    if (GetCurvedShapes().Curved(B.material())) {
      IndexType iq = B.of_quad();
      quadofcurved.insert(make_pair(iq, il));
    }
  }
}

/*------------------------------------------------------*/

IndexType
HierarchicalMesh2d::Vater(const IndexType i) const
{
  return quad(i).father();
}

/*------------------------------------------------------*/

IndexVector
HierarchicalMesh2d::Nachkommen(const IndexType i) const
{
  IndexVector k = Kinder(i);
  if (k.size() == 0)
    return k;
  IndexVector k1;
  IndexType ks = k.size();
  for (IndexType i = 0; i < ks; ++i) {
    IndexVector k1 = Nachkommen(k[i]);
    for (IndexType j = 0; j < k1.size(); ++j)
      k.push_back(k1[j]);
  }
  return k;
}

/*------------------------------------------------------*/

IndexType
HierarchicalMesh2d::nactivedescendants(IndexType i) const
{
  if (!quads[i].sleep())
    return 1;
  IndexType k = 0;
  for (IndexType j = 0; j < quads[i].nchilds(); ++j)
    k += nactivedescendants(quads[i].child(j));
  return k;
}

/*------------------------------------------------------*/

IndexVector
HierarchicalMesh2d::GetVertices(IndexType c) const
{
  IndexVector v;
  for (IndexType i = 0; i < quads[c].nvertexs(); ++i)
    v.push_back(quads[c][i]);
  return v;
}

/*------------------------------------------------------*/

IndexVector
HierarchicalMesh2d::Kinder(const IndexType i) const
{
  IndexVector k = quad(i).childs();
  return k;
}

/*------------------------------------------------------*/

IndexVector
HierarchicalMesh2d::Geschwister(const IndexType i) const
{
  const Quad& q = quad(i);
  int father = q.father();
  if (father == -1) {
    IndexVector n(1, i);
    return n;
  }
  return Kinder(father);
}

/*------------------------------------------------------*/

std::array<IndexType, 2>
HierarchicalMesh2d::ChildrenOfEdge(IndexType e) const
{
  IndexType s = edge(e).slave();
  IndexType is = edge(e).LocalSlaveIndex();

  assert(s >= 0);

  std::array<IndexType, 2> f;
  for (IndexType ii = 0; ii < 2; ii++) {
    IndexType ic = quad(s).child(QuadLaO.ChildsOfEdge(is, ii));
    f[ii] = quad(ic).edge(QuadLaO.ChildEdge(is));
  }

  return f;
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh2d::NodeOnEdge(IndexType e) const
{
  // only for real hanging nodes
  IndexType m = edge(e).master();
  IndexType im = edge(e).LocalMasterIndex();
  IndexType s = edge(e).slave();
  IndexType is = edge(e).LocalSlaveIndex();

  if (s < 0) {
    assert(quad(im).sleep());

    return QuadLaO.edge_vertex(quad(m), im);
  }

  return QuadLaO.edge_vertex(quad(s), is);
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::ghostglobalcoarse(HangContainer2d& hangset,
                                      const IndexSet& cellcoarse)
{
  EdgeVector lineglob;

  IndexVector coarse;
  Set2Vec(coarse, cellcoarse);
  LevelComparer2d lc(*this, coarse);

  IndexVector BTS(lc.size());
  iota(BTS.begin(), BTS.end(), 0);
  sort(BTS.begin(), BTS.end(), CompareObjectBigToSmall<LevelComparer2d>(lc));

  for (IndexType cp = 0; cp < coarse.size(); cp++) {
    IndexType f = coarse[BTS[cp]];
    for (IndexType edge = 0; edge < 4; ++edge) {
      QuadLaO.global_edge_unsorted(lineglob, quad(f), edge);

      IndexType ve = QuadLaO.edge_vertex(quad(f), edge);
      hangset.ghost_coarse(lineglob, f, ve);
    }
  }
  ghost_fill_neighbours2d();
  hangset.make_consistent();
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::ghost2d(HangContainer2d& hangset,
                            const IndexSet& cellref,
                            const IndexSet& cellcoarse)
{
  EdgeVector lineglob;

  for (IntSetIt cp = cellcoarse.begin(); cp != cellcoarse.end(); ++cp) {
    IndexType f = *cp;
    for (IndexType edge = 0; edge < 4; ++edge) {
      QuadLaO.global_edge_unsorted(lineglob, quad(f), edge);

      IndexType ve = QuadLaO.edge_vertex(quad(f), edge);
      hangset.ghost_coarse(lineglob, f, ve);
    }
  }
  for (IntSetIt cp = cellref.begin(); cp != cellref.end(); ++cp) {
    IndexType f = *cp;
    for (IndexType edge = 0; edge < 4; ++edge) {
      QuadLaO.global_edge_unsorted(lineglob, quad(f), edge);

      hangset.ghost_refine(lineglob, f);
    }
  }
  ghost_fill_neighbours2d();
  hangset.make_consistent();
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::prepare2d(const IndexVector& cell_ref,
                              const IndexVector& cell_coarse,
                              IndexSet& CellRefList,
                              IndexSet& CellCoarseList)
{
  /* copies cell_ref into CellRefList without duplets  */

  IndexVector::const_iterator cp = cell_ref.begin();
  while (cp != cell_ref.end()) {
    IndexType c = *cp;
    if ((c >= 0) && (c < quads.size())) {
      if (!quads[c].sleep()) {
        CellRefList.insert(c);
      }
    }
    ++cp;
  }

  /* copies cell_coarse into CellCoarseList without duplets
     checks if :  -- coarse cell in refine list
                  -- coarse cell is cneighbour of a hang */

  IndexSet help;

  for (IndexVector::const_iterator cp = cell_coarse.begin();
       cp != cell_coarse.end();
       ++cp) {
    IndexType ic = *cp;
    if ((ic < 0) || (ic >= quads.size()))
      continue;

    Quad& q = quads[ic];
    if (q.sleep())
      continue;
    if (!q.level())
      continue;
    if (CellRefList.find(ic) != CellRefList.end())
      continue;

    help.insert(ic);
  }

  LineHangList::const_iterator Lp;

  for (Lp = LineHang.begin(); Lp != LineHang.end(); Lp++) {
    IndexType cn = Lp->second.cneighbour();
    if (help.find(cn) != help.end())
      help.erase(cn);
  }

  multiset<IndexType> ff;

  for (IntSetCIt hp = help.begin(); hp != help.end(); ++hp) {
    ff.insert(quad(*hp).father());
  }

  for (multiset<IndexType>::iterator fp = ff.begin(); fp != ff.end(); ++fp) {
    if (ff.count(*fp) == 4) {
      CellCoarseList.insert(*fp);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::ghost_fill_neighbours2d()
{
  EdgeVector lineglob;
  for (IndexType ic = 0; ic < quads.size(); ic++) {
    for (IndexType i = 0; i < 4; i++) {
      QuadLaO.global_edge_unsorted(lineglob, quad(ic), i);
      LineHangList::iterator p = LineHang.find(lineglob);

      if (p != LineHang.end()) {
        int cn = p->second.cneighbour();
        int rn = p->second.rneighbour();

        // coarse neighbour

        if ((cn == -1) && (rn != ic)) {
          p->second.cneighbour() = ic;
        }

        // refined neighbour

        if ((rn == -1) && (cn != ic)) {
          p->second.rneighbour() = ic;
        }
      }
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::basic_fill_neighbours2d()
{
  EdgeVector lineglob;
  for (IndexType ic = 0; ic < quads.size(); ic++) {
    for (IndexType i = 0; i < 4; i++) {
      QuadLaO.global_edge_unsorted(lineglob, quad(ic), i);
      LineHangList::iterator p = LineHang.find(lineglob);

      if (p != LineHang.end()) {
        int cn = p->second.cneighbour();
        int rn = p->second.rneighbour();

        if ((cn == -1) && (!quads[ic].sleep())) {
          p->second.cneighbour() = ic;
        } else if (rn == -1) {
          p->second.rneighbour() = ic;
        }
      }
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::basic_refine2d(HangContainer2d& hangset,
                                   const IndexSet& CellRefList,
                                   const IndexSet& CellCoarseList)
{
  IndexType ov = nnodes();
  IndexType oc = quads.size();
  IndexType oe = edges.size();

  IndexType csub = 4 * CellCoarseList.size();
  IndexType cadd = 4 * CellRefList.size();

  IndexType vsub = hangset.NToBeDeleted() + CellCoarseList.size();
  IndexType vadd = hangset.NToBeCreated() + CellRefList.size();

  IndexType cdiff = cadd - csub;
  IndexType vdiff = vadd - vsub;

  IndexType nv = ov + vdiff;
  IndexType nc = oc + cdiff;

  //   cerr << "nv ov vdiff " << nv <<" "<< ov<<" "<< vdiff<<endl;

  clear_transfer_lists();

  IndexVector cdel, vdel, edel;

  for (IntSetCIt p = CellCoarseList.begin(); p != CellCoarseList.end(); p++) {
    const Quad& q = quad(*p);
    vdel.push_back(QuadLaO.middle_vertex(q));
    cdel.insert(cdel.end(), q.childs().begin(), q.childs().end());
  }
  hangset.load_elimination(vdel);
  hangset.NeighbourSwapper();

  //////
  EdgeManager EM(edges, quads, co2n, eo2n);
  //////

  EM.LoadEdgeElimination(edel, CellCoarseList, hangset);

  IndexSet LineRefList, LineCoarseList, ccdel;

  boundary_prepare2d(LineRefList, LineCoarseList, ccdel, hangset);

  transfer(oc, co2n, cdel);
  transfer(ov, vo2n, vdel);
  transfer(oe, eo2n, edel);

  IndexVector cnew(cadd), vnew(vadd);
  iota(cnew.begin(), cnew.end(), oc - csub);
  iota(vnew.begin(), vnew.end(), ov - vsub);

  delete_vertexs2d(vo2n);
  /////
  EM.DeleteEdges();
  /////
  delete_cells<Quad>(CellCoarseList, quads, co2n, vo2n);

  vertexs2d.reserve(nv);
  vertexs2d.resize(nv);
  quads.reserve(nc);
  quads.resize(nc);

  hangset.update_olds(vo2n, co2n);
  hangset.update_news(vnew, CellRefList.size());

  new_vertexs2d(hangset, vnew, CellRefList);
  new_quads(hangset, cnew, vnew, ov, CellRefList);

  basic_fill_neighbours2d();

  new_boundary2d(LineRefList, LineCoarseList, ccdel);

  boundary_newton2d();
  InitQuadOfCurved();
  inner_vertex_newton2d(vnew, CellRefList);

  EM.Build(CellRefList, hangset);
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::boundary_prepare2d(IndexSet& LineRefList,
                                       IndexSet& LineCoarseList,
                                       IndexSet& ccdel,
                                       const HangContainer2d& hangset)
{
  EdgeVector lineglob;
  for (IndexType i = 0; i < Blines.size(); i++) {
    const BoundaryLine& bl = Blines[i];

    lineglob[0] = bl.vertex(0);
    lineglob[1] = bl.vertex(1);
    sort(lineglob.begin(), lineglob.end());

    if (bl.sleep()) {
      if (hangset.ToBeDeleted(lineglob)) {
        LineCoarseList.insert(i);
        for (IndexType j = 0; j < 2; j++)
          ccdel.insert(bl.child(j));
      }
    } else {
      if (hangset.ToBeCreated(lineglob)) {
        LineRefList.insert(i);
      }
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::new_boundary2d(IndexSet& ref,
                                   IndexSet& coarse,
                                   IndexSet& ccdel)
{
  IndexType oc = Blines.size();
  IndexType csub = 2 * coarse.size();
  IndexType cadd = 2 * ref.size();
  IndexType nc = oc + cadd - csub;

  IndexVector lo2n;
  transfer(oc, lo2n, ccdel);
  delete_cells<BoundaryLine>(coarse, Blines, lo2n, vo2n);

  update_boundary_data2d(coarse);

  IndexVector cnew(cadd);
  iota(cnew.begin(), cnew.end(), oc - csub);

  Blines.reserve(nc);
  Blines.resize(nc);

  new_lines(lo2n, cnew, ref);
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::new_quads(const HangContainer2d& hangset,
                              const IndexVector& cnew,
                              const IndexVector& vnew,
                              IndexType nvold,
                              const IndexSet& CellRefList)
{
  // cerr << "new_quads()" << endl;
  // neue zellen erzeugen
  // eintragen der "Vater-Vertexs" in den kindern

  IndexType nci = 0;
  IndexType ivm = 0;

  IntSetIt cp = CellRefList.begin();
  EdgeVector lineglob;
  while (cp != CellRefList.end()) {
    IndexType father = co2n[*cp++];

    assert(father >= 0);

    vector<IndexType>& qc = quads[father].childs();
    qc.resize(4);
    IndexType material = quads[father].material();

    IndexType childlevel = quads[father].level() + 1;
    for (IndexType ic = 0; ic < 4; ic++) {
      IndexType inold = cnew[nci + ic];
      qc[ic] = inold;
      quads[inold].level() = childlevel;
      quads[inold].father() = father;
      quads[inold].material() = material;
      quads[inold].childs().resize(0);
      quads[inold].edges().fill(-1);
    }
    QuadLaO.fill_corner_vertex_in_childs(quads[father]);
    QuadLaO.fill_middle_vertex_in_childs(quads[father], vnew[ivm]);
    ivm++;
    nci += 4;

    // Edge Vertex -- linehanginfo schon ok (hanging) !
    IndexType ive(-1);
    for (IndexType i = 0; i < 4; i++) {
      QuadLaO.global_edge_unsorted(lineglob, quad(father), i);

      ive = hangset.vertex_index(lineglob);

      assert(ive >= 0);

      QuadLaO.fill_edge_vertex_in_childs(quads[father], i, ive);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::inner_vertex_newton2d(const IndexVector& vnew,
                                          const IndexSet& CellRefList)
{
  if (GetCurvedShapes().empty())
    return;

  std::array<IndexType, 2> w;
  IntSetIt cp = CellRefList.begin();

  for (IndexType i = 0; i < CellRefList.size(); i++) {
    IndexType oldind = *cp;
    IndexType ind = co2n[*cp];
    cp++;
    if (quadofcurved.find(oldind) == quadofcurved.end())
      continue;

    const Quad& f = quad(ind);

    // alter version
    //       for (IndexType j=0; j<4; j++) v[j] =
    //       QuadLawOrder().edge_vertex(f,j); new_face_vertex2d(vnew[i],v);

    // neue Version
    //
    IndexType bl = quadofcurved.find(oldind)->second;
    IndexType ei = bline(bl).edge_in_quad();
    IndexType ei2 = (ei + 2) % 4;
    w[0] = QuadLawOrder().edge_vertex(f, ei);
    w[1] = QuadLawOrder().edge_vertex(f, ei2);
    new_edge_vertex2d(vnew[i], w);
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::new_vertexs2d(HangContainer2d& hangset,
                                  const IndexVector& vnew,
                                  const IndexSet& CellRefList)
{
  IndexType nv1 = CellRefList.size();

  IntSetIt cp = CellRefList.begin();
  for (IndexType i = 0; i < nv1; i++) {
    IndexType f = co2n[*cp++];

    new_face_vertex2d(vnew[i], quads[f]);
  }
  HangList<2>::const_iterator p = hangset.Creating().begin();

  for (; p != hangset.Creating().end(); p++) {
    new_edge_vertex2d(p->second.hanging(), p->first);
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::init_line(BoundaryLine& newline)
{
  std::array<IndexType, 2> v;
  for (IndexType i = 0; i < quads.size(); i++) {
    for (IndexType edge = 0; edge < 4; edge++) {
      v[0] = quad(i).vertex(edge);
      v[1] = quad(i).vertex((edge + 1) % 4);
      if (newline == v) {
        newline.of_quad() = i;
        newline.edge_in_quad() = edge;
        return;
      }
      /*	  w[0] = v[1];
                w[1] = v[0];
                if (newline==w)
                  {
                    newline.of_quad()      = i;
                    newline.edge_in_quad() = edge;
                    newline[0] = v[0];
                    newline[1] = v[1];
                    return;
                  }
      */
    }
  }
  cerr << "Sophie im Brunnen !" << endl;
  cerr << newline;
  abort();
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::update_boundary_data2d(const IndexSet& LCoarse)
{
  IndexType no = Blines.size() - 2 * LCoarse.size();
  for (IndexType i = 0; i < no; ++i) {
    IndexType oq = Blines[i].of_quad();
    assert(co2n[oq] >= 0);

    Blines[i].of_quad() = co2n[oq];
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::boundary_newton2d()
{
  if (GetCurvedShapes().empty())
    return;

  for (IndexType i = 0; i < Blines.size(); i++) {
    BoundaryLine& bl = Blines[i];
    IndexType color = bl.material();

    if (GetCurvedShapes().Curved(color)) {
      GetCurvedShapes().newton(color, vertexs2d[bl.vertex(0)]);
      GetCurvedShapes().newton(color, vertexs2d[bl.vertex(1)]);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::new_lines(const IndexVector& lo2n,
                              const IndexVector& cnew,
                              const IndexSet& LRef)
{
  IndexType nci = 0;

  for (IntSetCIt cp = LRef.begin(); cp != LRef.end(); ++cp) {
    IndexType father = lo2n[*cp];

    assert(father >= 0);

    BoundaryLine& blf = Blines[father];
    // change father boundary
    vector<IndexType>& qc = blf.childs();
    qc.resize(2);

    IndexType edge = blf.edge_in_quad();
    IndexType iq = blf.of_quad();
    IndexType vm = QuadLaO.edge_vertex(quad(iq), edge);
    std::array<IndexType, 2> chvec;
    QuadLaO.childs_of_edge(chvec, quad(iq), edge);

    for (IndexType ic = 0; ic < 2; ic++) {
      IndexType inold = cnew[nci + ic];
      // set childs in father
      qc[ic] = inold;
      // set properties of childs
      BoundaryLine& bl = Blines[inold];
      bl.level() = blf.level() + 1;
      bl.father() = father;
      bl.material() = blf.material();
      bl.childs().resize(0);

      bl.vertex(ic) = blf.vertex(ic);
      bl.vertex((ic + 1) % 2) = vm;

      bl.of_quad() = chvec[ic];
      bl.edge_in_quad() = QuadLaO.ChildEdge(edge);
    }

    nci += 2;
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::refine(const IndexVector& cell_ref,
                           const IndexVector& cell_coarse)
{
  IndexSet CellRefList, CellCoarseList;
  _refine2d(CellRefList, CellCoarseList, cell_ref, cell_coarse);
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::_refine2d(IndexSet& CellRefList,
                              IndexSet& CellCoarseList,
                              const IndexVector& cell_ref,
                              const IndexVector& cell_coarse)
{
  prepare2d(cell_ref, cell_coarse, CellRefList, CellCoarseList);

  while (regular_grid2d_three_refine(CellRefList)) {
  }
  while (regular_grid2d_three_coarse(CellRefList, CellCoarseList)) {
  }

  HangContainer2d hangset(LineHang);

  ghost2d(hangset, CellRefList, CellCoarseList);

  basic_refine2d(hangset, CellRefList, CellCoarseList);

  post_refine2d();
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::delete_vertexs2d(const IndexVector& vo2n)
{
  for (unsigned oi = 0; oi < vo2n.size(); ++oi) {
    IndexType ni = vo2n[oi];
    if (ni >= 0) {
      vertexs2d[ni] = vertexs2d[oi];
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::post_refine2d()
{
  check_mesh2d();
  mnlevels = 0;
  for (IndexType i = 0; i < quads.size(); i++) {
    mnlevels = std::max(mnlevels, quads[i].level());
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::new_edge_vertex2d(IndexType nv, const EdgeVector& v)
{
  vertexs2d[nv].equ(0.5, vertexs2d[v[0]], 0.5, vertexs2d[v[1]]);
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::new_face_vertex2d(IndexType nv, const FaceVector& v)
{
  vertexs2d[nv].equ(0.25,
                    vertexs2d[v[0]],
                    0.25,
                    vertexs2d[v[1]],
                    0.25,
                    vertexs2d[v[2]],
                    0.25,
                    vertexs2d[v[3]]);
}

/*---------------------------------------------------*/

ostream&
operator<<(ostream& os, const pair<std::array<IndexType, 2>, Hang>& H)
{
  cerr << H.first << " -> " << H.second << endl;
  return os;
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::check_mesh2d() const
{
  // check quads

  IndexType cmin = 0, cmax = quads.size();
  IndexType vmin = 0, vmax = nnodes();

  for (vector<Quad>::const_iterator p = quads.begin(); p != quads.end(); p++) {
    const Quad& q = *p;

    // check vertex id
    for (IndexType i = 0; i < q.nvertexs(); i++) {
      if ((q.vertex(i) < vmin) || (q.vertex(i) > vmax)) {
        cerr << "Vertex invalid in Cell: " << *p << " : ";
        cerr << q.vertex(i) << endl;
        exit(1);
      }
    }
    // check child id
    for (IndexType i = 0; i < q.nchilds(); i++) {
      if ((q.child(i) < cmin) || (q.child(i) > cmax)) {
        cerr << "Chid invalid in Cell: " << *p << " : ";
        cerr << q.child(i) << endl;
        exit(1);
      }
    }
  }

  // check linehang
  LineHangList::const_iterator hp;

  for (hp = LineHang.begin(); hp != LineHang.end(); ++hp) {
    assert(hp->second.rneighbour() >= 0);
  }
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh2d::QuadNeighbour(const Quad& q, IndexType e) const
{
  IndexType ie = q.edge(e);
  IndexType m = edges[ie].master();
  if (q == quads[m])
    return edges[ie].slave();
  return m;
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::init_edges2d()
{
  EdgeManager EM(edges, quads, co2n, eo2n);
  EM.InitEdges();
  EM.SortHangings();
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::GetVertexesOfEdge(std::array<IndexType, 3>& v,
                                      IndexType e) const
{
  const Edge& E = edge(e);
  const Quad* Q = &quad(E.master());

  IndexType le = E.LocalMasterIndex();
  v[0] = (*Q)[le];
  v[1] = (*Q)[(le + 1) % 4];
  v[2] = -1;

  if (!Q->sleep()) {
    if (E.slave() == -1)
      return;
    Q = &quad(E.slave());
    le = E.LocalSlaveIndex();
  }
  v[2] = QuadLaO.edge_vertex(*Q, le);
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::GetVertexesOfEdge(std::array<IndexType, 2>& v,
                                      IndexType e) const
{
  const Edge& E = edge(e);
  const Quad* Q = &quad(E.master());

  IndexType le = E.LocalMasterIndex();
  v[0] = (*Q)[le];
  v[1] = (*Q)[(le + 1) % 4];
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::GetAwakeCells(set<IndexType>& v) const
{
  for (IndexType i = 0; i < ncells(); i++) {
    if (!quads[i].sleep())
      v.insert(i);
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::GetAwakePatchs(set<IndexType>& v) const
{
  v.clear();

  for (IndexType i = 0; i < ncells(); i++) {
    if (!quads[i].sleep()) {
      int f = quads[i].father();
      assert(f != -1);
      v.insert(f);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::ConstructQ2PatchMesh(IndexVector& q2patchmesh) const
{
  typedef IndexSet::iterator It;
  IndexSet patche;
  GetAwakePatchs(patche);
  vector<IndexSet> patch_on_level(nlevels());
  q2patchmesh.resize(0);
  for (It it = patche.begin(); it != patche.end(); ++it)
    patch_on_level[quads[*it].level()].insert(*it);
  // grobgitterzellen kommen ins patchmesh
  for (It it = patch_on_level[0].begin(); it != patch_on_level[0].end(); ++it)
    q2patchmesh.push_back(*it);
  // der Rest wird eins groeber
  for (IndexType l = 1; l < nlevels(); ++l) {
    It it = patch_on_level[l].begin();
    while (it != patch_on_level[l].end()) {
      int v = quads[*it].father();
      assert(v != -1);
      q2patchmesh.push_back(v);
      IndexVector nk = Nachkommen(v);
      for (IndexType i = 0; i < nk.size(); ++i)
        patch_on_level[quads[nk[i]].level()].erase(nk[i]);
      it = patch_on_level[l].begin();
    }
  }
}

/*---------------------------------------------------*/

IndexVector
HierarchicalMesh2d::ConstructQ4Patch(IndexType c) const
{
  IndexVector patch(25, -1);

  for (int i = 0; i < 25; i++) {
    // Vertex i steht an Position (x,y)
    int x = i % 5;
    int y = i / 5;

    // Position von erstem Kind
    int fcx = x / 3;
    int fcy = y / 3;
    // Index davon
    int fci = fcy * 2 + abs(fcx - fcy);

    // Position vom Kind im Kind
    int scx = (x - 2 * fcx) / 2;
    int scy = (y - 2 * fcy) / 2;
    // Index davon
    int sci = scy * 2 + abs(scx - scy);

    // Position des Vertex
    int vx = x - 2 * fcx - scx;
    int vy = y - 2 * fcy - scy;
    // Index davon
    int vi = vy * 2 + abs(vx - vy);

    patch[i] = static_cast<IndexType>(
      quads[quads[quads[c].child(fci)].child(sci)].vertex(vi));
  }
  return patch;
}

/*---------------------------------------------------*/

set<IndexType>
HierarchicalMesh2d::CellNeighbours(IndexType iq) const
{
  // soll nur wache Zellen zurueckgeben !!
  set<IndexType> neighbors;
  const Quad& q = quad(iq);
  if (q.sleep())
    return neighbors;

  for (IndexType i = 0; i < 4; i++) {
    IndexType ge = q.edge(i);
    const Edge& e = edge(ge);
    int m = e.master();
    int s = e.slave();
    if (m != -1) {
      if (!quad(m).sleep())
        neighbors.insert(m);
      else {
        for (IndexType ii = 0; ii < 4; ii++)
          neighbors.insert(quad(m).child(ii));
      }
    }
    // wenn slave groeber, koennen keine Nachbarn gefunden werden !!!
    if (s != -1) {
      if (!quad(s).sleep())
        neighbors.insert(s);
      else {
        for (IndexType ii = 0; ii < 4; ii++)
          neighbors.insert(quad(s).child(ii));
      }
    }
  }
  return neighbors;
}

/*------------------------------------------------------*/

IndexType
HierarchicalMesh2d::neighbour(IndexType c, IndexType le) const
{
  assert(le < 4);
  const Quad& Q = quad(c);
  assert(Q.edge(le) >= 0);
  const Edge& E = edge(Q.edge(le));
  IndexType m = E.master();
  IndexType nq = m;
  if (m == c)
    nq = E.slave();
  return nq;
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::FillAllBoundaryLines()
{
  HangBLList B;
  // erstmal alle rein
  EdgeVector edge;
  for (IndexType i = 0; i < nblines(); i++) {
    BoundaryLine b = bline(i);
    edge = b.vertex();
    sort(edge.begin(), edge.end());
    B.insert(make_pair(edge, b));
  }
  for (IndexType i = 0; i < ncells(); i++) {
    Quad q = quad(i);
    for (IndexType e = 0; e < 4; ++e) {
      QuadLawOrder().global_edge_unsorted(edge, q, e);
      sort(edge.begin(), edge.end());
      HangBLList::iterator p = B.find(edge);
      if (p == B.end()) {
        BoundaryLine nb;
        nb.of_quad() = i;
        nb.edge_in_quad() = e;
        nb.material() = -1;
        B.insert(make_pair(edge, nb));
      } else {
        BoundaryLine& b = p->second;
        if (b.material() == -1)
          B.erase(p);
        else {
          b.of_quad() = i;
          b.edge_in_quad() = e;
        }
      }
    }
  }
  IndexType n = B.size();
  Blines.resize(n);
  HangBLList::iterator p = B.begin();
  for (IndexType i = 0; i < n; i++) {
    Blines[i] = p->second;
    if (Blines[i].material() == -1) {
      Blines[i].vertex() = p->first;
    }
    p++;
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::WriteAll(const string& name) const
{
  ofstream out(name.c_str());

  out << dimension() << " dimension" << endl;
  out << nnodes() << " vertexs" << endl;
  out << mnlevels << " mnlevels" << endl;

  out << "vertexs2d\n" << vertexs2d << endl;

  out << "vo2n\n" << vo2n << endl;
  out << "co2n\n" << co2n << endl;
  out << "eo2n\n" << eo2n << endl;

  out << "quads\n" << quads << endl;
  out << "Blines\n" << Blines << endl;
  out << "LineHang\n" << LineHang << endl;

  out.close();
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::write_inp(const string& bname) const
{
  string name = bname;
  IndexType name_size = name.size();
  if (name_size < 4)
    name += ".inp";
  if (name.substr(name_size - 4, 4) != ".inp") {
    name += ".inp";
  }

  ofstream file(name.c_str());
  if (!file.is_open()) {
    cerr << "HierarchicalMesh2d::write_inp()\n";
    cerr << "cannot open file " << name << endl;
    abort();
  }

  IndexType nt = ncells() + nblines();
  file << nnodes() << " " << nt << " " << 0 << " " << 0 << " " << 0 << endl;

  for (IndexType i = 0; i < nnodes(); i++)
    file << i << " " << vertex2d(i) << " " << 0. << endl;

  for (IndexType i = 0; i < ncells(); i++) {
    file << i << " " << 0 << " quad " << quad(i).vertex() << endl;
  }
  for (IndexType i = 0; i < nblines(); i++) {
    file << i << " " << bline(i).material() << " line " << bline(i).vertex()
         << endl;
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::write_vtk(const string& name) const
{
  ofstream file(name.c_str());
  if (!file.is_open()) {
    cerr << "HierarchicalMesh2d::write_vtk()\n";
    cerr << "cannot open file " << name << endl;
    exit(1);
  }

  IndexType nn = nnodes();

  file << "# vtk DataFile Version 2.4 " << endl;
  file << "output from GascoigneStd" << endl;
  file << "ASCII" << endl;
  file << "DATASET UNSTRUCTURED_GRID" << endl;
  file << "POINTS " << nn << " FLOAT " << endl;

  IndexType ne = ncells();
  for (IndexType i = 0; i < nn; i++) {
    file << vertex2d(i) << " " << 0 << endl;
  }

  file << endl << "CELLS " << ne << " " << 5 * ne << endl;

  for (IndexType c = 0; c < ne; c++) {
    file << 4 << " ";
    for (IndexType i = 0; i < 4; i++) {
      file << vertex_of_cell(c, i) << " ";
    }
    file << endl;
  }
  file << endl << "CELL_TYPES " << ne << endl;
  for (IndexType i = 0; i < ne; i++)
    file << 9 << " ";
}

/*---------------------------------------------------*/

pair<bool, triple<IndexType, IndexType, IndexType>>
HierarchicalMesh2d::check_inp(const string& name)
{
  //  detect some errors in input-file... and more

  ifstream file(name.c_str());
  if (!file.is_open()) {
    cerr << "HierarchicalMesh2d::check_inp()\n";
    cerr << "cannot open file " << name << endl;
    exit(1);
  }

  bool first_one = 1;
  IndexType nv, nl, nq, nt;
  IndexType n_unkonwn;
  file >> nv >> nt >> n_unkonwn >> n_unkonwn >> n_unkonwn;

  Vertex3d c;
  IndexType ind;
  for (IndexType i = 0; i < nv; i++) {
    file >> ind >> c;
  }

  nq = 0;
  nl = 0;

  std::array<IndexType, 4> iq;
  std::array<IndexType, 2> il;
  for (IndexType i = 0; i < nt; i++) {
    string name;
    string mat;
    IndexType ii;
    file >> ii >> mat >> name;
    if (name == "quad") {
      file >> iq;
      nq++;
      if ((iq[0] == 0) || (iq[1] == 0) || (iq[2] == 0) || (iq[3] == 0)) {
        first_one = 0;
      }
    } else if (name == "line") {
      file >> il;
      nl++;
      if ((il[0] == 0) || (il[1] == 0)) {
        first_one = 0;
      }
    }
  }

  // fehlerabfragen ....
  if (nt != (nl + nq)) {
    cerr << "wrong number of cells: " << nt << endl;
    exit(1);
  }

  file.close();

  return make_pair(first_one,
                   make_triple<IndexType, IndexType, IndexType>(nl, nq, 0));
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::read_inp(const string& name)
{
  // check mesh.....

  pair<bool, tint> p = check_inp(name);
  bool first_one = p.first;
  tint n = p.second;

  IndexType nl = n.first;
  IndexType nq = n.second;

  ifstream file(name.c_str());
  if (!file.is_open()) {
    cerr << "HierarchicalMesh2d::read_inp()\n";
    cerr << "cannot open file " << name << endl;
    abort();
  }

  IndexType nv, nt, n_unkonwn;

  file >> nv >> nt >> n_unkonwn >> n_unkonwn >> n_unkonwn;

  if (_i_showoutput) {
    cout << "2D Mesh: " << nv << " nodes, ";
    cout << nl << " lines, ";
    cout << nq << " quads" << endl;
  }

  vertexs2d.reserve(nv);
  vertexs2d.resize(nv, Vertex2d());
  quads.reserve(nq);
  quads.resize(nq, Quad());
  Blines.reserve(nl);
  Blines.resize(nl);

  Vertex2d c;
  IndexType ind;
  double z;
  for (IndexType i = 0; i < nv; i++) {
    file >> ind >> c >> z;
    vertexs2d[i] = c;
  }

  std::array<IndexType, 4> iqv;
  std::array<IndexType, 2> ilv;
  IndexType iq = 0;
  IndexType il = 0;
  for (IndexType i = 0; i < nt; i++) {
    string name;
    IndexType unknown;
    string matstring;
    file >> unknown >> matstring >> name;
    if (name == "quad") {
      file >> iqv;
      if (first_one)
        for (IndexType iii = 0; iii < 4; iii++)
          iqv[iii]--;
      quads[iq].vertex() = iqv;
      quads[iq].material() = atoi(matstring.c_str());
      iq++;
    } else if (name == "line") {
      file >> ilv;
      if (first_one) {
        ilv[0]--;
        ilv[1]--;
      }

      BoundaryLine li;
      li.material() = atoi(matstring.c_str());
      li.vertex() = ilv;
      init_line(li);
      Blines[il++] = li;
    }
  }
  init_edges2d();
  IndexVector nix(0);
  refine(nix, nix);
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::write(const string& bname) const
{
  write_gup(bname);
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::write_gup(const string& bname) const
{
  string name = bname;
  IndexType name_size = name.size();
  if (name_size < 4)
    name += ".gup";
  if (name.substr(name_size - 4, 4) != ".gup") {
    name += ".gup";
  }

  ofstream out(name.c_str());

  out.precision(16);
  out << dimension() << " dimension" << endl;
  out << nnodes() << " vertexs" << endl;

  for (IndexType i = 0; i < nnodes(); i++) {
    out << " " << vertex2d(i) << endl;
  }
  out << quads.size() << " quads" << endl;

  for (IndexType i = 0; i < quads.size(); i++) {
    out << quad(i);
  }
  out << LineHang << endl;

  out << Blines.size() << " boundarylines" << endl;
  for (IndexType i = 0; i < Blines.size(); i++) {
    out << Blines[i].material() << " " << Blines[i];
  }
  out << endl << endl << edges.size() << " edges" << endl;
  for (IndexType i = 0; i < edges.size(); i++) {
    out << " " << edges[i];
  }
  out << endl;
  out.close();
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::write_gip(const string& bname) const
{
  string name = bname;
  IndexType name_size = name.size();
  if (name_size < 4)
    name += ".gip";
  if (name.substr(name_size - 4, 4) != ".gip") {
    name += ".gip";
  }

  ofstream out(name.c_str(), ios_base::out | ios_base::binary);

  out.precision(16);
  IndexType dim = dimension(), n = nnodes(), sizeInt = sizeof(IndexType);
  out.write(reinterpret_cast<const char*>(&dim), sizeInt);
  out.write(reinterpret_cast<const char*>(&n), sizeInt);

  for (IndexType i = 0; i < nnodes(); i++) {
    ArrayBinWrite(out, vertex2d(i));
  }

  IndexType nquads = quads.size();
  out.write(reinterpret_cast<const char*>(&nquads), sizeInt);

  for (IndexType i = 0; i < quads.size(); i++) {
    ArrayBinWrite(out, quad(i));
  }

  LineHang.BinWrite(out);
  IndexType nblines = Blines.size();
  out.write(reinterpret_cast<const char*>(&nblines), sizeInt);
  for (IndexType i = 0; i < Blines.size(); i++) {
    IndexType mat = Blines[i].material();
    out.write(reinterpret_cast<const char*>(&mat), sizeInt);
    ArrayBinWrite(out, Blines[i]);
  }
  IndexType nedges = edges.size();
  out.write(reinterpret_cast<const char*>(&nedges), sizeInt);
  for (IndexType i = 0; i < edges.size(); i++) {
    edges[i].BinWrite(out);
  }
  out.close();
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::read_gup(const string& bname)
{
  string name = bname;
  IndexType name_size = name.size();
  if (name_size < 4)
    name += ".gup";
  if (name.substr(name_size - 4, 4) != ".gup") {
    name += ".gup";
  }

  string symbol;

  ifstream file(name.c_str());

  if (!file.is_open()) {
    cerr << "HierarchicalMesh2d::read_gup(): error in file " << name << endl;
    abort();
  }

  Blines.clear();
  edges.clear();
  LineHang.clear();

  IndexType n, dim;
  file >> dim >> symbol >> n >> symbol;

  assert(dim == 2);

  if (_i_showoutput) {
    cout << "Mesh 2d  : ";
  }
  vertexs2d.reserve(n);
  vertexs2d.resize(n);

  if (symbol != "vertexs") {
    cout << "HierarchicalMesh2d::read error" << endl;
    exit(1);
  }
  for (IndexType i = 0; i < n; i++) {
    file >> vertexs2d[i];
  }
  if (_i_showoutput) {
    cout << n << " nodes, ";
  }

  file >> n >> symbol;

  if (symbol != "quads") {
    cout << "HierarchicalMesh2d::read error 2" << endl;
    exit(1);
  }
  quads.reserve(n);
  quads.resize(n);
  for (IndexType i = 0; i < quads.size(); i++) {
    file >> quads[i];
  }
  if (_i_showoutput) {
    cout << n << " quads, ";
  }
  file >> LineHang;
  if (_i_showoutput) {
    cout << LineHang.size() << " hangs, ";
  }

  file >> n >> symbol;
  IndexType number = 0;
  if (symbol == "boundary_materials") {
    /* old format */
    BoundaryLine bol;
    for (IndexType i = 0; i < n; i++) {
      file >> symbol;
      if (symbol != "material") {
        cout << "HierarchicalMesh2d::read error 4" << endl;
        exit(1);
      }
      IndexType nn;
      file >> bol.material() >> nn;
      for (IndexType j = 0; j < nn; j++) {
        file >> bol;
        Blines.push_back(bol);
      }
      number += nn;
    }
  } else if (symbol == "boundarylines") {
    /* new format */
    BoundaryLine bol;
    for (IndexType i = 0; i < n; i++) {
      file >> bol.material() >> bol;
      Blines.push_back(bol);
    }
    number = n;
  } else {
    cout << "HierarchicalMesh2d::read error 3" << endl;
    exit(1);
  }
  if (_i_showoutput) {
    cout << number << " lines, ";
  }
  file >> n >> symbol;
  if (_i_showoutput) {
    cout << n << " edges" << endl;
  }
  if (symbol == "edges") {
    Edge e;
    for (IndexType i = 0; i < n; i++) {
      file >> e;
      edges.push_back(e);
    }
  }
  if (edges.size() == 0)
    init_edges2d();
  post_refine2d();
  file.close();
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::read_gip(const string& bname)
{
  string name = bname;
  IndexType name_size = name.size();
  if (name_size < 4)
    name += ".gip";
  if (name.substr(name_size - 4, 4) != ".gip") {
    name += ".gip";
  }

  string symbol;

  ifstream file(name.c_str(), ios_base::in | ios_base::binary);

  if (!file.is_open()) {
    cerr << "HierarchicalMesh2d::read_gip(): error in file " << name << endl;
    abort();
  }

  Blines.clear();
  edges.clear();
  LineHang.clear();

  IndexType n, dim, sizeInt = sizeof(IndexType);
  file.read(reinterpret_cast<char*>(&dim), sizeInt);
  file.read(reinterpret_cast<char*>(&n), sizeInt);

  assert(dim == 2);

  if (_i_showoutput) {
    cout << "Mesh 2d  : ";
  }
  vertexs2d.reserve(n);
  vertexs2d.resize(n);

  for (IndexType i = 0; i < n; i++) {
    ArrayBinRead(file, vertexs2d[i]);
  }
  if (_i_showoutput) {
    cout << n << " nodes, ";
  }
  file.read(reinterpret_cast<char*>(&n), sizeInt);
  quads.reserve(n);
  quads.resize(n);
  for (IndexType i = 0; i < quads.size(); i++) {
    ArrayBinRead(file, quads[i]);
  }
  if (_i_showoutput) {
    cout << n << " quads, ";
  }
  LineHang.BinRead(file);
  if (_i_showoutput) {
    cout << LineHang.size() << " hangs, ";
  }
  file.read(reinterpret_cast<char*>(&n), sizeInt);
  IndexType number = 0;
  BoundaryLine bol;
  for (IndexType i = 0; i < n; i++) {
    file.read(reinterpret_cast<char*>(&bol.material()), sizeInt);
    ArrayBinRead(file, bol);
    Blines.push_back(bol);
  }
  number = n;
  if (_i_showoutput) {
    cout << number << " lines, ";
  }
  file.read(reinterpret_cast<char*>(&n), sizeInt);
  if (_i_showoutput) {
    cout << n << " edges" << endl;
  }
  Edge e;
  for (IndexType i = 0; i < n; i++) {
    e.BinRead(file);
    edges.push_back(e);
  }
  if (edges.size() == 0)
    init_edges2d();
  post_refine2d();
  file.close();
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh2d::regular_grid2d_one(IndexSet& celllist,
                                       IndexVector& coarsesub,
                                       IndexSet& CellRefList,
                                       IndexSet& CellCoarseList) const
{
  /* detects jump over two levels across LineHangs */

  IndexType n = 0;
  LineHangList::const_iterator hp = LineHang.begin();

  for (; hp != LineHang.end(); ++hp) {
    int cr = hp->second.rneighbour();
    int cn = hp->second.cneighbour();

    assert(cr >= 0);

    if (cn != -1) {
      if (quad(cr).childs().size() == 0)
        continue;

      std::array<IndexType, 2> f;
      QuadLaO.childs_of_global_edge(f, quad(cr), hp->first);

      for (unsigned i = 0; i < f.size(); ++i) {
        IndexType c = f[i];
        if ((CellRefList.find(c) != CellRefList.end()) || (quads[c].sleep())) {
          if (CellCoarseList.find(cn) != CellCoarseList.end()) {
            coarsesub.push_back(cn);
            break;
          } else {
            assert(!quad(cn).sleep());

            pair<IntSetIt, bool> p = celllist.insert(cn);
            if (p.second) {
              n++;
              break;
            }
          }
        }
      }
    }
  }
  // unique(coarsesub.begin(),coarsesub.end());
  return n + coarsesub.size();
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh2d::regular_grid2d_two(IndexSet& celllist,
                                       IndexSet& CellRefList) const
{
  /* detects more than 2 HangLines on one Quad */

  IndexType nto = 0;
  IndexVector nh(quads.size());
  LineHangList::const_iterator p = LineHang.begin();
  for (; p != LineHang.end(); p++) {
    IndexType i = p->second.cneighbour();
    if (i < 0)
      continue;
    if (quad(i).sleep())
      continue;
    if (CellRefList.find(i) != CellRefList.end())
      continue;

    nh[i]++;
  }
  for (IndexType i = 0; i < quads.size(); i++) {
    if (nh[i] > 2) {
      pair<IntSetIt, bool> pp = celllist.insert(i);
      if (pp.second)
        nto++;
    }
  }
  return nto;
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh2d::smooth_edges()
{
  IndexSet CellRefList, CellCoarseList, coarseadd;
  HangContainer2d hangset(LineHang);

  IndexType r = 0; // regular_grid2d_three(CellRefList);

  ghost2d(hangset, CellRefList, coarseadd);

  basic_refine2d(hangset, CellRefList, CellCoarseList);
  post_refine2d();

  return r;
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh2d::regular_grid2d_three(IndexSet& CellRef,
                                         IndexSet& CellCoarse) const
{
  vector<IndexType> maxlevel(nnodes());
  vector<IndexType> minlevel(nnodes(), 1000);

  // set maximal levels for vertices
  for (IndexType i = 0; i < quads.size(); i++) {
    const Quad& q = quad(i);
    if (q.sleep())
      continue;

    IndexType lev = q.level();
    if (CellRef.find(i) != CellRef.end())
      lev++;

    IndexType father = q.father();
    if (father >= 0) {
      if (CellCoarse.find(father) != CellCoarse.end()) {
        lev--;
      }
    }
    for (IndexType j = 0; j < 4; j++) {
      IndexType k = q[j];
      maxlevel[k] = std::max(maxlevel[k], lev);
      minlevel[k] = std::min(minlevel[k], lev);
    }
  }
  set<IndexType> cand;
  for (IndexType i = 0; i < nnodes(); i++) {
    if (maxlevel[i] >= minlevel[i] + 2) {
      cand.insert(i);
    }
  }

  IndexType ref = 0;
  for (IndexType i = 0; i < quads.size(); i++) {
    const Quad& q = quad(i);
    if (q.sleep())
      continue;
    IndexType lev = q.level();
    for (IndexType j = 0; j < 4; j++) {
      IndexType k = q[j];
      set<IndexType>::const_iterator p = cand.find(k);
      if (p != cand.end()) {
        if (CellRef.find(i) == CellRef.end()) {
          if (maxlevel[k] >= lev + 2)
            CellRef.insert(i);
          ref++;
        }
      }
    }
  }
  IndexSet coarsesub;
  IndexSet::const_iterator p = CellCoarse.begin();
  for (; p != CellCoarse.end(); p++) {
    const Quad& q = quad(*p);
    IndexType lev = q.level();
    for (IndexType j = 0; j < 4; j++) {
      IndexType k = q[j];
      if (lev <= minlevel[k] - 1) {
        coarsesub.insert(*p);
      }
    }
  }
  IndexType coarse = 0;
  p = coarsesub.begin();
  for (; p != coarsesub.end(); p++) {
    if (CellCoarse.find(*p) != CellCoarse.end()) {
      CellCoarse.erase(*p);
      coarse++;
    }
  }
  // cout << "(" << ref << " " << coarse << ")\n";
  return ref + coarse;
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::GetMinMaxLevels(IndexVector& maxi,
                                    IndexVector& mini,
                                    const IndexSet& CellRef) const
{
  // set maximal levels for vertices
  //
  maxi.resize(nnodes());
  mini.resize(nnodes());
  mini = 1000;
  for (IndexType i = 0; i < quads.size(); i++) {
    const Quad& q = quad(i);
    if (q.sleep())
      continue;

    IndexType lev = q.level();
    if (CellRef.find(i) != CellRef.end())
      lev++;
    for (IndexType j = 0; j < 4; j++) {
      IndexType k = q[j];
      maxi[k] = std::max(maxi[k], lev);
      mini[k] = std::min(mini[k], lev);
    }
  }
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh2d::regular_grid2d_three_refine(IndexSet& CellRef) const
{
  IndexVector maxlevel, minlevel;

  GetMinMaxLevels(maxlevel, minlevel, CellRef);

  set<IndexType> cand;
  for (IndexType i = 0; i < nnodes(); i++) {
    if (maxlevel[i] >= minlevel[i] + 2) {
      cand.insert(i);
    }
  }

  IndexType ref = 0;
  for (IndexType i = 0; i < quads.size(); i++) {
    const Quad& q = quad(i);
    if (q.sleep())
      continue;
    IndexType lev = q.level();
    for (IndexType j = 0; j < 4; j++) {
      IndexType k = q[j];
      set<IndexType>::const_iterator p = cand.find(k);
      if (p != cand.end()) {
        if (CellRef.find(i) == CellRef.end()) {
          if (maxlevel[k] >= lev + 2)
            CellRef.insert(i);
          ref++;
        }
      }
    }
  }
  // cout << "(" << ref << ")";
  return ref;
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh2d::regular_grid2d_three_coarse(IndexSet& CellRef,
                                                IndexSet& CellCoarse) const
{
  IndexType maxl = 0;
  vector<IndexSet> LevelCellCoarse;
  {
    IndexSet::const_iterator p = CellCoarse.begin();
    for (; p != CellCoarse.end(); p++) {
      const Quad& q = quad(*p);
      IndexType lev = q.level();
      maxl = std::max(maxl, lev);
    }
    LevelCellCoarse.resize(maxl + 1);
    p = CellCoarse.begin();
    for (; p != CellCoarse.end(); p++) {
      const Quad& q = quad(*p);
      IndexType lev = q.level();
      LevelCellCoarse[lev].insert(*p);
    }
  }
  IndexType coarse = 0;
  for (int i = maxl; i >= 0; i--) {
    IndexVector maxlevel, minlevel;

    GetMinMaxLevels(maxlevel, minlevel, CellRef);

    IndexSet coarsesub;
    // IndexSet::const_iterator p = LevelCellCoarse[i].begin();
    // for (; p != LevelCellCoarse[i].end(); p++) {
    for (const IndexType& p : LevelCellCoarse[i]) {
      const Quad& q = quad(p);
      IndexType lev = q.level();
      for (IndexType j = 0; j < 4; j++) {
        IndexType k = q[j];
        if (lev + 1 < maxlevel[k]) {
          coarsesub.insert(p);
          continue;
        }
      }
    }
    // // CellCoarse.erase(coarsesub.begin(),coarsesub.end());
    // p = coarsesub.begin();
    // for (; p != coarsesub.end(); p++) {
    for (const IndexType& p : coarsesub) {
      if (CellCoarse.find(p) != CellCoarse.end()) {
        CellCoarse.erase(p);
        coarse++;
      }
    }
  }
  // cout << "(" << coarse << ")";
  return coarse;
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::global_coarse()
{
  IndexSet RefList, CoarseList;

  for (IndexType i = 0; i < quads.size(); i++) {
    const Quad& q = quads[i];
    if (q.sleep())
      continue;
    if (!q.level())
      continue;
    CoarseList.insert(q.father());
  }

  HangContainer2d hangset(LineHang);

  ghostglobalcoarse(hangset, CoarseList);

  basic_refine2d(hangset, RefList, CoarseList);

  post_refine2d();
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::LoadFathers(IndexVector& v) const
{
  IndexSet fathers;
  for (IndexType i = 0; i < v.size(); i++) {
    if (quad(v[i]).level() == 0)
      continue;
    IndexType f = quad(v[i]).father();
    if (f < 0) {
      cerr << "HierarchicalMesh2d::LoadFathers no father\n";
      abort();
    }
    fathers.insert(f);
  }
  Set2Vec(v, fathers);
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::FillVertexLevels(IndexVector& dst) const
{
  dst.resize(nnodes(), 100);
  for (IndexType i = 0; i < ncells(); i++) {
    const Quad& Q = quads[i];
    if (Q.sleep())
      continue;

    IndexType level = Q.level();
    for (IndexType j = 0; j < 4; j++) {
      IndexType k = Q[j];
      dst[k] = std::min(dst[k], level);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::RefineCoarseNodes(IndexSet& dst,
                                      const IndexVector& refnodes,
                                      const IndexVector& vertexlevel) const
{
  IndexSet h;

  Vec2Set(h, refnodes);

  dst.clear();
  for (IndexType i = 0; i < ncells(); i++) {
    const Quad& Q = quads[i];
    if (Q.sleep())
      continue;

    IndexType f = Q.father();
    if (f < 0)
      continue;

    const Quad& QF = quads[f];

    for (IndexType j = 0; j < 4; j++) {
      IndexType k = Q[j];
      if (h.find(k) == h.end())
        continue;

      IndexType minlevel = vertexlevel[QF[0]];
      for (IndexType v = 1; v < 4; v++) {
        minlevel = std::min(minlevel, vertexlevel[QF[v]]);
      }
      for (IndexType v = 0; v < 4; v++) {
        IndexType w = QF[v];
        assert(w < vertexlevel.size());
        if (vertexlevel[w] == minlevel) {
          dst.insert(w);
        }
      }
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::VertexToCells(IndexVector& dst,
                                  const IndexSet& src,
                                  const IndexVector& vertexlevel) const
{
  for (IndexType i = 0; i < ncells(); i++) {
    const Quad& Q = quads[i];
    if (Q.sleep())
      continue;
    IndexType level = Q.level();
    for (IndexType j = 0; j < 4; j++) {
      IndexType k = Q[j];
      if (vertexlevel[k] == level) {
        if (src.find(k) != src.end()) {
          dst.push_back(i);
          break;
        }
      }
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::VertexToCellsCoarsening(
  IndexVector& dst,
  const IndexSet& src,
  const IndexVector& vertexlevel) const
{
  IndexType limit = 2;

  for (IndexType i = 0; i < ncells(); i++) {
    const Quad& Q = quads[i];
    if (Q.sleep())
      continue;
    IndexType count = 0;
    for (IndexType j = 0; j < 4; j++) {
      IndexType k = Q[j];
      if (src.find(k) != src.end()) {
        count++;
      }
    }
    if (count >= limit)
      dst.push_back(i);
  }
}

/*---------------------------------------------------*/

// void HierarchicalMesh2d::double_patch_refine
// (IndexVector& cell_ref, IndexVector& cell_coarse)
// {
//   LoadFathers(cell_ref);
//   LoadFathers(cell_ref);
//   LoadFathers(cell_coarse);
//   LoadFathers(cell_coarse);

//   DoubleCoarseHierarchicalMesh2d CM(*this);

//   CM.BasicInit();
//   CM.refine(cell_ref,cell_coarse);
//   CM.GetRefinedList(cell_ref);
//   CM.GetCoarsedList(cell_coarse);

//   IndexVector ref(0), coarse(0);

//   for (IndexType i=0; i<cell_ref.size(); i++)
//     {
//       const Quad& q = quad(cell_ref[i]);
//       if (!q.sleep())
// 	{
// 	  cout << "HierarchicalMesh2d::patchrefine am dampfen1 !";
// 	  abort();
// 	}
//       for (IndexType j=0; j<4; j++)
// 	{
// 	  IndexType child = q.child(j);
// 	  for (IndexType jj=0; jj<4; jj++)
// 	    {
// 	      IndexType grandchild = quad(child).child(jj);
// 	      ref.push_back(grandchild);
// 	    }
// 	}
//     }
//   cell_coarse.resize(0);
//   refine(ref,cell_coarse);
// }

/*---------------------------------------------------*/

void
HierarchicalMesh2d::recursive_childs(IndexType q,
                                     IndexVector& ref,
                                     IndexType d) const
{
  if (d > 0) {
    const Quad& Q = quad(q);
    assert(Q.sleep());
    for (IndexType j = 0; j < 4; j++) {
      IndexType child = Q.child(j);
      recursive_childs(child, ref, d - 1);
    }
  } else {
    ref.push_back(q);
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::patch_refine(IndexVector& cell_ref,
                                 IndexVector& cell_coarse)
{
  for (IndexType i = 0; i < pdepth; i++) {
    LoadFathers(cell_ref);
    LoadFathers(cell_coarse);
  }

  CoarseHierarchicalMesh2d CM(*this);

  CM.BasicInit(pdepth);
  CM.refine(cell_ref, cell_coarse);
  CM.GetRefinedList(cell_ref);
  CM.GetCoarsedList(cell_coarse);

  IndexVector ref(0), coarse(0);

  for (IndexType i = 0; i < cell_ref.size(); i++) {
    recursive_childs(cell_ref[i], ref, pdepth);
  }
  for (IndexType i = 0; i < cell_coarse.size(); i++) {
    recursive_childs(cell_coarse[i], coarse, pdepth + 1);
  }
  refine(ref, coarse);
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::FillVolumes(DoubleVector& vol) const
{
  vol.resize(ncells(), 0.);
  for (IndexType i = 0; i < ncells(); i++) {
    const Quad& Q = quad(i);
    Vertex2d V = vertex2d(Q[2]);
    V -= vertex2d(Q[0]);
    vol[i] = V * V;
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh2d::writeq2(const IndexVector& a,
                            const vector<IndexType>& b,
                            IndexType np) const
{
  char s[29];
  for (IndexType p = 0; p < np; ++p) {
    sprintf(s, "n_%d_%d", ncells(), p);
    ofstream aus(s);

    for (IndexType i = 0; i < a.size(); ++i) {
      const Quad& Q = quads[a[i]];
      if (b[i] == p) {
        aus << vertex2d(Q[0]) << endl;
        aus << vertex2d(Q[1]) << endl;
        aus << vertex2d(Q[2]) << endl;
        aus << vertex2d(Q[3]) << endl;
        aus << vertex2d(Q[0]) << endl << endl;
      }
    }
    aus.close();
  }
}
} // namespace Gascoigne

/*---------------------------------------------------*/
