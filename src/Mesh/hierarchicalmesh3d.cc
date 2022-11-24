/**
 *
 * Copyright (C) 2004, 2005, 2006, 2008, 2010, 2011 by the Gascoigne 3D authors
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

#include "hierarchicalmesh3d.h"

#include <fstream>
#include <numeric>

#include "../Common/set2vec.h"
#include "../Common/stlio.h"
#include "../Common/vecalgo.h"

#include "coarsehierarchicalmesh3d.h"
#include "deletecells.h"
#include "facemanager.h"
#include "levelcomparer3d.h"
#include "regular_update.h"

using namespace std;

namespace Gascoigne {
typedef triple<IndexType, IndexType, IndexType> tint;

/*------------------------------------------------------*/

HierarchicalMesh3d::HierarchicalMesh3d()
  : HierarchicalMesh()
  , HexLaO(hexs)
{}

/*------------------------------------------------------*/

HierarchicalMesh3d::HierarchicalMesh3d(const ParamFile& paramfile)
  : HierarchicalMesh()
  , HexLaO(hexs)
{
  BasicInit(paramfile);
}

/*------------------------------------------------------*/

HierarchicalMesh3d::HierarchicalMesh3d(const HierarchicalMesh3d& H)
  : HierarchicalMesh()
  , HexLaO(hexs)
{
  *this = H;
}

/*------------------------------------------------------*/

HierarchicalMesh3d&
HierarchicalMesh3d::operator=(const HierarchicalMesh3d& H)
{
  HierarchicalMesh::operator=(H);
  // copy all data
  vertexs3d = H.vertex3d();
  hexs = H.hex();
  Bquads = H.bquad();
  LineHang = H.linehang();
  QuadHang = H.quadhang();
  hexofcurved = H.GetHexOfCurved();

  return *this;
}

/*------------------------------------------------------*/

pair<IndexType, IndexType>
HierarchicalMesh3d::GetBoundaryInformation(IndexType i) const
{
  IndexType material = -1;
  IndexType le = -1;
  IndexType ib = GetBoundaryCellOfCurved(i);
  if (ib >= 0) {
    material = bquad(ib).material();
    le = bquad(ib).edge_in_quad();
  }
  return make_pair(material, le);
}

/*------------------------------------------------------*/

const BoundaryFunction<3>*
HierarchicalMesh3d::quad_shape(IndexType i) const
{
  if (GetCurvedShapes().empty())
    return NULL;

  if (GetCurvedShapes().Curved(i))
    return &(GetCurvedShapes().GetShape(i));

  return 0;
}

/*------------------------------------------------------*/

set<IndexType>
HierarchicalMesh3d::GetColors() const
{
  set<IndexType> coleur;

  for (IndexType i = 0; i < nbquads(); i++) {
    coleur.insert(bquad(i).material());
  }
  return coleur;
}

/*------------------------------------------------------*/

IndexType
HierarchicalMesh3d::FindPatchDepth() const
{
  // simple version, sucht nur p=1, p=0
  for (IndexType i = 0; i < ncells(); i++) {
    const Hex& q = hex(i);
    if (q.sleep())
      continue;
    IndexType father = q.father();
    if (father == -1)
      return 0;
    const Hex& qf = hex(father);
    for (IndexType ii = 0; ii < 8; ii++) {
      if (hex(qf.child(ii)).sleep())
        return 0;
    }
  }
  return 1;
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::InitHexOfCurved()
{
  hexofcurved.clear();

  if (GetCurvedShapes().empty())
    return;

  for (IndexType il = 0; il < nbquads(); ++il) {
    const BoundaryQuad& B = bquad(il);
    if (GetCurvedShapes().Curved(B.material())) {
      IndexType iq = B.of_quad();
      hexofcurved.insert(make_pair(iq, il));
    }
  }
}

/*------------------------------------------------------*/

void
HierarchicalMesh3d::prepare3d(const IndexVector& cell_ref,
                              const IndexVector& cell_coarse,
                              IndexSet& CellRefList,
                              IndexSet& CellCoarseList)
{
  /* copies cell_ref into CellRefList without duplets  */

  for (IndexVector::const_iterator cp = cell_ref.begin(); cp != cell_ref.end();
       ++cp) {
    IndexType c = *cp;
    if ((c >= 0) && (c < hexs.size())) {
      if (!hexs[c].sleep()) {
        CellRefList.insert(c);
      }
    }
  }

  /* copies cell_coarse into CellCoarseList without duplets
     checks if coarse cell in refine list */

  IndexSet help;

  for (IndexType i = 0; i < cell_coarse.size(); i++) {
    IndexType ic = cell_coarse[i];
    if ((ic < 0) || (ic >= hexs.size()))
      continue;

    if (hexs[ic].sleep())
      continue;
    if (!hexs[ic].level())
      continue;
    if (CellRefList.find(ic) != CellRefList.end())
      continue;

    help.insert(ic);
  }

  /* checks if coarse cell is cneighbour of a hang */

  for (HangList<4>::const_iterator Lp = QuadHang.begin(); Lp != QuadHang.end();
       Lp++) {
    IndexType cn = Lp->second.cneighbour();
    if (help.find(cn) != help.end())
      help.erase(cn);
  }

  /* marks the father */

  multiset<IndexType> ff;

  for (IndexSet::const_iterator hp = help.begin(); hp != help.end(); ++hp) {
    ff.insert(hex(*hp).father());
  }

  for (multiset<IndexType>::iterator fp = ff.begin(); fp != ff.end(); ++fp) {
    if (ff.count(*fp) == 8) {
      CellCoarseList.insert(*fp);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::FaceCoarse(HangContainer3d& hangset,
                               const IndexSet& cellcoarse) const
{
  std::array<IndexType, 4> quadglob;
  for (IntSetIt cp = cellcoarse.begin(); cp != cellcoarse.end(); ++cp) {
    IndexType f = *cp;
    for (IndexType face = 0; face < 6; ++face) {
      HexLaO.global_face_unsorted(quadglob, hex(f), face);

      IndexType ve = HexLaO.face_vertex(hex(f), face);
      hangset.face_coarse(quadglob, f, ve);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::FaceRefine(HangContainer3d& hangset,
                               const IndexSet& cellref) const
{
  std::array<IndexType, 4> quadglob;
  for (IntSetIt cp = cellref.begin(); cp != cellref.end(); ++cp) {
    IndexType f = *cp;
    for (IndexType face = 0; face < 6; ++face) {
      HexLaO.global_face_unsorted(quadglob, hex(f), face);
      hangset.face_refine(quadglob, f);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::UpdateHangingEdges(HangContainer3d& hangset,
                                       const IndexSet& cellref,
                                       const IndexSet& cellcoarse) const
{
  HangList<2> oldhangs(LineHang);
  hangset.build_hanging_lines(oldhangs);
  std::array<IndexType, 2> lineglob;

  for (IntSetIt cp = cellcoarse.begin(); cp != cellcoarse.end(); ++cp) {
    IndexType f = *cp;
    for (IndexType edge = 0; edge < 12; ++edge) {
      HexLaO.global_edge_unsorted(lineglob, hex(f), edge);

      IndexType ve = HexLaO.edge_vertex(hex(f), edge);
      hangset.line_coarse(lineglob, f, ve);
    }
  }
  for (IntSetIt cp = cellref.begin(); cp != cellref.end(); ++cp) {
    IndexType f = *cp;
    for (IndexType edge = 0; edge < 12; ++edge) {
      HexLaO.global_edge_unsorted(lineglob, hex(f), edge);

      hangset.line_refine(lineglob, f, oldhangs);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::UpdateHangs(HangContainer3d& hangset,
                                const IndexSet& cellref,
                                const IndexSet& cellcoarse)
{
  // faces
  FaceCoarse(hangset, cellcoarse);
  FaceRefine(hangset, cellref);
  ghost_fill_neighbours3d();

  // edges
  UpdateHangingEdges(hangset, cellref, cellcoarse);
  //  ghost_fill_neighbours2d();
  hangset.make_consistent();
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::ghost_fill_neighbours3d()
{
  FaceVector quadglob;
  for (IndexType ic = 0; ic < hexs.size(); ic++) {
    for (IndexType i = 0; i < 6; i++) {
      HexLaO.global_face_unsorted(quadglob, hex(ic), i);
      HangList<4>::iterator p = QuadHang.find(quadglob);

      if (p != QuadHang.end()) {
        IndexType cn = p->second.cneighbour();
        IndexType rn = p->second.rneighbour();

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
HierarchicalMesh3d::ghost_fill_neighbours2d()
{
  EdgeVector lineglob;
  for (IndexType ic = 0; ic < hexs.size(); ic++) {
    for (IndexType i = 0; i < 12; i++) {
      HexLaO.global_edge_unsorted(lineglob, hex(ic), i);
      LineHangList::iterator p = LineHang.find(lineglob);

      if (p != LineHang.end()) {
        IndexType cn = p->second.cneighbour();
        IndexType rn = p->second.rneighbour();

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
HierarchicalMesh3d::basic_refine3d(HangContainer3d& hangset,
                                   const IndexSet& CellRefList,
                                   const IndexSet& CellCoarseList)
{
  IndexType ov = nnodes();
  IndexType oc = hexs.size();
  IndexType oe = edges.size();

  IndexType csub = 8 * CellCoarseList.size();
  IndexType cadd = 8 * CellRefList.size();

  IndexType vsub = hangset.nDel() + CellCoarseList.size();
  IndexType vadd = hangset.nNew() + CellRefList.size();

  IndexType cdiff = cadd - csub;
  IndexType vdiff = vadd - vsub;

  IndexType nv = ov + vdiff;
  IndexType nc = oc + cdiff;

  clear_transfer_lists();

  /* delete vertex */

  IndexVector vdel, cdel, edel;

  for (IntSetCIt p = CellCoarseList.begin(); p != CellCoarseList.end(); p++) {
    const Hex& q = hex(*p);
    vdel.push_back(HexLaO.middle_vertex(q));
    cdel.insert(cdel.end(), q.childs().begin(), q.childs().end());
  }
  hangset.load_elimination(vdel);
  hangset.NeighbourSwapper();

  //////
  FaceManager EM(edges, hexs, co2n, eo2n);
  //////

  if (withfaces)
    EM.LoadFaceElimination(edel, CellCoarseList, hangset);

  IndexSet QuadRefList, QuadCoarseList, ccdel;

  boundary_prepare3d(QuadRefList, QuadCoarseList, ccdel, hangset);

  transfer(oc, co2n, cdel);
  transfer(ov, vo2n, vdel);
  transfer(oe, eo2n, edel);

  IndexVector cnew(cadd), vnew(vadd);
  iota(cnew.begin(), cnew.end(), oc - csub);
  iota(vnew.begin(), vnew.end(), ov - vsub);

  delete_vertexs3d(vo2n);

  if (withfaces)
    EM.DeleteFaces();

  delete_cells<Hex>(CellCoarseList, hexs, co2n, vo2n);

  vertexs3d.reserve(nv);
  vertexs3d.resize(nv);
  hexs.reserve(nc);
  hexs.resize(nc);

  hangset.update_olds(vo2n, co2n);
  hangset.update_news(vnew, CellRefList.size());

  new_vertexs3d(hangset, vnew, CellRefList);
  new_hexs(hangset, cnew, vnew, ov, CellRefList);

  if (withfaces)
    EM.Check(hangset);

  basic_fill_neighbours3d();

  new_boundary3d(QuadRefList, QuadCoarseList, ccdel);

  //   Lukas ?????
  //
  IndexSet adjustvertex;
  //   write_inp("refined0.inp");
  boundary_newton3d(adjustvertex);
  //   write_inp("refined1.inp");
  //   boundary_newton3d(adjustvertex);
  //   write_inp("refined2.inp");
  //   boundary_newton3d(adjustvertex);
  //   write_inp("refined3.inp");
  //   boundary_newton3d(adjustvertex);
  //   write_inp("refined4.inp");
  //   boundary_newton3d(adjustvertex);
  //   write_inp("refined5.inp");
  //   boundary_newton3d(adjustvertex);
  //   write_inp("refined6.inp");
  inner_vertex_newton3d(vnew, CellRefList, adjustvertex);

  if (withfaces)
    EM.Build(CellRefList, hangset);
  Testing();
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::Testing()
{
  for (IndexType i = 0; i < edges.size(); i++) {
    IndexType m = edges[i].master();
    IndexType s = edges[i].slave();

    if (s < 0)
      continue;
    if (!hex(m).sleep() && hex(s).sleep()) {
      swap(edges[i].master(), edges[i].slave());
      swap(edges[i].LocalMasterIndex(), edges[i].LocalSlaveIndex());
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::boundary_newton3d(IndexSet& adjustvertex)
{
  if (GetCurvedShapes().empty())
    return;

  for (IndexType i = 0; i < Bquads.size(); i++) {
    BoundaryQuad& bl = Bquads[i];
    IndexType color = bl.material();

    if (GetCurvedShapes().Curved(color)) {
      for (IndexType j = 0; j < 4; j++) {
        IndexType k = bl.vertex(j);

        GetCurvedShapes().newton(color, vertexs3d[k]);
        adjustvertex.insert(k);
      }
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::boundary_prepare3d(IndexSet& QuadRefList,
                                       IndexSet& QuadCoarseList,
                                       IndexSet& ccdel,
                                       const HangContainer3d& hangset)
{
  for (IndexType i = 0; i < Bquads.size(); i++) {
    const BoundaryQuad& bl = Bquads[i];

    FaceVector lineglob;
    lineglob[0] = bl.vertex(0);
    lineglob[1] = bl.vertex(1);
    lineglob[2] = bl.vertex(2);
    lineglob[3] = bl.vertex(3);
    sort(lineglob.begin(), lineglob.end());

    if (bl.sleep()) {
      if (hangset.ToBeDeleted(lineglob)) {
        QuadCoarseList.insert(i);
        for (IndexType j = 0; j < 4; j++)
          ccdel.insert(bl.child(j));
      }
    } else {
      if (hangset.ToBeCreated(lineglob)) {
        QuadRefList.insert(i);
      }
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::new_boundary3d(IndexSet& ref,
                                   IndexSet& coarse,
                                   IndexSet& ccdel)
{
  IndexType oc = Bquads.size();
  IndexType csub = 4 * coarse.size();
  IndexType cadd = 4 * ref.size();
  IndexType nc = oc + cadd - csub;

  IndexVector lo2n;
  transfer(oc, lo2n, ccdel);
  delete_cells<BoundaryQuad>(coarse, Bquads, lo2n, vo2n);

  update_boundary_data3d(coarse);

  IndexVector cnew(cadd);
  iota(cnew.begin(), cnew.end(), oc - csub);

  Bquads.reserve(nc);
  Bquads.resize(nc);

  new_bquads(lo2n, cnew, ref);
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::update_boundary_data3d(const IndexSet& LCoarse)
{
  IndexType no = Bquads.size() - 4 * LCoarse.size();
  for (IndexType i = 0; i < no; ++i) {
    IndexType oq = Bquads[i].of_quad();

    assert(co2n[oq] >= 0);

    Bquads[i].of_quad() = co2n[oq];
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::new_bquads(const IndexVector& lo2n,
                               const IndexVector& cnew,
                               const IndexSet& LRef)
{
  IndexType nci = 0;

  vector<Quad> emptyq;

  for (IntSetIt cp = LRef.begin(); cp != LRef.end(); ++cp) {
    IndexType father = lo2n[*cp];

    assert(father != -1);

    BoundaryQuad& bqf = Bquads[father];
    // change father boundary
    vector<IndexType>& qc = bqf.childs();
    qc.resize(4);

    IndexType face = bqf.edge_in_quad(); // face index in hex
    IndexType hexi = bqf.of_quad();      // hex index of father
    FaceVector chvec;
    HexLaO.childs_of_face(chvec, hex(hexi), face);

    for (IndexType ic = 0; ic < 4; ic++) {
      IndexType childhex = chvec[ic];
      IndexType inold = cnew[nci + ic];
      // set childs in father
      qc[ic] = inold;
      // set properties of childs
      BoundaryQuad& bl = Bquads[inold];
      bl.level() = bqf.level() + 1;
      bl.father() = father;
      bl.material() = bqf.material();
      bl.childs().resize(0);

      HexLaO.load_face(bl.vertex(), hex(childhex), face);
      bl.of_quad() = childhex;
      bl.edge_in_quad() = face;
    }
    nci += 4;
  }
  // nboundary += nci;
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::new_vertexs3d(HangContainer3d& hangset,
                                  const IndexVector& vnew,
                                  const IndexSet& CellRefList)
{
  IndexType nv1 = CellRefList.size();

  IntSetIt cp = CellRefList.begin();
  for (IndexType i = 0; i < nv1; i++) {
    IndexType f = co2n[*cp++];
    new_middle_vertex3d(vnew[i], f);
  }
  /* new vertexes on faces */
  HangList<4>::const_iterator f = hangset.FaceCreating().begin();

  for (; f != hangset.FaceCreating().end(); f++) {
    new_face_vertex3d(f->second.hanging(), f->first);
  }
  /* new vertexes on edges */
  HangList<2>::const_iterator p = hangset.Creating().begin();

  for (; p != hangset.Creating().end(); p++) {
    new_edge_vertex3d(p->second.hanging(), p->first);
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::new_middle_vertex3d(IndexType i, IndexType f)
{
  IndexType v0 = hexs[f].vertex(0);
  IndexType v1 = hexs[f].vertex(1);
  IndexType v2 = hexs[f].vertex(2);
  IndexType v3 = hexs[f].vertex(3);
  IndexType v4 = hexs[f].vertex(4);
  IndexType v5 = hexs[f].vertex(5);
  IndexType v6 = hexs[f].vertex(6);
  IndexType v7 = hexs[f].vertex(7);

  Vertex3d w1, w2;

  w1.equ(0.25,
         vertexs3d[v0],
         0.25,
         vertexs3d[v1],
         0.25,
         vertexs3d[v2],
         0.25,
         vertexs3d[v3]);
  w2.equ(0.25,
         vertexs3d[v4],
         0.25,
         vertexs3d[v5],
         0.25,
         vertexs3d[v6],
         0.25,
         vertexs3d[v7]);

  vertexs3d[i].equ(0.5, w1, 0.5, w2);
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::new_hexs(const HangContainer3d& hangset,
                             const IndexVector& cnew,
                             const IndexVector& vnew,
                             IndexType nvold,
                             const IndexSet& CellRefList)
{
  // neue zellen erzeugen
  // eintragen der "Vater-Vertexs" in den kindern

  IndexType nci = 0;
  IndexType ivm = 0;

  for (IntSetCIt cp = CellRefList.begin(); cp != CellRefList.end(); cp++) {
    IndexType father = co2n[*cp];

    assert(father != -1);

    vector<IndexType>& qc = hexs[father].childs();
    qc.resize(8);
    IndexType material = hexs[father].material();
    IndexType material_Vanka = hexs[father].material_Vanka();
    std::array<Vertex3d, 3> basis_Vanka = hexs[father].basis_Vanka();

    IndexType childlevel = hexs[father].level() + 1;
    for (IndexType ic = 0; ic < 8; ic++) {
      IndexType inold = cnew[nci + ic];
      qc[ic] = inold;
      hexs[inold].level() = childlevel;
      hexs[inold].father() = father;
      hexs[inold].material() = material;
      hexs[inold].material_Vanka() = material_Vanka;
      hexs[inold].basis_Vanka() = basis_Vanka;
      hexs[inold].childs().resize(0);
      hexs[inold].edges().fill(-1);
    }
    HexLaO.fill_corner_vertex_in_childs(hexs[father]);

    HexLaO.fill_middle_vertex_in_childs(hexs[father], vnew[ivm]);
    ivm++;
    nci += 8;

    // Edge Vertex -- linehanginfo schon ok (hanging) !
    IndexType ive(-1);
    FaceVector faceglob;
    for (IndexType i = 0; i < 6; i++) {
      HexLaO.global_face_unsorted(faceglob, hex(father), i);

      ive = hangset.vertex_index(faceglob);

      assert(ive >= 0);

      HexLaO.fill_face_vertex_in_childs(hexs[father], i, ive);
    }
    std::array<IndexType, 2> lineglob;
    for (IndexType i = 0; i < 12; i++) {
      HexLaO.global_edge_unsorted(lineglob, hex(father), i);

      ive = hangset.vertex_index(lineglob);

      if (ive < 0) {
        cout << "Line " << lineglob << endl;
        cout << "new vertex " << ive << endl;
        cout << "father hex " << father << endl;
      }

      assert(ive >= 0);

      HexLaO.fill_edge_vertex_in_childs(hexs[father], i, ive);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::inner_vertex_newton3d(const IndexVector& vnew,
                                          const IndexSet& CellRefList,
                                          const IndexSet& adjustvertex)
{
  if (GetCurvedShapes().empty())
    return;

  // baue lokalen set auf um spaeter nicht alle hex zu justieren
  IndexSet Hexset;
  for (IndexType i = 0; i < Bquads.size(); i++) {
    IndexType hi = Bquads[i].of_quad();
    Hexset.insert(hi);
  }

  IntSetIt cp = CellRefList.begin();

  for (IndexType i = 0; i < CellRefList.size(); i++) {
    IndexType hi = co2n[*cp++];
    const Hex& h = hex(hi);

    if (Hexset.find(hi) == Hexset.end())
      continue;

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEU
    // edges
    for (IndexType e = 0; e < 12; ++e) {
      IndexType ev = HexLaO.edge_vertex(h, e);
      if (adjustvertex.find(ev) != adjustvertex.end())
        continue;

      std::array<IndexType, 2> fe;
      HexLaO.globalvertices_of_edge(h, fe, e);
      vertexs3d[ev] = vertexs3d[fe[0]];
      vertexs3d[ev] += vertexs3d[fe[1]];
      vertexs3d[ev] *= 0.5;
    }
    // faces
    for (IndexType f = 0; f < 6; ++f) {
      IndexType fv = HexLaO.face_vertex(h, f);
      if (adjustvertex.find(fv) != adjustvertex.end())
        continue;

      std::array<IndexType, 4> fe;
      HexLaO.LoadEdgeVerticesOfFace(h, f, fe);
      vertexs3d[fv].equ(0.25,
                        vertexs3d[fe[0]],
                        0.25,
                        vertexs3d[fe[1]],
                        0.25,
                        vertexs3d[fe[2]],
                        0.25,
                        vertexs3d[fe[3]]);
    }
    // middle
    /*     IndexType cv = HexLaO.middle_vertex(h);
         assert (adjustvertex.find(cv)==adjustvertex.end());

         // weighted sum of edge vertices
         double edgewgt = 0.3;
         // .. face vertices
         double facewgt = 0.3;
         // ... and nodes
         double nodewgt = 0.4;
         assert(fabs(edgewgt + facewgt + nodewgt - 1.0)<1.e-10);

         // All vertices of the patch
         vector<IndexType> pv;
         pv.push_back(h.vertex(0));
         pv.push_back(HexLaO.edge_vertex(h,0));
         pv.push_back(h.vertex(1));
         pv.push_back(HexLaO.edge_vertex(h,3));
         pv.push_back(HexLaO.face_vertex(h,0));
         pv.push_back(HexLaO.edge_vertex(h,1));
         pv.push_back(h.vertex(3));
         pv.push_back(HexLaO.edge_vertex(h,2));
         pv.push_back(h.vertex(2));
         //
         pv.push_back(HexLaO.edge_vertex(h,8));
         pv.push_back(HexLaO.face_vertex(h,4));
         pv.push_back(HexLaO.edge_vertex(h,9));
         pv.push_back(HexLaO.face_vertex(h,3));
         pv.push_back(HexLaO.middle_vertex(h));
         pv.push_back(HexLaO.face_vertex(h,1));
         pv.push_back(HexLaO.edge_vertex(h,11));
         pv.push_back(HexLaO.face_vertex(h,2));
         pv.push_back(HexLaO.edge_vertex(h,10));
         //
         pv.push_back(h.vertex(4));
         pv.push_back(HexLaO.edge_vertex(h,4));
         pv.push_back(h.vertex(5));
         pv.push_back(HexLaO.edge_vertex(h,7));
         pv.push_back(HexLaO.face_vertex(h,5));
         pv.push_back(HexLaO.edge_vertex(h,5));
         pv.push_back(h.vertex(7));
         pv.push_back(HexLaO.edge_vertex(h,6));
         pv.push_back(h.vertex(6));
         assert(pv.size()==27);

         // measure hx/hy/hz in local coordinate system
         double hx = (vertexs3d[pv[14]]-vertexs3d[pv[12]]).norm();
         double hy = (vertexs3d[pv[16]]-vertexs3d[pv[10]]).norm();
         double hz = (vertexs3d[pv[22]]-vertexs3d[pv[4]]).norm();
         assert(hz>hx);
         assert(hz>hy);


         for (IndexType z=0;z<3;++z)
           {
             IndexType i0=9*z;

             // edges orthogonal to z
             IndexType ed[4][3]= {{1,0,2},{3,0,6},{5,2,8},{7,6,8}};
             for (IndexType e=0;e<4;++e)
               {
                 IndexType ev = pv[i0+ed[e][0]];
                 if (adjustvertex.find(ev)!=adjustvertex.end()) continue;
                 vertexs3d[ev].equ(0.5,vertexs3d[pv[i0+ed[e][1]]],0.5,vertexs3d[pv[i0+ed[e][2]]]);
               }

             IndexType mv = pv[i0+4];
             if (adjustvertex.find(mv)!=adjustvertex.end()) continue;
             vertexs3d[mv].equ(0.25,vertexs3d[pv[i0+1]],
                               0.25,vertexs3d[pv[i0+3]],
                               0.25,vertexs3d[pv[i0+5]],
                               0.25,vertexs3d[pv[i0+7]]);
           }
      */

    IndexType fv = HexLaO.middle_vertex(h);
    assert(adjustvertex.find(fv) == adjustvertex.end());
    std::array<IndexType, 6> fe;
    HexLaO.LoadFaceVertices(h, fe);
    vertexs3d[fv] = 0;
    for (IndexType i = 0; i < 6; ++i)
      vertexs3d[fv] += vertexs3d[fe[i]];
    vertexs3d[fv] *= 1. / 6.;
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEU

    //       // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ALT
    //       for (IndexType face=0; face<6; face++)
    // 	{
    // 	  HexLaO.LoadEdgeVerticesOfFace(h,face,v);
    // 	  IndexType mv = HexLaO.face_vertex(h,face);
    // 	  if (adjustvertex.find(mv)!=adjustvertex.end()) continue;
    // 	  new_face_vertex3d(mv,v);
    // 	}
    // //       std::array<IndexType,6> w;
    // //       IndexType mv = HexLaO.middle_vertex(h);
    // //       HexLaO.LoadFaceVertices(h,w);
    // //       new_vertex3d(mv,w);
    //       // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ALT
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::basic_fill_neighbours3d()
{
  IndexType n = 0;
  FaceVector faceglob;
  for (IndexType ic = 0; ic < hexs.size(); ic++) {
    for (IndexType i = 0; i < 6; i++) {
      HexLaO.global_face_unsorted(faceglob, hex(ic), i);
      HangList<4>::iterator p = QuadHang.find(faceglob);

      if (p != QuadHang.end()) {
        IndexType cn = p->second.cneighbour();
        IndexType rn = p->second.rneighbour();

        if ((cn == -1) && (!hexs[ic].sleep())) {
          n++;
          p->second.cneighbour() = ic;
        } else if (rn == -1) {
          n++;
          p->second.rneighbour() = ic;
        }
      }
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::init_quad(BoundaryQuad& newquad)
{
  FaceVector v;
  IndexSet newquad_set;
  newquad_set.insert(newquad[0]);
  newquad_set.insert(newquad[1]);
  newquad_set.insert(newquad[2]);
  newquad_set.insert(newquad[3]);
  for (IndexType i = 0; i < hexs.size(); i++) {
    for (IndexType face = 0; face < 6; face++) {
      HexLaO.global_face_unsorted(v, hex(i), face);
      IndexSet v_set;
      v_set.insert(v[0]);
      v_set.insert(v[1]);
      v_set.insert(v[2]);
      v_set.insert(v[3]);
      if (newquad_set == v_set) {
        HexLaO.global_face_unsorted(v, hex(i), face);
        newquad.vertex() = v;
        newquad.of_quad() = i;
        newquad.edge_in_quad() = face;
        return;
      }
    }
  }
  cerr << "BoundaryQuad not found !" << endl;
  cerr << newquad;
  abort();
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::init_edges3d()
{
  if (withfaces) {
    FaceManager EM(edges, hexs, co2n, eo2n);
    EM.InitFaces();
    EM.SortHangings();
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::GetVertexesOfFace(std::array<IndexType, 5>& v,
                                      IndexType e) const
{
  const Edge& E = edge(e);
  const Hex* Q = &hex(E.master());

  IndexType le = E.LocalMasterIndex();
  v[0] = (*Q)[le];
  v[1] = (*Q)[(le + 1) % 8];
  v[2] = (*Q)[(le + 2) % 8];
  v[3] = (*Q)[(le + 3) % 8];
  v[4] = -1;

  if (!Q->sleep()) {
    if (E.slave() == -1)
      return;
    Q = &hex(E.slave());
    le = E.LocalSlaveIndex();
  }
  v[4] = HexLaO.face_vertex(*Q, le);
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::GetVertexesOfFace(std::array<IndexType, 4>& v,
                                      IndexType e) const
{
  const Edge& E = edge(e);
  const Hex* Q = &hex(E.master());

  IndexType le = E.LocalMasterIndex();
  v[0] = (*Q)[le];
  v[1] = (*Q)[(le + 1) % 8];
  v[2] = (*Q)[(le + 2) % 8];
  v[3] = (*Q)[(le + 3) % 8];
}

/*------------------------------------------------------*/

IndexType
HierarchicalMesh3d::Vater(const IndexType i) const
{
  return hex(i).father();
}

/*------------------------------------------------------*/

IndexType
HierarchicalMesh3d::nactivedescendants(IndexType i) const
{
  if (!hexs[i].sleep())
    return 1;
  IndexType k = 0;
  for (IndexType j = 0; j < hexs[i].nchilds(); ++j)
    k += nactivedescendants(hexs[i].child(j));
  return k;
}

/*------------------------------------------------------*/

IndexVector
HierarchicalMesh3d::GetVertices(IndexType c) const
{
  IndexVector v;
  for (IndexType i = 0; i < hexs[c].nvertexs(); ++i)
    v.push_back(hexs[c][i]);
  return v;
}

/*------------------------------------------------------*/

IndexVector
HierarchicalMesh3d::Nachkommen(const IndexType i) const
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

IndexVector
HierarchicalMesh3d::Kinder(const IndexType i) const
{
  IndexVector k = hex(i).childs();
  return k;
}

/*------------------------------------------------------*/

IndexVector
HierarchicalMesh3d::Geschwister(const IndexType i) const
{
  const Hex& q = hex(i);
  IndexType father = q.father();
  if (father == -1) {
    IndexVector n(1, i);
    return n;
  }
  return Kinder(father);
}

/*------------------------------------------------------*/

std::array<IndexType, 4>
HierarchicalMesh3d::ChildrenOfFace(IndexType e) const
{
  IndexType s = edge(e).slave();
  IndexType is = edge(e).LocalSlaveIndex();

  assert(s >= 0);

  std::array<IndexType, 4> f;
  for (IndexType ii = 0; ii < 4; ii++) {
    IndexType ic = hex(s).child(HexLaO.ChildsOfFace(is, ii));
    f[ii] = hex(ic).edge(HexLaO.ChildFace(is));
  }
  return f;
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::GetAwakeCells(set<IndexType>& v) const
{
  for (IndexType i = 0; i < ncells(); i++) {
    if (!hexs[i].sleep())
      v.insert(i);
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::GetAwakePatchs(set<IndexType>& v) const
{
  for (IndexType i = 0; i < ncells(); i++) {
    if (!hexs[i].sleep()) {
      IndexType f = hexs[i].father();
      assert(f != -1);
      v.insert(f);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::ConstructQ2PatchMesh(IndexVector& q2patchmesh) const
{
  typedef set<IndexType>::iterator It;
  set<IndexType> patche;
  GetAwakePatchs(patche);
  vector<set<IndexType>> patch_on_level(nlevels());
  q2patchmesh.resize(0);
  for (It it = patche.begin(); it != patche.end(); ++it)
    patch_on_level[hexs[*it].level()].insert(*it);
  // grobgitterzellen kommen ins patchmesh
  for (It it = patch_on_level[0].begin(); it != patch_on_level[0].end(); ++it)
    q2patchmesh.push_back(*it);
  // der Rest wird eins groeber
  for (IndexType l = 1; l < nlevels(); ++l) {
    It it = patch_on_level[l].begin();
    while (it != patch_on_level[l].end()) {
      IndexType v = hexs[*it].father();
      assert(v != -1);
      q2patchmesh.push_back(v);
      IndexVector nk = Nachkommen(v);
      for (IndexType i = 0; i < nk.size(); ++i)
        patch_on_level[hexs[nk[i]].level()].erase(nk[i]);
      it = patch_on_level[l].begin();
    }
  }
}

/*---------------------------------------------------*/

IndexVector
HierarchicalMesh3d::ConstructQ4Patch(IndexType c) const
{
  IndexVector patch(125, -1);
  for (IndexType i = 0; i < 125; i++) {
    // Vertex i steht an Position (x,y,z)
    IndexType x = i % 5;
    IndexType y = (i % 25) / 5;
    IndexType z = i / 25;

    // Position von erstem Kind
    IndexType fcx = x / 3;
    IndexType fcy = y / 3;
    IndexType fcz = z / 3;
    // Index davon
    IndexType fci =
      fcz * 4 + fcy * 2 + abs(static_cast<long>(fcx) - static_cast<long>(fcy));

    // Position von Kind im Kind
    IndexType scx = (x - 2 * fcx) / 2;
    IndexType scy = (y - 2 * fcy) / 2;
    IndexType scz = (z - 2 * fcz) / 2;
    // Index davon
    IndexType sci =
      scz * 4 + scy * 2 + abs(static_cast<long>(scx) - static_cast<long>(scy));

    // Position des Vertex
    IndexType vx = x - 2 * fcx - scx;
    IndexType vy = y - 2 * fcy - scy;
    IndexType vz = z - 2 * fcz - scz;
    // Index davon
    IndexType vi =
      vz * 4 + vy * 2 + abs(static_cast<long>(vx) - static_cast<long>(vy));

    patch[i] = hexs[hexs[hexs[c].child(fci)].child(sci)].vertex(vi);
  }
  return patch;
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh3d::NodeOnFace(IndexType e) const
{
  // only for real hanging nodes
  IndexType m = edge(e).master();
  IndexType im = edge(e).LocalMasterIndex();
  IndexType s = edge(e).slave();
  IndexType is = edge(e).LocalSlaveIndex();

  if (s < 0) {
    assert(hex(im).sleep());
    return HexLaO.face_vertex(hex(m), im);
  }

  return HexLaO.face_vertex(hex(s), is);
}

/*------------------------------------------------------*/

IndexType
HierarchicalMesh3d::neighbour(IndexType c, IndexType le) const
{
  const Hex& Q = hex(c);
  const Edge& E = edge(Q.edge(le));
  IndexType m = E.master();
  IndexType nq = m;
  if (m == c)
    nq = E.slave();
  return nq;
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh3d::neighbour_neighbour(IndexType c, IndexType le) const
{
  assert(le < 6);
  IndexType n = neighbour(c, le);
  assert(n >= 0);

  IndexType nn = 0;
  for (nn = 0; nn < 6; ++nn)
    if (c == neighbour(n, nn))
      break;
  assert(nn < 6);
  return nn;
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::FillAllBoundaryLines()
{
  cout << "HierarchicalMesh3d::FillAllBoundaryLines() not written" << endl;
  abort();
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::WriteAll(const string& name) const
{
  ofstream out(name.c_str());

  out << dimension() << " dimension" << endl;
  out << nnodes() << " vertexs" << endl;
  out << mnlevels << " mnlevels" << endl;

  cerr << "HierarchicalMesh3d::WriteAll()\n";
  cerr << "not written 3d\n";

  out.close();
  abort();
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::write_inp(const string& name) const
{
  ofstream file(name.c_str());
  if (!file.is_open()) {
    cerr << "HierarchicalMesh3d::write_inp()\n";
    cerr << "cannot open file " << name << endl;
    abort();
  }

  IndexType nt = ncells() + nbquads();
  file << nnodes() << " " << nt << " " << 0 << " " << 0 << " " << 0 << endl;

  for (IndexType i = 0; i < nnodes(); i++)
    file << i << " " << vertex3d(i) << " " << endl;

  for (IndexType i = 0; i < ncells(); i++) {
    file << i << " " << 0 << " hex " << hex(i).vertex() << endl;
  }
  for (IndexType i = 0; i < nbquads(); i++) {
    file << i << " " << bquad(i).material() << " quad " << bquad(i).vertex()
         << endl;
  }
}

/*---------------------------------------------------*/

pair<bool, tint>
HierarchicalMesh3d::check_inp(const string& name)
{
  //  detect some errors in input-file... and more

  ifstream file(name.c_str());
  if (!file.is_open()) {
    cerr << "HierarchicalMesh3d::check_inp()\n";
    cerr << "cannot open file " << name << endl;
    abort();
  }

  bool first_one = 1;
  IndexType nv, nl, nq, nh, nt;
  IndexType n_unkonwn;
  file >> nv >> nt >> n_unkonwn >> n_unkonwn >> n_unkonwn;

  Vertex3d c;
  IndexType ind;
  for (IndexType i = 0; i < nv; i++) {
    file >> ind >> c;
  }

  nh = 0;
  nq = 0;
  nl = 0;
  std::array<IndexType, 8> ih;
  std::array<IndexType, 4> iq;
  std::array<IndexType, 2> il;
  for (IndexType i = 0; i < nt; i++) {
    string name;
    string mat;
    IndexType ii;
    file >> ii >> mat >> name;
    if (name == "hex") {
      file >> ih;
      nh++;
      if ((ih[0] == 0) || (ih[1] == 0) || (ih[2] == 0) || (ih[3] == 0) ||
          (ih[4] == 0) || (ih[5] == 0) || (ih[6] == 0) || (ih[7] == 0)) {
        first_one = 0;
      }
    } else if (name == "quad") {
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
  if (nt != (nl + nq + nh)) {
    cerr << "wrong number of cells: " << nt << endl;
    cerr << "lines quads hexs: " << nl << " " << nq << " " << nh << endl;
    abort();
  }

  file.close();

  return make_pair(first_one, make_triple(nl, nq, nh));
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::read_inp(const string& name)
{
  // check mesh.....

  pair<bool, tint> p = check_inp(name);
  bool first_one = p.first;
  tint n = p.second;

  IndexType nl = n.first;
  IndexType nq = n.second;
  IndexType nh = n.third;

  ifstream file(name.c_str());
  if (!file.is_open()) {
    cerr << "HierarchicalMesh3d::read_inp()\n";
    cerr << "cannot open file " << name << endl;
    abort();
  }

  IndexType nv, nt, n_unkonwn;

  file >> nv >> nt >> n_unkonwn >> n_unkonwn >> n_unkonwn;

  cout << "3D Mesh:  " << nv << " nodes, ";
  cout << nl << " lines, ";
  cout << nq << " quads, ";
  cout << nh << " hexs" << endl;

  vertexs3d.reserve(nv);
  vertexs3d.resize(nv, Vertex3d());
  hexs.reserve(nh);
  hexs.resize(nh, Hex());
  Bquads.reserve(nq);
  Bquads.resize(nq);

  Vertex3d c;
  IndexType ind;
  for (IndexType i = 0; i < nv; i++) {
    file >> ind >> c;
    vertexs3d[i] = c;
  }
  std::array<IndexType, 8> ihv;
  std::array<IndexType, 4> iqv;

  IndexType ih = 0;
  IndexType iq = 0;
  for (IndexType i = 0; i < nt; i++) {
    string name;
    IndexType unknown;
    string matstring;
    file >> unknown >> matstring >> name;
    if (name == "hex") {
      file >> ihv;
      if (first_one)
        for (IndexType iii = 0; iii < 8; iii++)
          ihv[iii]--;
      hexs[ih].vertex() = ihv;
      hexs[ih].material() = atoi(matstring.c_str());

      hexs[ih].material_Vanka() = ih;

      /*----------------------------------------------*/
      vector<pair<Vertex3d, double>> v_orient;
      Vertex3d v_zw = vertexs3d[ihv[1]];
      v_zw += vertexs3d[ihv[2]];
      v_zw += vertexs3d[ihv[5]];
      v_zw += vertexs3d[ihv[6]];
      v_zw -= vertexs3d[ihv[0]];
      v_zw -= vertexs3d[ihv[3]];
      v_zw -= vertexs3d[ihv[4]];
      v_zw -= vertexs3d[ihv[7]];
      // v_zw-=vertexs3d[ihv[0]];
      v_orient.push_back(make_pair(v_zw, v_zw.norm()));
      v_zw = vertexs3d[ihv[4]];
      v_zw += vertexs3d[ihv[5]];
      v_zw += vertexs3d[ihv[6]];
      v_zw += vertexs3d[ihv[7]];
      v_zw -= vertexs3d[ihv[0]];
      v_zw -= vertexs3d[ihv[1]];
      v_zw -= vertexs3d[ihv[2]];
      v_zw -= vertexs3d[ihv[3]];

      // v_zw=vertexs3d[ihv[4]]	;
      // v_zw-=vertexs3d[ihv[0]];
      v_orient.push_back(make_pair(v_zw, v_zw.norm()));
      v_zw = vertexs3d[ihv[2]];
      v_zw += vertexs3d[ihv[3]];
      v_zw += vertexs3d[ihv[6]];
      v_zw += vertexs3d[ihv[7]];
      v_zw -= vertexs3d[ihv[0]];
      v_zw -= vertexs3d[ihv[1]];
      v_zw -= vertexs3d[ihv[5]];
      v_zw -= vertexs3d[ihv[4]];

      // v_zw=vertexs3d[ihv[3]]	;
      // v_zw-=vertexs3d[ihv[0]];
      v_orient.push_back(make_pair(v_zw, v_zw.norm()));

      std::sort(v_orient.begin(), v_orient.end(), sort_pred());
      for (IndexType i = 0; i < 3; i++) {
        hexs[ih].basis_Vanka()[i] = v_orient[i].first;
      }
      /*----------------------------------------------*/
      ih++;
    } else if (name == "quad") {
      file >> iqv;
      if (first_one)
        for (IndexType iii = 0; iii < 4; iii++)
          iqv[iii]--;

      BoundaryQuad li;
      li.material() = atoi(matstring.c_str());
      li.vertex() = iqv;
      init_quad(li);
      Bquads[iq++] = li;
    }
  }
  init_edges3d();
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::write(const string& bname) const
{
  write_gup(bname);
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::write_gup(const string& bname) const
{
  string name = bname;
  name += ".gup";

  ofstream out(name.c_str());

  out << dimension() << " dimension" << endl;
  out << nnodes() << " vertexs" << endl;

  for (IndexType i = 0; i < nnodes(); i++) {
    out << " " << vertex3d(i) << endl;
  }
  out << hexs.size() << " hexs" << endl;
  for (IndexType i = 0; i < hexs.size(); i++) {
    out << hex(i);
  }
  out << QuadHang << endl;
  out << LineHang << endl;
  out << Bquads.size() << " boundaryquads" << endl;
  for (IndexType i = 0; i < Bquads.size(); i++) {
    out << Bquads[i].material() << " " << Bquads[i] << endl;
  }
  out << endl << edges.size() << " edges" << endl;
  for (IndexType i = 0; i < edges.size(); i++) {
    out << " " << edges[i];
  }
  out.close();
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::write_gip(const string& bname) const
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
    ArrayBinWrite(out, vertex3d(i));
  }

  IndexType nhexs = hexs.size();
  out.write(reinterpret_cast<const char*>(&nhexs), sizeInt);

  for (IndexType i = 0; i < hexs.size(); i++) {
    ArrayBinWrite(out, hex(i));
  }

  QuadHang.BinWrite(out);
  LineHang.BinWrite(out);
  IndexType nbquads = Bquads.size();
  out.write(reinterpret_cast<const char*>(&nbquads), sizeInt);
  for (IndexType i = 0; i < Bquads.size(); i++) {
    IndexType mat = Bquads[i].material();
    out.write(reinterpret_cast<const char*>(&mat), sizeInt);
    Bquads[i].BinWrite(out);
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
HierarchicalMesh3d::read_gup(const string& name)
{
  string symbol;

  ifstream file(name.c_str());

  assert(file.is_open());

  Bquads.clear();
  edges.clear();
  QuadHang.clear();
  LineHang.clear();

  IndexType n;
  IndexType dim;
  file >> dim >> symbol >> n >> symbol;

  assert(dim == 3);

  if (_i_showoutput) {
    cout << "Mesh 3d  : ";
  }
  vertexs3d.reserve(n);
  vertexs3d.resize(n);

  assert(symbol == "vertexs");

  for (IndexType i = 0; i < n; i++) {
    file >> vertexs3d[i];
  }
  if (_i_showoutput) {
    cout << n << " nodes, ";
  }

  file >> n >> symbol;

  assert(symbol == "hexs");

  hexs.reserve(n);
  hexs.resize(n);
  for (IndexType i = 0; i < hexs.size(); i++) {
    file >> hexs[i];
  }
  if (_i_showoutput) {
    cout << n << " hexs, ";
  }
  file >> QuadHang;
  if (_i_showoutput) {
    cout << QuadHang.size() << " quadhangs, ";
  }
  file >> LineHang;
  if (_i_showoutput) {
    cout << LineHang.size() << " linehangs, ";
  }

  file >> n >> symbol;
  IndexType number = 0;
  assert(symbol == "boundaryquads");

  BoundaryQuad bol;
  for (IndexType i = 0; i < n; i++) {
    file >> bol.material() >> bol;
    Bquads.push_back(bol);
  }
  number = n;

  if (_i_showoutput) {
    cout << number << " boundaryquads, ";
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
    init_edges3d();
  post_refine3d();
  file.close();
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::read_gip(const string& bname)
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
    cerr << "HierarchicalMesh3d::read_gip(): error in file " << name << endl;
    abort();
  }

  Bquads.clear();
  edges.clear();
  QuadHang.clear();
  LineHang.clear();

  IndexType n, dim, sizeInt = sizeof(IndexType);
  file.read(reinterpret_cast<char*>(&dim), sizeInt);
  file.read(reinterpret_cast<char*>(&n), sizeInt);

  assert(dim == 3);

  if (_i_showoutput) {
    cout << "Mesh 3d  : ";
  }
  vertexs3d.reserve(n);
  vertexs3d.resize(n);

  for (IndexType i = 0; i < n; i++) {
    ArrayBinRead(file, vertexs3d[i]);
  }
  if (_i_showoutput) {
    cout << n << " nodes, ";
  }
  file.read(reinterpret_cast<char*>(&n), sizeInt);
  hexs.reserve(n);
  hexs.resize(n);
  for (IndexType i = 0; i < hexs.size(); i++) {
    hexs[i].BinRead(file);
  }
  if (_i_showoutput) {
    cout << n << " hexs, ";
  }
  QuadHang.BinRead(file);
  if (_i_showoutput) {
    cout << QuadHang.size() << " quadhangs, ";
  }
  LineHang.BinRead(file);
  if (_i_showoutput) {
    cout << LineHang.size() << " linehangs, ";
  }
  file.read(reinterpret_cast<char*>(&n), sizeInt);
  IndexType number = 0;
  BoundaryQuad bol;
  for (IndexType i = 0; i < n; i++) {
    file.read(reinterpret_cast<char*>(&bol.material()), sizeInt);
    bol.BinRead(file);
    Bquads.push_back(bol);
  }
  number = n;
  if (_i_showoutput) {
    cout << number << " boundaryquads, ";
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
    init_edges3d();
  post_refine3d();
  file.close();
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh3d::regular_grid3d_one(IndexSet& celllist,
                                       IndexVector& coarsesub,
                                       const IndexSet& CellRefList,
                                       const IndexSet& CellCoarseList)
{
  /* detects jump over two levels across LineHangs */

  IndexType n = 0;
  HangList<4>::const_iterator hp;

  for (hp = QuadHang.begin(); hp != QuadHang.end(); ++hp) {
    IndexType cr = hp->second.rneighbour();
    IndexType cn = hp->second.cneighbour();

    assert(cr >= 0);

    if (cn != -1) {
      if (hex(cr).childs().size() == 0)
        continue;

      FaceVector f;
      HexLaO.childs_of_global_face(f, hex(cr), hp->first);

      for (unsigned i = 0; i < f.size(); ++i) {
        IndexType c = f[i];

        if ((CellRefList.find(c) != CellRefList.end()) || (hexs[c].sleep())) {
          if (CellCoarseList.find(cn) != CellCoarseList.end()) {
            coarsesub.push_back(cn);
            break;
          } else {
            assert(!hex(cn).sleep());

            pair<IntSetIt, bool> p = celllist.insert(cn);
            if (p.second) {
              n++;
              break;
            }
            // celllist.push_back(cn);
            // break;
          }
        }
      }
    }
  }
  return n + coarsesub.size();
  // unique(coarsesub.begin(),coarsesub.end());
  // unique(celllist.begin(),celllist.end());
  // return celllist.size()+coarsesub.size();
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh3d::regular_grid3d_one(IndexVector& celllist,
                                       IndexVector& coarsesub,
                                       const IndexSet& CellRefList,
                                       const IndexSet& CellCoarseList)
{
  IndexSet h;
  Vec2Set(h, celllist);
  IndexType r = regular_grid3d_one(h, coarsesub, CellRefList, CellCoarseList);
  Set2Vec(celllist, h);
  return r;

  /* detects jump over two levels across LineHangs */

  // IndexType n = 0;
  HangList<4>::const_iterator hp;

  for (hp = QuadHang.begin(); hp != QuadHang.end(); ++hp) {
    int cr = hp->second.rneighbour();
    int cn = hp->second.cneighbour();
    assert(cr >= 0);

    if (cn != -1) {
      if (hex(cr).childs().size() == 0)
        continue;

      FaceVector f;
      HexLaO.childs_of_global_face(f, hex(cr), hp->first);

      for (unsigned i = 0; i < f.size(); ++i) {
        IndexType c = f[i];

        if ((CellRefList.find(c) != CellRefList.end()) || (hexs[c].sleep())) {
          if (CellCoarseList.find(cn) != CellCoarseList.end()) {
            coarsesub.push_back(cn);
            break;
          } else {
            assert(!hex(cn).sleep());
            IndexType number = 0;
            for (IndexType kk = 0; kk < celllist.size(); kk++) {
              if (celllist[kk] == cn)
                number++;
            }
            if (!number) {
              celllist.push_back(cn);
              // n++;
              break;
            }
          }
        }
      }
    }
  }
  return celllist.size() + coarsesub.size();
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh3d::regular_grid3d_two(IndexVector& celllist,
                                       const IndexSet& CellRefList)
{
  /* detects more than 4 HangFaces on one Hex */

  IndexVector nh(hexs.size());
  for (HangList<4>::const_iterator p = QuadHang.begin(); p != QuadHang.end();
       p++) {
    IndexType i = p->second.cneighbour();
    if (i < 0)
      continue;
    if (hex(i).sleep())
      continue;
    if (CellRefList.find(i) != CellRefList.end())
      continue;

    nh[i]++;
  }
  for (IndexType i = 0; i < hexs.size(); i++) {
    if (nh[i] > 4) {
      // pair<IntSetIt,bool> pp = celllist.insert(i);
      // if(pp.second)  nto++;
      celllist.push_back(i);
      // nto++;
    }
  }
  return celllist.size();
  // return nto;
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::GetMinMaxLevels(IndexVector& maxi,
                                    IndexVector& mini,
                                    const IndexSet& CellRef) const
{
  // set maximal levels for vertices
  //
  maxi.resize(nnodes());
  mini.resize(nnodes());
  mini = 1000;
  for (IndexType i = 0; i < hexs.size(); i++) {
    const Hex& q = hex(i);
    if (q.sleep())
      continue;

    IndexType lev = q.level();
    if (CellRef.find(i) != CellRef.end())
      lev++;
    for (IndexType j = 0; j < 8; j++) {
      IndexType k = q[j];
      maxi[k] = std::max(maxi[k], lev);
      mini[k] = std::min(mini[k], lev);
    }
  }
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh3d::regular_grid3d_three_refine(IndexSet& CellRef) const
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
  for (IndexType i = 0; i < hexs.size(); i++) {
    const Hex& q = hex(i);
    if (q.sleep())
      continue;
    IndexType lev = q.level();
    for (IndexType j = 0; j < 8; j++) {
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
  return ref;
}

/*---------------------------------------------------*/

IndexType
HierarchicalMesh3d::regular_grid3d_three_coarse(IndexSet& CellRef,
                                                IndexSet& CellCoarse) const
{
  IndexType maxl = 0;
  vector<IndexSet> LevelCellCoarse;
  {
    IndexSet::const_iterator p = CellCoarse.begin();
    for (; p != CellCoarse.end(); p++) {
      const Hex& q = hex(*p);
      IndexType lev = q.level();
      maxl = std::max(maxl, lev);
    }
    LevelCellCoarse.resize(maxl + 1);
    p = CellCoarse.begin();
    for (; p != CellCoarse.end(); p++) {
      const Hex& q = hex(*p);
      IndexType lev = q.level();
      LevelCellCoarse[lev].insert(*p);
    }
  }
  IndexType coarse = 0;
  for (IndexType i = maxl; i >= 0; i--) {
    IndexVector maxlevel, minlevel;

    GetMinMaxLevels(maxlevel, minlevel, CellRef);

    IndexSet coarsesub;
    IndexSet::const_iterator p = LevelCellCoarse[i].begin();
    for (; p != LevelCellCoarse[i].end(); p++) {
      const Hex& q = hex(*p);
      IndexType lev = q.level();
      for (IndexType j = 0; j < 8; j++) {
        IndexType k = q[j];
        if (lev + 1 < maxlevel[k]) {
          coarsesub.insert(*p);
          continue;
        }
      }
    }
    p = coarsesub.begin();
    for (; p != coarsesub.end(); p++) {
      if (CellCoarse.find(*p) != CellCoarse.end()) {
        CellCoarse.erase(*p);
        coarse++;
      }
    }
  }
  return coarse;
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::global_coarse3d()
{
  IndexSet RefList, CoarseList;

  for (IndexType i = 0; i < hexs.size(); i++) {
    const Hex& q = hexs[i];
    if (q.sleep())
      continue;
    if (!q.level())
      continue;
    CoarseList.insert(q.father());
  }

  HangContainer3d hangset(LineHang, QuadHang);

  UpdateHangs(hangset, RefList, CoarseList);

  basic_refine3d(hangset, RefList, CoarseList);

  post_refine3d();
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::LoadFathers3d(IndexVector& v) const
{
  IndexSet fathers;

  for (IndexType i = 0; i < v.size(); i++) {
    if (hex(v[i]).level() == 0)
      continue;
    IndexType f = hex(v[i]).father();

    assert(f >= 0);

    fathers.insert(f);
  }
  Set2Vec(v, fathers);
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::FillVertexLevels(IndexVector& dst) const
{
  dst.resize(nnodes(), 100);
  for (IndexType i = 0; i < ncells(); i++) {
    const Hex& Q = hexs[i];
    if (Q.sleep())
      continue;

    IndexType level = Q.level();
    for (IndexType j = 0; j < 8; j++) {
      IndexType k = Q[j];
      dst[k] = std::min(dst[k], level);
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::RefineCoarseNodes(IndexSet& dst,
                                      const IndexVector& refnodes,
                                      const IndexVector& vertexlevel) const
{
  IndexSet h;

  Vec2Set(h, refnodes);

  dst.clear();
  for (IndexType i = 0; i < ncells(); i++) {
    const Hex& Q = hexs[i];
    if (Q.sleep())
      continue;

    IndexType f = Q.father();
    if (f < 0)
      continue;

    const Hex& QF = hexs[f];

    for (IndexType j = 0; j < 8; j++) {
      IndexType k = Q[j];
      if (h.find(k) == h.end())
        continue;

      IndexType minlevel = vertexlevel[QF[0]];
      for (IndexType v = 1; v < 8; v++) {
        minlevel = std::min(minlevel, vertexlevel[QF[v]]);
      }
      for (IndexType v = 0; v < 8; v++) {
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
HierarchicalMesh3d::VertexToCells(IndexVector& dst,
                                  const IndexSet& src,
                                  const IndexVector& vertexlevel) const
{
  for (IndexType i = 0; i < ncells(); i++) {
    const Hex& Q = hexs[i];
    if (Q.sleep())
      continue;
    IndexType level = Q.level();
    for (IndexType j = 0; j < 8; j++) {
      IndexType k = Q[j];
      if (vertexlevel[k] == level) {
        if (src.find(k) != src.end()) {
          dst.push_back(i);
        }
      }
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::VertexToCellsCoarsening(
  IndexVector& dst,
  const IndexSet& src,
  const IndexVector& vertexlevel) const
{
  for (IndexType i = 0; i < ncells(); i++) {
    const Hex& Q = hexs[i];
    if (Q.sleep())
      continue;
    IndexType level = Q.level();
    IndexType count = 0;
    for (IndexType j = 0; j < 8; j++) {
      IndexType k = Q[j];
      if (vertexlevel[k] == level) {
        if (src.find(k) != src.end()) {
          count++;
        } else
          break;
      }
    }
    if (count == 8)
      dst.push_back(i);
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::recursive_childs(IndexType q,
                                     IndexVector& ref,
                                     IndexType d) const
{
  if (d > 0) {
    const Hex& Q = hex(q);
    assert(Q.sleep());
    for (IndexType j = 0; j < 8; j++) {
      IndexType child = Q.child(j);
      recursive_childs(child, ref, d - 1);
    }
  } else {
    ref.push_back(q);
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::patch_refine(IndexVector& cell_ref,
                                 IndexVector& cell_coarse)
{
  for (IndexType i = 0; i < pdepth; i++) {
    LoadFathers3d(cell_ref);
    LoadFathers3d(cell_coarse);
  }

  CoarseHierarchicalMesh3d CM(*this);

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
HierarchicalMesh3d::refine(const IndexVector& cell_ref,
                           const IndexVector& cell_coarse)
{
  IndexSet CellRefList, CellCoarseList;
  _refine3d(CellRefList, CellCoarseList, cell_ref, cell_coarse);
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::_refine3d(IndexSet& CellRefList,
                              IndexSet& CellCoarseList,
                              const IndexVector& cell_ref,
                              const IndexVector& cell_coarse)
{
  prepare3d(cell_ref, cell_coarse, CellRefList, CellCoarseList);

  while (regular_grid3d_three_refine(CellRefList)) {
  }
  while (regular_grid3d_three_coarse(CellRefList, CellCoarseList)) {
  }

  HangContainer3d hangset(LineHang, QuadHang);

  UpdateHangs(hangset, CellRefList, CellCoarseList);

  basic_refine3d(hangset, CellRefList, CellCoarseList);

  post_refine3d();
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::delete_vertexs3d(const IndexVector& vo2n)
{
  for (unsigned oi = 0; oi < vo2n.size(); ++oi) {
    IndexType ni = vo2n[oi];
    if (ni >= 0) {
      vertexs3d[ni] = vertexs3d[oi];
    }
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::post_refine3d()
{
  // check_mesh3d();
  mnlevels = 0;
  for (IndexType i = 0; i < hexs.size(); i++) {
    mnlevels = std::max(mnlevels, hexs[i].level());
  }
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::new_edge_vertex3d(IndexType nv, const EdgeVector& v)
{
  vertexs3d[nv].equ(0.5, vertexs3d[v[0]], 0.5, vertexs3d[v[1]]);
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::new_face_vertex3d(IndexType nv, const FaceVector& v)
{
  vertexs3d[nv].equ(0.25,
                    vertexs3d[v[0]],
                    0.25,
                    vertexs3d[v[1]],
                    0.25,
                    vertexs3d[v[2]],
                    0.25,
                    vertexs3d[v[3]]);
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::new_vertex3d(IndexType nv,
                                 const std::array<IndexType, 6>& v)
{
  cout << "@ " << vertexs3d[nv] << "\t";
  double d = 1. / 6.;
  vertexs3d[nv].equ(d,
                    vertexs3d[v[0]],
                    d,
                    vertexs3d[v[1]],
                    d,
                    vertexs3d[v[2]],
                    d,
                    vertexs3d[v[3]],
                    d,
                    vertexs3d[v[4]],
                    d,
                    vertexs3d[v[5]]);
  cout << " -> " << vertexs3d[nv] << "\n";
}

/*---------------------------------------------------*/

void
HierarchicalMesh3d::check_mesh3d() const
{

  // check quads

  IndexType cmin = 0, cmax = hexs.size();
  IndexType vmin = 0, vmax = nnodes();

  for (vector<Hex>::const_iterator p = hexs.begin(); p != hexs.end(); p++) {
    const Hex& q = *p;

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

  // check quadhang
  for (HangList<4>::const_iterator hp = QuadHang.begin(); hp != QuadHang.end();
       ++hp) {
    IndexType cr = hp->second.rneighbour();
    if (cr == -1) {
      cerr << "Refine Neighbour invalid in hang: ";
      // cerr << *hp << endl;
      exit(1);
    }
  }
}
} // namespace Gascoigne

/*---------------------------------------------------*/
