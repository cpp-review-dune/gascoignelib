/**
 *
 * Copyright (C) 2004, 2006, 2007, 2008, 2011 by the Gascoigne 3D authors
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

#include "levelmesh3d.h"
#include "gascoignehash.h"
#include "leveljumper.h"
#include "levelsorter3d.h"
#include "nmatrix.h"
#include "set2vec.h"
#include <algorithm>

using namespace std;

/*---------------------------------------------------*/

namespace Gascoigne {

LevelMesh3d::LevelMesh3d(const HierarchicalMesh* hmp)
  : Index()
{
  HMP = dynamic_cast<const HierarchicalMesh3d*>(hmp);
}

/*---------------------------------------------------*/

LevelMesh3d::~LevelMesh3d() {}

/*---------------------------------------------------*/

void
LevelMesh3d::BasicInit(const IndexSet& newh, const IndexSet& oldh)
{
  IndexType n = newh.size() + oldh.size();
  Index::Hexl2g().memory(n);

  IndexVector::const_iterator p = set_union(newh.begin(),
                                            newh.end(),
                                            oldh.begin(),
                                            oldh.end(),
                                            Index::Hexl2g().begin());
  n = p - Index::Hexl2g().begin();

  InitCells(n);
  InitNodes(n);
  InitEdges(n);
}

/*-----------------------------------------*/

void
LevelMesh3d::InitCells(IndexType n)
{
  Index::Hexl2g().memory(n);

  sort(Index::Hexl2g().begin(), Index::Hexl2g().end(), LevelSorter3d(*HMP));

  Index::InitHexs();
}

/*-----------------------------------------*/

void
LevelMesh3d::InitNodes(IndexType n)
{
  IndexSet nodes;

  for (IndexType i = 0; i < n; i++) {
    IndexType ind = Index::Hexl2g()[i];
    for (IndexType ii = 0; ii < 8; ii++) {
      nodes.insert(HMP->vertex_of_cell(ind, ii));
    }
  }
  Index::InitNodes(nodes);
}

/*-----------------------------------------*/

void
LevelMesh3d::InitEdges(IndexType n)
{
  IndexSet edges;
  for (IndexType i = 0; i < n; i++) {
    const Hex& Q = HMP->hex(Index::Hexl2g()[i]);
    for (IndexType ii = 0; ii < 6; ii++) {
      edges.insert(Q.edge(ii));
    }
  }
  Index::InitEdges(edges);

  // stable_sort(Edgel2g().begin(),Edgel2g().end(),HangEdgeSort5(*this));

  Edgeg2l().clear();
  for (IndexType i = 0; i < Edgel2g().size(); i++) {
    Edgeg2l().insert(make_pair(Edgel2g()[i], i));
  }
}

/*-----------------------------------------*/

bool
LevelMesh3d::BuildFathers(IndexSet& Vaeter) const
{
  for (IndexType i = 0; i < ncells(); i++) {
    const Hex& h = hex(i);
    IndexType findex = h.father();
    if (findex == -1) {
      return 0;
    }

    const Hex& hf = HMP->hex(findex);
    for (IndexType ii = 0; ii < hf.nchilds(); ii++) {
      IndexType cindex = hf.child(ii);
      if (Hexg2lCheck(cindex) == -2) {
        return 0;
      }
    }
    Vaeter.insert(findex);
  }
  return 1;
}

/*---------------------------------------------------*/

void
LevelMesh3d::ConstructIndOfPatch(nvector<IndexVector>& dst) const
{
  set<IndexType> Vaeter;
  BuildFathers(Vaeter);

  IndexType nh = ncells() / 8;
  dst.reserve(nh);
  dst.resize(nh, IndexVector(27));

  IndexType count = 0;
  set<IndexType>::const_iterator pf = Vaeter.begin();
  nmatrix<IndexType> A(27, 2);

  A(0, 0) = 0;
  A(0, 1) = 0;
  A(1, 0) = 0;
  A(1, 1) = 1;
  A(2, 0) = 1;
  A(2, 1) = 1;
  A(3, 0) = 0;
  A(3, 1) = 3;
  A(4, 0) = 0;
  A(4, 1) = 2;
  A(5, 0) = 1;
  A(5, 1) = 2;
  A(6, 0) = 3;
  A(6, 1) = 3;
  A(7, 0) = 3;
  A(7, 1) = 2;
  A(8, 0) = 2;
  A(8, 1) = 2;

  A(9, 0) = 4;
  A(9, 1) = 0;
  A(10, 0) = 4;
  A(10, 1) = 1;
  A(11, 0) = 5;
  A(11, 1) = 1;
  A(12, 0) = 4;
  A(12, 1) = 3;
  A(13, 0) = 4;
  A(13, 1) = 2;
  A(14, 0) = 5;
  A(14, 1) = 2;
  A(15, 0) = 7;
  A(15, 1) = 3;
  A(16, 0) = 7;
  A(16, 1) = 2;
  A(17, 0) = 6;
  A(17, 1) = 2;

  A(18, 0) = 4;
  A(18, 1) = 4;
  A(19, 0) = 4;
  A(19, 1) = 5;
  A(20, 0) = 5;
  A(20, 1) = 5;
  A(21, 0) = 4;
  A(21, 1) = 7;
  A(22, 0) = 4;
  A(22, 1) = 6;
  A(23, 0) = 5;
  A(23, 1) = 6;
  A(24, 0) = 7;
  A(24, 1) = 7;
  A(25, 0) = 7;
  A(25, 1) = 6;
  A(26, 0) = 6;
  A(26, 1) = 6;

  while (pf != Vaeter.end()) {
    IndexType findex = *pf;
    const Hex& hf = HMP->hex(findex);
    std::array<IndexType, 8> FineHexs;
    for (IndexType ii = 0; ii < hf.nchilds(); ii++) {
      FineHexs[ii] = hf.child(ii);
    }
    for (IndexType i = 0; i < 27; i++) {
      IndexType fh = FineHexs[A(i, 0)];
      IndexType gi = Vertexg2l(HMP->hex(fh).vertex(A(i, 1)));
      dst[count][i] = gi;
    }
    count++;
    pf++;
  }
}

/*---------------------------------------------------*/

bool
LevelMesh3d::ConstructCellIndOfPatch(IndexVector& dst) const
{
  set<IndexType> Vaeter;
  BuildFathers(Vaeter);

  IndexType nh = ncells() / 8;
  dst.reservesize(nh);

  IndexType count = 0;
  set<IndexType>::const_iterator pf = Vaeter.begin();

  while (pf != Vaeter.end()) {
    IndexType findex = *pf;
    dst[count] = findex;

    count++;
    pf++;
  }
  return 1;
}

/*---------------------------------------------------*/

void
LevelMesh3d::ConstructHangingStructureQuadratic(
  QuadraticHNStructure3& hnq2,
  QuadraticHNStructure9& hnq2face) const
{
  hnq2.clear();
  hnq2face.clear();
  set<IndexType> habschon;

  const HexLawAndOrder& HLaO = HMP->HexLawOrder();

  for (IndexType i = 0; i < ncells(); i++) {
    const Hex& q = HMP->hex(Hexl2g(i));
    IndexType father = q.father();
    if (father < 0)
      continue;
    if (habschon.find(father) != habschon.end())
      continue;
    habschon.insert(father);
    for (IndexType in = 0; in < 6; in++) {
      IndexType neighbour = HMP->neighbour(father, in);
      if (neighbour < 0)
        continue;
      const Hex& qfn = HMP->hex(neighbour);
      if (qfn.nchilds() == 0)
        continue;

      IndexType ne = in;
      {
        IndexType start = 0;
        IndexType neighbourneighbour = HMP->neighbour(neighbour, ne);
        while ((neighbourneighbour != father) && (start < 6)) {
          start++;
          ne = (ne + 1) % 6;
          neighbourneighbour = HMP->neighbour(neighbour, ne);
        }
        assert(neighbourneighbour == father);
      }

      std::array<IndexType, 4> childs;
      HLaO.childs_of_face(childs, qfn, ne);

      {
        IndexType child = childs[0];
        if (Hexg2lCheck(child) >= 0)
          continue;

        const Hex& qfc = HMP->hex(child);

        if (qfc.nchilds() == 0)
          continue;

        if (Hexg2lCheck(qfc.child(0)) < 0)
          continue;
      }

      // jetzt haengt er

      std::array<IndexType, 9> F = HLaO.PatchVerticesOfFace(neighbour, ne);
      IndexType nec = HLaO.ChildFace(ne);

      // ordne F;
      for (IndexType j = 0; j < 4; j++) {
        const Hex& qfcc = HMP->hex(childs[j]);

        IndexType hnl = Vertexg2l(HLaO.face_vertex(qfcc, nec));

        // permutiere F
        std::array<IndexType, 9> G = HLaO.GiveOrdering(F, qfcc);

        std::array<IndexType, 4> face;
        face[0] = F[G[0]];
        face[1] = F[G[1]];
        face[2] = F[G[3]];
        face[3] = F[G[4]];

        for (IndexType k = 0; k < 9; k++)
          G[k] = Vertexg2l(F[G[k]]);

        hnq2face.insert(make_pair(hnl, G));

        std::array<IndexType, 4> childface;
        HLaO.GetFace(childface, childs[j], nec);
        // nun die hanging lines
        for (IndexType e = 0; e < 4; e++) {
          IndexType hne = Vertexg2l(HLaO.EdgeVertexOfFace(qfcc, face, e));
          if (hnq2.find(hne) == hnq2.end()) {
            std::array<IndexType, 3> line;
            IndexType e0 = childface[e];
            IndexType e1 = childface[(e + 1) % 4];

            if (e0 == F[4])
              swap(e1, e0);
            if (e1 != F[4]) {
              if ((e0 == F[1]) || (e0 == F[3]) || (e0 == F[5]) ||
                  (e0 == F[7])) {
                swap(e1, e0);
              }
            }

            line[0] = Vertexg2l(e0);
            line[1] = Vertexg2l(e1);

            IndexType last;
            if ((e0 == F[0]) && (e1 == F[1]))
              last = F[2];
            else if ((e0 == F[2]) && (e1 == F[1]))
              last = F[0];
            else if ((e0 == F[3]) && (e1 == F[4]))
              last = F[5];
            else if ((e0 == F[5]) && (e1 == F[4]))
              last = F[3];
            else if ((e0 == F[6]) && (e1 == F[7]))
              last = F[8];
            else if ((e0 == F[8]) && (e1 == F[7]))
              last = F[6];
            else if ((e0 == F[0]) && (e1 == F[3]))
              last = F[6];
            else if ((e0 == F[6]) && (e1 == F[3]))
              last = F[0];
            else if ((e0 == F[1]) && (e1 == F[4]))
              last = F[7];
            else if ((e0 == F[7]) && (e1 == F[4]))
              last = F[1];
            else if ((e0 == F[2]) && (e1 == F[5]))
              last = F[8];
            else if ((e0 == F[8]) && (e1 == F[5]))
              last = F[2];
            else
              abort();

            line[2] = Vertexg2l(last);
            hnq2.insert(make_pair(hne, line));
          }
        }
      }
    }
  }
  //  cout << "****** LevelM 3d " << hnq2face.size() << " " << hnq2 .size() <<
  //  endl;
}

/*---------------------------------------------------*/

void
LevelMesh3d::ConstructHangingStructureQuartic(
  QuarticHNStructure5& hnq4,
  QuarticHNStructure25& hnq4face) const
{
  hnq4.clear();
  hnq4face.clear();
  assert(HMP->patchdepth() >= 2);

  std::set<IndexType> habschon;
  std::array<IndexType, 81> nodesonface;

  for (IndexType hexl = 0; hexl < ncells(); ++hexl) {
    IndexType hex = Hexl2g(hexl);
    // Opa holen
    IndexType vater = HMP->hex(hex).father();
    assert(vater >= 0);
    IndexType opa = HMP->hex(vater).father();
    assert(opa >= 0);
    // Jeden Opa nur einmal betrachten
    if (habschon.find(opa) != habschon.end())
      continue;
    habschon.insert(opa);
    // alle nachbarn durchlaufen
    for (IndexType ni = 0; ni < 6; ++ni) {
      IndexType neighbour = HMP->neighbour(opa, ni);
      if (neighbour < 0)
        continue;
      // Face haengt, wenn Nachbar einmal mehr, also 3-mal verfeinert ist.
      // Wenn Nachbar groeber ist, dann haengt die Face eben von der anderen
      // Seite
      if (refine_level(neighbour) < 3)
        continue;
      // Liste von allen 81 Knoten auf dieser Face aufbauen
      ConstructNodesOnFaceQ4(nodesonface, opa, ni);
      // und die haengende face einfuegen
      InsertHangingFacesQ4(hnq4face, nodesonface);
      // und gegebenenfalls die edges
      InsertHangingEdgesQ4(hnq4, nodesonface);
    }
  }
}

/*---------------------------------------------------*/

// tiefe, bis Kinder aktiv sind
IndexType
LevelMesh3d::refine_level(IndexType n) const
{
  IndexType l = 0;
  while (Hexg2lCheck(n) < 0) {
    ++l;
    const Hex& h = HMP->hex(n);
    assert(h.nchilds() > 0);
    n = h.child(0);
  }
  return l;
}

/*---------------------------------------------------*/

// Sucht alle 81 Knoten, die an einer haengenden
// Flaeche einer groben Q4-Zelle liegen
// Dazu wird die entsprechende Funktion der Q2-Face verwendet
void
LevelMesh3d::ConstructNodesOnFaceQ4(std::array<IndexType, 81>& nodesonface,
                                    IndexType opa,
                                    IndexType ni) const
{
  const HexLawAndOrder& HLaO = HMP->HexLawOrder();
  assert(opa >= 0);
  IndexType neighbour = HMP->neighbour(opa, ni);
  assert(neighbour >= 0);
  // nachbar-nachbar index suchen
  IndexType nni = HMP->neighbour_neighbour(opa, ni);

  // Die vier Kinder der groben Zelle an dieser Face
  std::array<IndexType, 4> childs;
  HLaO.childs_of_face(childs, HMP->hex(opa), ni);

  // die 25 Patchknoten auf einer Kindesface
  std::vector<std::array<IndexType, 25>> F(4);
  for (IndexType i = 0; i < 4; ++i)
    ConstructNodesOnFace(F[i], childs[i], ni);
  // der Mittelknoten der Vaters
  IndexType middle_node = HLaO.face_vertex(HMP->hex(neighbour), nni);
  // Die zellen selbst koennen  rotiert sein. Alle 4
  // haben aber die gleiche Orientierung
  // die zellen so sortieren, dass sie von links unten
  // nach recht oben lexikogr. liegen
  // Dazu wird jeweils die Zelle gesucht, so dass der Mittelknoten richtig liegt
  IndexType edges[4] = { 24, 20, 4, 0 };
  std::array<IndexType, 4> s;
  for (IndexType c = 0; c < 4;
       ++c) { // zelle finden, die an Position c liegt, also
              // zelle, von der edges[c] die mittelzelle ist.
    for (s[c] = 0; s[c] < 4; ++s[c])
      if (F[s[c]][edges[c]] == middle_node)
        break;
  }

  // Index zusammenbasteln
  for (IndexType y = 0; y < 5; ++y)
    for (IndexType x = 0; x < 5; ++x)
      nodesonface[y * 9 + x] = F[s[0]][y * 5 + x];
  for (IndexType y = 0; y < 5; ++y)
    for (IndexType x = 5; x < 9; ++x)
      nodesonface[y * 9 + x] = F[s[1]][y * 5 + x - 4];
  for (IndexType y = 5; y < 9; ++y)
    for (IndexType x = 0; x < 5; ++x)
      nodesonface[y * 9 + x] = F[s[2]][(y - 4) * 5 + x];
  for (IndexType y = 5; y < 9; ++y)
    for (IndexType x = 5; x < 9; ++x)
      nodesonface[y * 9 + x] = F[s[3]][(y - 4) * 5 + x - 4];
}

/*---------------------------------------------------*/

void
LevelMesh3d::InsertHangingFacesQ4(QuarticHNStructure25& hnq4face,
                                  const std::array<IndexType, 81>& nof) const
{
  std::array<IndexType, 25> I;

  for (IndexType y = 0; y < 5; ++y)
    for (IndexType x = 0; x < 5; ++x)
      I[y * 5 + x] = nof[y * 18 + x * 2]; // links unten
  InsertHangingFaceQ4(hnq4face, nof, 10, 12, 28, 30, I);
  for (IndexType y = 0; y < 5; ++y)
    for (IndexType x = 0; x < 5; ++x)
      I[y * 5 + x] = nof[x * 18 + (4 - y) * 2]; // rechts unten
  InsertHangingFaceQ4(hnq4face, nof, 16, 34, 14, 32, I);
  for (IndexType y = 0; y < 5; ++y)
    for (IndexType x = 0; x < 5; ++x)
      I[y * 5 + x] = nof[(4 - y) * 18 + (4 - x) * 2]; // rechts oben
  InsertHangingFaceQ4(hnq4face, nof, 70, 68, 52, 50, I);
  for (IndexType y = 0; y < 5; ++y)
    for (IndexType x = 0; x < 5; ++x)
      I[y * 5 + x] = nof[(4 - x) * 18 + y * 2]; // links oben
  InsertHangingFaceQ4(hnq4face, nof, 64, 46, 66, 48, I);
}

/*---------------------------------------------------*/

void
LevelMesh3d::InsertHangingEdgesQ4(QuarticHNStructure5& hnq4,
                                  const std::array<IndexType, 81>& nof) const
{
  // horizontal
  for (IndexType y = 0; y < 5; ++y)
    InsertHangingEdgeQ4(hnq4,
                        nof,
                        18 * y + 1,
                        18 * y + 3,
                        18 * y + 5,
                        18 * y + 7,
                        18 * y + 0,
                        18 * y + 2,
                        18 * y + 4,
                        18 * y + 6,
                        18 * y + 8);
  // vertikal
  for (IndexType x = 0; x < 5; ++x)
    InsertHangingEdgeQ4(hnq4,
                        nof,
                        2 * x + 9,
                        2 * x + 27,
                        2 * x + 45,
                        2 * x + 63,
                        2 * x + 0,
                        2 * x + 18,
                        2 * x + 36,
                        2 * x + 54,
                        2 * x + 72);
}

/*---------------------------------------------------*/

// Sucht alle 25 Knoten, die an einer haengenden
// Flaeche einer groben Q2-Zelle liegen
void
LevelMesh3d::ConstructNodesOnFace(std::array<IndexType, 25>& nodesonface,
                                  IndexType vater,
                                  IndexType ni) const
{
  const HexLawAndOrder& HLaO = HMP->HexLawOrder();
  assert(vater >= 0);
  IndexType neighbour = HMP->neighbour(vater, ni);
  assert(neighbour >= 0);
  // nachbar-nachbar index suchen
  IndexType nni = HMP->neighbour_neighbour(vater, ni);

  // Die vier Kinder an dieser Face
  // Diese vier Kinder sind gegen den Uhrzeigersinn angeordnet.
  std::array<IndexType, 4> childs;
  HLaO.childs_of_face(childs, HMP->hex(neighbour), nni);

  // die 9 Patchknoten auf einer Kindesface
  std::vector<std::array<IndexType, 9>> F(4);
  for (IndexType i = 0; i < 4; ++i)
    F[i] = HLaO.PatchVerticesOfFace(childs[i], nni);
  // der Mittelknoten der Vaters
  IndexType middle_node = HLaO.face_vertex(HMP->hex(neighbour), nni);
  // Die zellen selbst koennen  rotiert sein. Alle 4
  // haben aber die gleiche Orientierung
  // die zellen so sortieren, dass sie von links unten
  // nach recht oben lexikogr. liegen
  // Dazu wird jeweils die Zelle gesucht, so dass der Mittelknoten richtig liegt
  IndexType edges[4] = { 8, 6, 2, 0 };
  std::array<IndexType, 4> s;
  for (IndexType c = 0; c < 4;
       ++c) { // zelle finden, die an Position c liegt, also
              // zelle, von der edges[c] die mittelzelle ist.
    for (s[c] = 0; s[c] < 4; ++s[c])
      if (F[s[c]][edges[c]] == middle_node)
        break;
  }
  // Index zusammenbasteln
  for (IndexType y = 0; y < 3; ++y)
    for (IndexType x = 0; x < 3; ++x)
      nodesonface[y * 5 + x] = F[s[0]][y * 3 + x];
  for (IndexType y = 0; y < 3; ++y)
    for (IndexType x = 3; x < 5; ++x)
      nodesonface[y * 5 + x] = F[s[1]][y * 3 + x - 2];
  for (IndexType y = 3; y < 5; ++y)
    for (IndexType x = 0; x < 3; ++x)
      nodesonface[y * 5 + x] = F[s[2]][(y - 2) * 3 + x];
  for (IndexType y = 3; y < 5; ++y)
    for (IndexType x = 3; x < 5; ++x)
      nodesonface[y * 5 + x] = F[s[3]][(y - 2) * 3 + x - 2];
}

/*---------------------------------------------------*/

void
LevelMesh3d::InsertHangingFaceQ4(QuarticHNStructure25& hnq4face,
                                 const std::array<IndexType, 81>& nof,
                                 IndexType n1,
                                 IndexType n2,
                                 IndexType n3,
                                 IndexType n4,
                                 const std::array<IndexType, 25>& I) const
{
  // eine haengende face hat 4 haengende. Die vier typen dieser Knoten
  // sind so wie die 9 Stuetzknoten lexikografisch vergeben.
  assert(hnq4face.find(Vertexg2l(nof[n1])) == hnq4face.end());
  std::array<IndexType, 26> G;
  for (IndexType i = 0; i < 25; ++i)
    G[i] = Vertexg2l(I[i]);
  G[25] = 0;
  hnq4face.insert(std::make_pair(Vertexg2l(nof[n1]), G));
  G[25] = 1;
  hnq4face.insert(std::make_pair(Vertexg2l(nof[n2]), G));
  G[25] = 2;
  hnq4face.insert(std::make_pair(Vertexg2l(nof[n3]), G));
  G[25] = 3;
  hnq4face.insert(std::make_pair(Vertexg2l(nof[n4]), G));
}

/*---------------------------------------------------*/

void
LevelMesh3d::InsertHangingEdgeQ4(QuarticHNStructure5& hnq4,
                                 const std::array<IndexType, 81>& nof,
                                 IndexType n1,
                                 IndexType n2,
                                 IndexType n3,
                                 IndexType n4,
                                 IndexType i1,
                                 IndexType i2,
                                 IndexType i3,
                                 IndexType i4,
                                 IndexType i5) const
{
  // die ganze Kante muss neu haengen. Es gibt 2 mal 2 haengende Knoten.
  // Knoten am Rand haben jeweils den Typ 0, die in der Mitte Typ 1
  if (hnq4.find(Vertexg2l(nof[n1])) != hnq4.end()) {
    assert(hnq4.find(Vertexg2l(nof[n2])) != hnq4.end());
    assert(hnq4.find(Vertexg2l(nof[n3])) != hnq4.end());
    assert(hnq4.find(Vertexg2l(nof[n4])) != hnq4.end());
    return;
  }
  std::array<IndexType, 6> G;
  G[0] = Vertexg2l(nof[i1]);
  G[1] = Vertexg2l(nof[i2]);
  G[2] = Vertexg2l(nof[i3]);
  G[3] = Vertexg2l(nof[i4]);
  G[4] = Vertexg2l(nof[i5]);
  G[5] = 0;
  hnq4.insert(std::make_pair(Vertexg2l(nof[n1]), G));
  G[5] = 1;
  hnq4.insert(std::make_pair(Vertexg2l(nof[n2]), G));
  std::swap(G[0], G[4]);
  std::swap(G[1], G[3]);
  G[5] = 1;
  hnq4.insert(std::make_pair(Vertexg2l(nof[n3]), G));
  G[5] = 0;
  hnq4.insert(std::make_pair(Vertexg2l(nof[n4]), G));
}

/*---------------------------------------------------*/

bool
LevelMesh3d::EnkelUniform(const Hex& Q) const
{
  for (IndexType ii = 0; ii < Q.nchilds(); ii++) {
    IndexType qccindex = Q.child(ii);
    const Hex& qcc = HMP->hex(qccindex);
    for (IndexType iii = 0; iii < qcc.nchilds(); iii++) {
      IndexType qcindex = qcc.child(iii);
      if (Hexg2lCheck(qcindex) == -2) {
        return 0;
      }
    }
  }
  return 1;
}

/*---------------------------------------------------*/

void
LevelMesh3d::fill_opis(IndexSet& dst, IndexSet& oldhexs) const
{
  for (IndexType i = 0; i < ncells(); i++) {
    const Hex& Q = hex(i);

    IndexType f = Q.father();
    assert(f >= 0);

    IndexType opi = HMP->hex(f).father();

    if (opi < 0) {
      IndexType j = Hexl2g(i);
      oldhexs.insert(j);
    } else {
      dst.insert(opi);
    }
  }
}

/*---------------------------------------------------*/

void
LevelMesh3d::fill_childs(IndexSet& dst, const Hex& Q) const
{
  for (IndexType i = 0; i < Q.nchilds(); i++) {
    dst.insert(Q.child(i));
  }
}

/*---------------------------------------------------*/

void
LevelMesh3d::fill_enkel(IndexSet& dst, const Hex& Q) const
{
  for (IndexType i = 0; i < Q.nchilds(); i++) {
    const Hex& C = HMP->hex(Q.child(i));
    for (IndexType j = 0; j < C.nchilds(); j++) {
      IndexType cc = C.child(j);
      if (Hexg2lCheck(cc) >= 0) {
        dst.insert(cc);
      }
    }
  }
}

/*---------------------------------------------------*/

void
LevelMesh3d::check_leveljump() const
{
  LevelJumper Phi;
  for (IndexType c = 0; c < ncells(); c++) {
    Phi.update(hex(c));
  }
  assert(!Phi.check());
}

/*---------------------------------------------------*/

void
LevelMesh3d::construct_lists(IndexSet& newhexs, IndexSet& oldhexs) const
{
  newhexs.clear();
  oldhexs.clear();

  check_leveljump();

  set<IndexType> Opis;
  fill_opis(Opis, oldhexs);
  for (set<IndexType>::const_iterator p = Opis.begin(); p != Opis.end(); p++) {
    const Hex& Q = HMP->hex(*p);

    if (EnkelUniform(Q)) {
      fill_childs(newhexs, Q);
    } else {
      fill_enkel(oldhexs, Q);
    }
  }

  // Iteration zum Regulaer machen (LevelJump)
  while (1) {
    LevelJumper Phi;
    for (set<IndexType>::const_iterator p = newhexs.begin(); p != newhexs.end();
         p++) {
      Phi.update(HMP->hex(*p));
    }
    for (set<IndexType>::const_iterator p = oldhexs.begin(); p != oldhexs.end();
         p++) {
      Phi.update(HMP->hex(*p));
    }
    if (!Phi.check())
      break;

    IndexType rep = 0;
    IndexSet help(newhexs);
    for (set<IndexType>::const_iterator p = newhexs.begin(); p != newhexs.end();
         p++) {
      const Hex& q = HMP->hex(*p);
      if (!Phi.VertexOK(q)) {
        rep++;
        const Hex& qf = HMP->hex(q.father());
        for (IndexType ii = 0; ii < 8; ii++) {
          help.erase(qf.child(ii));
        }
        fill_enkel(oldhexs, qf);
      }
    }
    newhexs = help;
    //       cerr << "\t Regular Iteration\t" << count++ << " " << rep << endl;
  }
}

/*---------------------------------------------------*/

void
LevelMesh3d::InitBoundaryHandler(BoundaryIndexHandler& BI,
                                 const PatchIndexHandler& PIH) const
{
  // bquads
  // probably slowly (could be in multigridmesh !)
  // since we cannnot go from quad -> bline
  IndexSet bquads;
  for (IndexType i = 0; i < HMP->nbquads(); i++) {
    IndexType q = HMP->bquad(i).of_quad();
    if (Hexg2l().find(q) != Hexg2l().end()) {
      bquads.insert(i);
    }
  }

  // which colors are there ?
  BI.clear();
  for (IndexSet::const_iterator p = bquads.begin(); p != bquads.end(); p++) {
    const BoundaryQuad& bl = HMP->bquad(*p);
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
  vector<set<IndexType>> H1(nc); // for verteces
  // for cells and local indices
  vector<set<std::array<IndexType, 2>>> H2(nc);

  for (IndexSet::const_iterator q = bquads.begin(); q != bquads.end(); q++) {
    const BoundaryQuad& bl = HMP->bquad(*q);
    IndexType col = bl.material();

    map<IndexType, IndexType>::const_iterator p = inv.find(col);
    if (p == inv.end()) {
      cout << "BoundaryIndexHandler::init3d" << endl;
      abort();
    }
    IndexType pos = p->second;

    for (IndexType ii = 0; ii < 4; ii++) {
      IndexType vindex = Vertexg2l(bl.vertex(ii));
      H1[pos].insert(vindex);
    }
    std::array<IndexType, 2> ind;
    ind[0] = Hexg2l(bl.of_quad());
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

  nvector<IndexType> cell2patch(PIH.npatches() << 3);
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
      if (habschon.find((ip << 3) + ile) != habschon.end())
        continue;
      habschon.insert((ip << 3) + ile);
      p1.push_back(ip);
      p2.push_back(ile);
    }
    BI.GetPatch().insert(make_pair(col, p1));
    BI.GetLocalPatch().insert(make_pair(col, p2));
  }
}
} // namespace Gascoigne
