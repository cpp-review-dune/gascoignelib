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

#include "coarsehierarchicalmesh3d.h"
#include "deletecells.h"
#include "facemanager.h"
#include "hierarchicalmesh3d.h"
#include "levelcomparer3d.h"
#include "quad.h"
#include "regular_update.h"
#include "set2vec.h"
#include "stlio.h"
#include "vecalgo.h"

#include "giota.h"
#include "localhierarchicalmesh3d.h"
#include <fstream>

using namespace std;

namespace Gascoigne {

void
LocalHierarchicalMesh3d::write_gup(const string& bname,
                                   const GlobalVector& u) const
{
  /*ofstream file(name.c_str());
  if(!file.is_open())
    {
      cerr << "HierarchicalMesh3d::write_inp()\n";
      cerr << "cannot open file " << name << endl;
      abort();
    }

  int nt = ncells()+nbquads();
  file <<nnodes()<<" "<<nt<<" "<<0<<" "<<0<<" "<<0<<endl;


  for(int i=0;i<nnodes();i++)
  {
    std::vector<double> new_coord;
    new_coord.resize(3);
                for(int c=0; c<u.ncomp(); c++)
                        {
                                new_coord[c]=u(i,c)+vertex3d(i)[c];
                        }

   file <<i<<" "<<new_coord << " " << endl;
   }

  for(int i=0;i<ncells();i++)
    {
      file << i<<" "<<0<<" hex "<<hex(i).vertex()<<endl;
    }
  for(int i=0;i<nbquads();i++)
    {
      file << i<<" "<< bquad(i).material() <<" quad "<<bquad(i).vertex()<<endl;
    }

*/

  string name = bname;
  name += ".gup";

  ofstream out(name.c_str());

  out << dimension() << " dimension" << endl;
  out << nnodes() << " vertexs" << endl;

  if (u.n() != nnodes()) {
    cout << "DDDDDDDDDDDDDDDDDDDDDDDDDD" << endl;
    abort();
  }
  for (int i = 0; i < nnodes(); i++) {
    std::vector<double> new_coord;
    new_coord.resize(3);
    for (int c = 0; c < u.ncomp(); c++) {
      new_coord[c] = u(i, c) + vertex3d(i)[c];
    }
    out << " " << new_coord << endl;
  }
  out << hexs.size() << " hexs" << endl;
  for (int i = 0; i < hexs.size(); i++) {
    out << hex(i);
  }
  out << QuadHang << endl;
  out << LineHang << endl;
  out << Bquads.size() << " boundaryquads" << endl;
  for (int i = 0; i < Bquads.size(); i++) {
    out << Bquads[i].material() << " " << Bquads[i] << endl;
  }
  out << endl << edges.size() << " edges" << endl;
  for (int i = 0; i < edges.size(); i++) {
    out << " " << edges[i];
  }
  out.close();
}

void
LocalHierarchicalMesh3d::inner_vertex_newton3d(const IntVector& vnew,
                                               const IntSet& CellRefList,
                                               const IntSet& adjustvertex)
{

  if (GetCurvedShapes().empty())
    return;

  // baue lokalen set auf um spaeter nicht alle hex zu justieren
  // Urspruengliche intention nur Randzellen anzupassen. Will ich hier ja gerade
  // nicht!!! IntSet Hexset; for (int i=0; i<Bquads.size(); i++)
  //  {
  //    int hi = Bquads[i].of_quad();
  //    Hexset.insert(hi);
  //  }

  std::array<int, 4> v;
  IntSetIt cp = CellRefList.begin();

  // Schleife uber alle Zellen, zaehlen der Nachbarzellen eines
  // Kantenmittelpunktes
  // Im Gebiet normalerweise 4 am Rand 2 (weitere Werte moegloch!)
  vector<int> RealBoundaryFaces;
  RealBoundaryFaces.clear();
  RealBoundaryFaces.resize(vertexs3d.size(), 0);

  for (int i = 0; i < CellRefList.size(); i++) {
    int hi = co2n[*cp++];
    const Hex& h = hex(hi);
    for (int e = 0; e < 12; ++e) {
      int ev = HexLaO.edge_vertex(h, e);
      RealBoundaryFaces[ev]++;
    }
  }

  // for (int i=0; i<Bquads.size(); i++)
  //  {
  //    int hi = Bquads[i].of_quad();
  //    Hexset.insert(hi);
  //  }
  //////////////////////////////////////////////////////////
  vector<vector<int>> Knotengegenueber;
  Knotengegenueber.resize(12);
  Knotengegenueber[0].push_back(2);
  Knotengegenueber[0].push_back(4);
  Knotengegenueber[1].push_back(3);
  Knotengegenueber[1].push_back(5);
  Knotengegenueber[2].push_back(0);
  Knotengegenueber[2].push_back(6);
  Knotengegenueber[3].push_back(7);
  Knotengegenueber[3].push_back(1);
  Knotengegenueber[4].push_back(0);
  Knotengegenueber[4].push_back(6);
  Knotengegenueber[5].push_back(1);
  Knotengegenueber[5].push_back(7);
  Knotengegenueber[6].push_back(4);
  Knotengegenueber[6].push_back(2);
  Knotengegenueber[7].push_back(3);
  Knotengegenueber[7].push_back(5);
  Knotengegenueber[8].push_back(9);
  Knotengegenueber[8].push_back(11);
  Knotengegenueber[9].push_back(8);
  Knotengegenueber[9].push_back(10);
  Knotengegenueber[10].push_back(9);
  Knotengegenueber[10].push_back(11);
  Knotengegenueber[11].push_back(8);
  Knotengegenueber[11].push_back(10);

  IntSet adjustadditionalvertex;
  IntSet localadjustvertex = adjustvertex;
  /////////////////////////////////////////////////////////
  int iter = 0;
  while (adjustadditionalvertex.size() != 0 || iter == 0) {
    iter++;
    adjustadditionalvertex.clear();
    cp = CellRefList.begin();
    for (int i = 0; i < CellRefList.size(); i++) {
      int hi = co2n[*cp++];
      const Hex& h = hex(hi);

      // if (Hexset.find(hi)==Hexset.end()) continue;
      // Abfrage ob es sich um eine Randzelle handelt

      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEU
      // edges
      for (int e = 0; e < 12; ++e) {
        int ev = HexLaO.edge_vertex(h, e);
        if (localadjustvertex.find(ev) != localadjustvertex.end())
          continue;

        if (adjustadditionalvertex.find(ev) == adjustadditionalvertex.end()) {
          std::array<int, 2> fe;
          HexLaO.globalvertices_of_edge(h, fe, e);
          vertexs3d[ev] = vertexs3d[fe[0]];
          vertexs3d[ev] += vertexs3d[fe[1]];
          vertexs3d[ev] *= 0.5;
        }
        // uberpruefen ob die Gegeneuberliegendende Kanten Knoten verschoben
        // worden sind
        for (int ii = 0; ii < 2; ii++) {
          int ev_gegenueber = HexLaO.edge_vertex(h, Knotengegenueber[e][ii]);
          if (localadjustvertex.find(ev_gegenueber) !=
              localadjustvertex.end()) {
            // berechnen wie weit der Knoten ev_gegenueber verschoben worden ist
            std::array<int, 2> fe;
            HexLaO.globalvertices_of_edge(h, fe, Knotengegenueber[e][ii]);
            Vertex3d verschiebungsvector = vertexs3d[ev_gegenueber];
            verschiebungsvector *= 2.0;
            verschiebungsvector -= vertexs3d[fe[0]];
            verschiebungsvector -= vertexs3d[fe[1]];
            verschiebungsvector *= 0.5;
            // Verhaltnis verschiebungsvector zu Zellhoehe:
            double lengthverschiebungsvector = 0;
            double lengthcell = 0;
            for (int iii = 0; iii < 3; iii++) {
              lengthverschiebungsvector +=
                verschiebungsvector[iii] * verschiebungsvector[iii];
              lengthcell += (0.5 * vertexs3d[fe[0]][iii] +
                             0.5 * vertexs3d[fe[1]][iii] - vertexs3d[ev][iii]) *
                            (0.5 * vertexs3d[fe[0]][iii] +
                             0.5 * vertexs3d[fe[1]][iii] - vertexs3d[ev][iii]);
            }
            // addieren von 0.25* verschiebungsvector auf den Knoten ev(knoten
            // grenzt im normalvall an 4 patches!!
            // bei Randknoten 0.5
            if (4 * lengthverschiebungsvector > lengthcell)
              verschiebungsvector *= 1. / RealBoundaryFaces[ev];
            else
              verschiebungsvector *=
                1. / RealBoundaryFaces[ev] *
                sqrt(sqrt(lengthverschiebungsvector) / sqrt(lengthcell));

            vertexs3d[ev] += verschiebungsvector;
            adjustadditionalvertex.insert(ev);
          }
        }
      }
    }
    localadjustvertex.insert(adjustadditionalvertex.begin(),
                             adjustadditionalvertex.end());
    cout << "iter: " << iter
         << "new adjust nodes: " << adjustadditionalvertex.size() << endl;
  }
  cp = CellRefList.begin();

  for (int i = 0; i < CellRefList.size(); i++) {
    int hi = co2n[*cp++];
    const Hex& h = hex(hi);

    // if (Hexset.find(hi)==Hexset.end()) { cout<<"hi not found"<<endl;
    // continue;} Abfrage ob es sich um eine Randzelle handelt

    // faces
    for (int f = 0; f < 6; ++f) {

      int fv = HexLaO.face_vertex(h, f);

      if (localadjustvertex.find(fv) != localadjustvertex.end())
        continue;

      std::array<int, 4> fe;
      HexLaO.LoadEdgeVerticesOfFace(h, f, fe);
      std::array<int, 4> fe_corner;
      HexLaO.globalvertices_of_face(h, fe_corner, f);

      ///////////////////////////////////////////////////////////////////////////////
      // If face is very anisotropic (||f[3]-f[1]||<<<||f[0]-f[2]|| ) the
      // algorithm will fail for slightly convex quads
      //*------------------*f[3]--------------*
      //|		       				 |
      //| *f[0]--------------*midface-------f[2]* |
      //|		 						  |
      //*------------------*f[1]--------------*
      // Then use for midface the midpoint between f[3] and f[1] which is always
      // in the face!

      /*fixarray<3,double> v1,v2;
      v1=vertexs3d[fe[0]]-vertexs3d[fe[2]];
      v2=vertexs3d[fe[1]]-vertexs3d[fe[3]];
      vector<double> nv(2,0);
      for(int i=0;i<4;i++)
      {
              nv[0]+=v1[i]*v1[i];
              nv[1]+=v2[i]*v2[i];
      }
      //Test if face is anisotropic
        if (nv[0]>25*nv[1])
              {
                      vertexs3d[fv].equ( 0.5, vertexs3d[fe[1]], 0.5,
      vertexs3d[fe[3]]);
              }
              else if (nv[1]>25*nv[0])
              {
                      vertexs3d[fv].equ( 0.5, vertexs3d[fe[0]],0.5,
      vertexs3d[fe[2]]);
              }
      /////////////////////////////////////////////////////////////////////////////
      else
              {
                      vertexs3d[fv].equ(0.25, vertexs3d[fe[0]],
                                      0.25, vertexs3d[fe[1]],
                                      0.25, vertexs3d[fe[2]],
                                      0.25, vertexs3d[fe[3]]);
       }*/
      Vertex3d zwischen1, zwischen2;
      zwischen1.equ(0.5,
                    vertexs3d[fe[0]],
                    0.5,
                    vertexs3d[fe[1]],
                    0.5,
                    vertexs3d[fe[2]],
                    0.5,
                    vertexs3d[fe[3]]);
      zwischen2.equ(-0.25,
                    vertexs3d[fe_corner[0]],
                    -0.25,
                    vertexs3d[fe_corner[1]],
                    -0.25,
                    vertexs3d[fe_corner[2]],
                    -0.25,
                    vertexs3d[fe_corner[3]]);
      vertexs3d[fv].equ(1.0, zwischen1, 1.0, zwischen2);
    }

    // middle
    int fv = HexLaO.middle_vertex(h);
    assert(localadjustvertex.find(fv) == localadjustvertex.end());
    std::array<int, 6> fe;
    HexLaO.LoadFaceVertices(h, fe);
    vector<Vertex3d> zwischen;
    zwischen.resize(10);
    for (int f = 0; f < 6; ++f) {
      std::array<int, 4> fe_corner;
      HexLaO.globalvertices_of_face(h, fe_corner, f);
      zwischen[f].equ(-0.25 / 12.0,
                      vertexs3d[fe_corner[0]],
                      -0.25 / 12.0,
                      vertexs3d[fe_corner[1]],
                      -0.25 / 12.0,
                      vertexs3d[fe_corner[2]],
                      -0.25 / 12.0,
                      vertexs3d[fe_corner[3]]);

      int fface = HexLaO.face_vertex(h, f);

      // Falls FaceMittelpunkt projeziert -> Korrektur um diesen Vektor
      if (localadjustvertex.find(fface) != localadjustvertex.end()) {
        std::array<int, 4> feedge;
        HexLaO.LoadEdgeVerticesOfFace(h, f, feedge);

        Vertex3d zwischen1, zwischen2;
        zwischen1.equ(0.5,
                      vertexs3d[feedge[0]],
                      0.5,
                      vertexs3d[feedge[1]],
                      0.5,
                      vertexs3d[feedge[2]],
                      0.5,
                      vertexs3d[feedge[3]]);
        zwischen2.equ(-0.25,
                      vertexs3d[fe_corner[0]],
                      -0.25,
                      vertexs3d[fe_corner[1]],
                      -0.25,
                      vertexs3d[fe_corner[2]],
                      -0.25,
                      vertexs3d[fe_corner[3]]);
        zwischen[f].equ(1.0,
                        zwischen[f],
                        -1.0 / 4.0,
                        zwischen1,
                        -1.0 / 4.0,
                        zwischen2,
                        +1.0 / 4.0,
                        vertexs3d[fe[f]]);
      }
    }

    zwischen[6].equ(
      1.0, zwischen[0], 1.0, zwischen[1], 1.0, zwischen[2], 1.0, zwischen[3]);
    zwischen[7].equ(1.0, zwischen[6], 1.0, zwischen[4], 1.0, zwischen[5]);

    zwischen[8].equ(0.25,
                    vertexs3d[fe[0]],
                    0.25,
                    vertexs3d[fe[1]],
                    0.25,
                    vertexs3d[fe[2]],
                    0.25,
                    vertexs3d[fe[3]]);
    zwischen[9].equ(
      1.0, zwischen[8], 0.25, vertexs3d[fe[4]], 0.25, vertexs3d[fe[5]]);

    vertexs3d[fv].equ(1.0, zwischen[7], 1.0, zwischen[9]);
  }
  ///////////////////////////////////////////////////////////////////////////////
  // If face is very anisotropic the algorithm will fail for slightly convex
  // quads
  /*fixarray<3,double> v1,v2,v3;
  v1=vertexs3d[fe[0]]-vertexs3d[fe[5]];
  v2=vertexs3d[fe[1]]-vertexs3d[fe[3]];
  v3=vertexs3d[fe[4]]-vertexs3d[fe[2]];
  vector<double> nv(3,0);
    int sort[]={0,1,2};
  for(int i=0;i<4;i++)
  {
          nv[0]+=v1[i]*v1[i];
          nv[1]+=v2[i]*v2[i];
    nv[2]+=v3[i]*v3[i];
  }

// Sort of the distance
if(nv[0]>nv[1]) {double zw=nv[0]; nv[0]=nv[1]; nv[1]=zw; int zwsort=sort[0];
sort[0]=sort[1];sort[1]=zwsort;} if(nv[1]>nv[2]) {double zw=nv[1]; nv[1]=nv[2];
nv[2]=zw; int zwsort=sort[1]; sort[1]=sort[2];sort[2]=zwsort;} if(nv[0]>nv[1])
{double zw=nv[0]; nv[0]=nv[1]; nv[1]=zw; int zwsort=sort[0];
sort[0]=sort[1];sort[1]=zwsort;}

Vertex3d zw1;
//Test if  quad is anisotropic
vector<double> weight(3,0);
if(nv[0]*25<nv[2] && nv[1]*25<nv[2])
  {
          for(int ii=0;ii<3;ii++) weight[ii]=0.25;
          weight[sort[2]]=0;
  }
else if(nv[0]*25<nv[2])
  {
          for(int ii=0;ii<3;ii++) weight[ii]=0.5;
          weight[sort[1]]=0;
          weight[sort[2]]=0;
  }
else
  {
          for(int ii=0;ii<3;ii++) weight[ii]=1./6.;
  }

  vertexs3d[fv].equ(weight[0],vertexs3d[fe[0]],
                          weight[0],vertexs3d[fe[5]],
                          weight[1],vertexs3d[fe[1]],
                          weight[1],vertexs3d[fe[3]],
                          weight[2],vertexs3d[fe[4]],
                          weight[2],vertexs3d[fe[2]]);
*/

  //////////////////////////////////////////////////////////////////
  /*else
  {
    for (int i=0;i<6;++i) vertexs3d[fv]+=vertexs3d[fe[i]];
    vertexs3d[fv]*=1./6.;
  }*/

  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEU

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ALT
  //       for (int face=0; face<6; face++)
  // 	{
  // 	  HexLaO.LoadEdgeVerticesOfFace(h,face,v);
  // 	  int mv = HexLaO.face_vertex(h,face);
  // 	  if (adjustvertex.find(mv)!=adjustvertex.end()) continue;
  // 	  new_face_vertex3d(mv,v);
  // 	}
  // //       fixarray<6,int> w;
  // //       int mv = HexLaO.middle_vertex(h);
  // //       HexLaO.LoadFaceVertices(h,w);
  // //       new_vertex3d(mv,w);
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ALT
}

} // namespace Gascoigne

/*---------------------------------------------------*/
