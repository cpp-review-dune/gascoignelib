/**
*
* Copyright (C) 2004, 2005, 2006, 2007 by the Gascoigne 3D authors
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


#ifndef __hierarchicalmesh3d_h
#define __hierarchicalmesh3d_h

#include "boundaryfunction.h"
#include "boundaryquad.h"
#include "hangcontainer3d.h"
#include "hex.h"
#include "hexlawandorder.h"
#include "hierarchicalmesh.h"
#include "metis.h"
#include "vertex.h"
#include <fstream>


/// Hier muessen wir was machen. Das ist zuviel in der h-Datei...
 #include <../../boost/config.hpp>
 #include <vector>
 #include <iostream>
 #include <../../boost/graph/adjacency_list.hpp>
 #include <../../boost/graph/cuthill_mckee_ordering.hpp>
 #include <../../boost/graph/properties.hpp>
 #include <../../boost/graph/bandwidth.hpp>

/*---------------------------------------------------*/

namespace Gascoigne
{
  class HierarchicalMesh3d : public HierarchicalMesh
  {
  protected:
    /*  typedef  */

    typedef std::vector<Vertex3d> VertexVec3d;
    typedef BoundaryCell<4> BoundaryQuad;
    typedef std::vector<Hex> HexVec;
    typedef std::vector<BoundaryQuad> BQuadVec;
    typedef HangList<2> LineHangList;
    typedef HangList<4> QuadHangList;
    typedef BoundaryFunction<3> BoundaryFunction3d;
    typedef std::map<int, std::array<int, 8>> HexChilds;

    /*  Data  */

    CurvedShapes<3> _curvedshapes;

    VertexVec3d vertexs3d;

    /* info fuer interpolation auf neues gitter */
    HexChilds hexchildsofdeleted;
    HexVec hexs;
    BQuadVec Bquads;
    LineHangList LineHang;
    QuadHangList QuadHang;
    HexLawAndOrder HexLaO;
    std::map<int, int> hexofcurved;

    /*  Functionen  */
    int Vater(const int i) const;
    IntVector Nachkommen(const int i) const;
    IntVector Geschwister(const int i) const;
    IntVector Kinder(const int i) const;

    void post_refine3d();

    void delete_vertexs3d(const IntVector &);

    void new_edge_vertex3d(int, const EdgeVector &);
    void new_face_vertex3d(int, const FaceVector &);
    void new_vertex3d(int, const std::array<int, 6> &);

    void check_mesh3d() const;

    std::pair<bool, tint> check_inp(const std::string &);
    std::pair<int, int> GetBoundaryInformation(int i) const;

    void init_quad(BoundaryQuad &);

    void build_neighbours() const;

    void prepare3d(const IntVector &, const IntVector &, IntSet &, IntSet &);
    void new_hexs(const HangContainer3d &,
                  const IntVector &,
                  const IntVector &,
                  int,
                  const IntSet &);
    void ghost_fill_neighbours2d();
    void ghost_fill_neighbours3d();
    void UpdateHangs(HangContainer3d &hangset,
                     const IntSet &cellref,
                     const IntSet &cellcoarse);
    void FaceCoarse(HangContainer3d &, const IntSet &) const;
    void FaceRefine(HangContainer3d &, const IntSet &) const;
    void UpdateHangingEdges(HangContainer3d &hangset,
                            const IntSet &cellref,
                            const IntSet &cellcoarse) const;
    void
    boundary_prepare3d(IntSet &, IntSet &, IntSet &, const HangContainer3d &);
    void new_boundary3d(IntSet &, IntSet &, IntSet &);
    void new_vertexs3d(HangContainer3d &, const IntVector &, const IntSet &);
    void basic_refine3d(HangContainer3d &, const IntSet &, const IntSet &);
    void basic_fill_neighbours3d();
    void boundary_newton3d(IntSet &);
    virtual void
    inner_vertex_newton3d(const IntVector &, const IntSet &, const IntSet &);
    void update_boundary_data3d(const IntSet &);
    void new_bquads(const IntVector &, const IntVector &, const IntSet &);
    void new_middle_vertex3d(int, int);

    int
    regular_grid3d_one(IntSet &, IntVector &, const IntSet &, const IntSet &);
    int regular_grid3d_one(IntVector &,
                           IntVector &,
                           const IntSet &,
                           const IntSet &);
    int regular_grid3d_two(IntVector &, const IntSet &);
    int regular_grid3d_three_refine(IntSet &) const;
    int regular_grid3d_three_coarse(IntSet &, IntSet &) const;

    void GetMinMaxLevels(IntVector &maxi,
                         IntVector &mini,
                         const IntSet &CellRef) const;

    void init_edges3d();
    void LoadFathers3d(IntVector &v) const;

    void _refine3d(IntSet &, IntSet &, const IntVector &, const IntVector &);
    void
    FillNeighbourFaces(const Hex &father, const FaceVector &Face, int rneigh);
    void FillNeighbourFaces(int M, int S, const FaceVector &Face);
    void InitHexOfCurved();
    int FindPatchDepth() const;
    void FillVertexLevels(IntVector &dst) const;
    void RefineCoarseNodes(IntSet &dst,
                           const IntVector &refnodes,
                           const IntVector &vertexlevel) const;
    void VertexToCells(IntVector &dst,
                       const IntSet &src,
                       const IntVector &vertexlevel) const;
    void VertexToCellsCoarsening(IntVector &dst,
                                 const IntSet &src,
                                 const IntVector &vertexlevel) const;
    void recursive_childs(int q, IntVector &ref, int d) const;

    struct sort_pred
    {
      bool operator()(const std::pair<Vertex3d, double> &left,
                      const std::pair<Vertex3d, double> &right)
      {
        return left.second < right.second;
      }
    };

  public:
    HierarchicalMesh3d();
    HierarchicalMesh3d(const HierarchicalMesh3d &H);
    HierarchicalMesh3d &operator=(const HierarchicalMesh3d &H);
    HierarchicalMesh3d(const ParamFile *paramfile);
    ~HierarchicalMesh3d()
    {
      GetCurvedShapes().clear();
    }

    std::string GetName() const
    {
      return "HierarchicalMesh3d";
    }

    /*  Zugriff  */

    int dimension() const
    {
      return 3;
    }

    int nnodes() const
    {
      return vertexs3d.size();
    }
    int ncells() const
    {
      return hexs.size();
    }
    int nbquads() const
    {
      return Bquads.size();
    }

    int nodes_per_cell(int i) const
    {
      return 8;
    }
    int VtkType(int i) const
    {
      return 12;
    }

    const CurvedShapes<3> &GetCurvedShapes() const
    {
      return _curvedshapes;
    }
    CurvedShapes<3> &GetCurvedShapes()
    {
      return _curvedshapes;
    }

    const VertexVec3d &GetVertexVector() const
    {
      return vertexs3d;
    }
    VertexVec3d &GetVertexVector()
    {
      return vertexs3d;
    }

    const Vertex3d &vertex3d(int i) const
    {
      return vertexs3d[i];
    }

    const Hex &hex(int i) const
    {
      return hexs[i];
    }
    const BoundaryQuad &bquad(int i) const
    {
      return Bquads[i];
    }

    int vertex_of_cell(int i, int ii) const
    {
      return hexs[i].vertex(ii);
    }
    int vertex_of_bquad(int i, int ii) const
    {
      return Bquads[i].vertex(ii);
    }
    int face_of_hex(int i, int ii) const
    {
      return hexs[i].edge(ii);
    }
    int level(int i) const
    {
      return hexs[i].level();
    }
    bool sleep(int i) const
    {
      return hexs[i].sleep();
    }

    int child(int i, int ii) const
    {
      return hexs[i].child(ii);
    }
    int nchilds(int i) const
    {
      return hexs[i].nchilds();
    }

    const HexLawAndOrder &HexLawOrder() const
    {
      return HexLaO;
    }
    const LineHangList &linehang() const
    {
      return LineHang;
    }
    const QuadHangList &quadhanglist() const
    {
      return QuadHang;
    }
    const BoundaryFunction3d *quad_shape(int i) const;

    const std::vector<BoundaryQuad> &quad_list() const
    {
      return Bquads;
    }

    const VertexVec3d &vertex3d() const
    {
      return vertexs3d;
    }
    const HexVec &hex() const
    {
      return hexs;
    }
    const BQuadVec &bquad() const
    {
      return Bquads;
    }
    const QuadHangList &quadhang() const
    {
      return QuadHang;
    }
    const std::map<int, int> &GetHexOfCurved() const
    {
      return hexofcurved;
    }

    /*  Functionen  */

    void write(const std::string &) const;
    void write_gup(const std::string &) const;
    void write_gip(const std::string &) const;

    void WriteAll(const std::string &) const;

    void write_inp(const std::string &) const;
    void read_inp(const std::string &);
    void read_gup(const std::string &);
    void read_gip(const std::string &);

    void global_coarse3d();

    void refine(const IntVector &, const IntVector &);
    void patch_refine(IntVector &, IntVector &);
    //  int    smooth_edges();
    void FillAllBoundaryLines();

    pint EdgeNeighbour(int i, int e) const;

    int NodeOnFace(int e) const;
    std::array<int, 4> ChildrenOfFace(int e) const;

    void GetVertexesOfFace(std::array<int, 4> &, int) const;
    void GetVertexesOfFace(std::array<int, 5> &, int) const;
    void GetAwakePatchs(std::set<int> &) const;
    void GetAwakeCells(std::set<int> &) const;
    void ConstructQ2PatchMesh(IntVector &pm) const;
    IntVector ConstructQ4Patch(int c) const;
    std::set<int> GetColors() const;


    int nactivedescendants(int i) const;
    IntVector GetVertices(int c) const;

    int GetBoundaryCellOfCurved(int iq) const
    {
      std::map<int, int>::const_iterator p = hexofcurved.find(iq);
      if (p != hexofcurved.end())
        return p->second;
      return -1;
    }
    void Testing();
    int neighbour(int c, int le) const;
    int neighbour_neighbour(int c, int le) const;

    void AddShape(int col, BoundaryFunction<3> *f)
    {
      GetCurvedShapes().AddShape(col, f);
    }
    void Sort()
    {
      //abort();
      
      // connectivity - patch based
       std::set<int> patches;
       GetAwakePatchs(patches);
       idx_t nn = nnodes();
       std::vector<std::set<int>> n2n_set(nn);

       for (auto it : patches)
       {
         // list of vertices in patch
         const Hex &H = hex(it);

         std::set<int> nh; // vertices in patch
         assert(H.nchilds() == 8);
         for (int c = 0; c < H.nchilds(); ++c)
         {
           const Hex &C = hex(H.child(c));
           for (int n = 0; n < 8; ++n)
             nh.insert(C[n]);
         }
         assert(nh.size() == 27);
         for (auto i1 : nh)
           for (auto i2 : nh)
             n2n_set[i1].insert(i2);
       }

       for (int i = 0; i < n2n_set.size(); ++i)
         assert(n2n_set[i].size() > 1);

       std::ofstream matrix_log("unsortedmatrix.txt");
       for (auto itn2n_set : n2n_set)
         matrix_log << itn2n_set << std::endl;
       matrix_log.close();

      // copy coupling structure to METIS format
       std::vector<idx_t> adj(nn + 1, 0);
       std::vector<idx_t> adjncy;

       int count = 0;
       adj[0] = 0;
       for (int r = 0; r < n2n_set.size(); ++r) // loop over 'rows'
       {
         for (auto it : n2n_set[r]) // loop over 'cols'
         {
           if (it == r)
             continue;
           ++count;
           adjncy.push_back(it);
         }
         adj[r + 1] = count;
       }

    //   /*
    //       std::vector<idx_t> iperm(nn);
    //       std::vector<idx_t> perm(nn);
    //       for (int i=0;i<nn;++i)
    //         perm[i]=i;


    //       //////////// METIS!!!!
    //       idx_t options[METIS_NOPTIONS];
    //       METIS_SetDefaultOptions(options);
    //       options[METIS_OPTION_NUMBERING] = 0;
    //       //////////// METIS!!!!


    //       METIS_NodeND(&nn,&adj[0],&adjncy[0],NULL,
    //          &options[0], &perm[0], &iperm[0]);
    //       //////////// CMC
    //       //abort();
    //   */
       using namespace boost;
       using namespace std;
       typedef adjacency_list<vecS,
                              vecS,
                              undirectedS,
                              property<vertex_color_t,
                                       default_color_type,
                                       property<vertex_degree_t, int>>>
           Graph;
       typedef graph_traits<Graph>::vertex_descriptor Vertex;
       typedef graph_traits<Graph>::vertices_size_type size_type;

      Graph G(nn);
       for (int i = 0; i < nn; i++)
       {
         for (auto j : n2n_set[i])
           if (j <= i)
             add_edge(i, j, G);
       }

       graph_traits<Graph>::vertex_iterator ui, ui_end;

      property_map<Graph, vertex_degree_t>::type deg = get(vertex_degree, G);
       for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
         deg[*ui] = degree(*ui, G);

       property_map<Graph, vertex_index_t>::type

           index_map = get(vertex_index, G);


       std::vector<Vertex> perm(num_vertices(G));
       std::vector<size_type> iperm(num_vertices(G));

    //   //    Vertex s = vertex(6, G);
    //   // reverse cuthill_mckee_ordering
    //   // cuthill_mckee_ordering(G, s, perm.rbegin(), get(vertex_color, G),
    //   //                       get(vertex_degree, G));
       cuthill_mckee_ordering(
           G, perm.rbegin(), get(vertex_color, G), get(vertex_degree, G));
       for (size_type c = 0; c != perm.size(); ++c)
        iperm[index_map[perm[c]]] = c;


       assert(nnodes() == perm.size());
       // vertices
       VertexVec3d newVV(nnodes());
       VertexVec3d &origVV = GetVertexVector();
       for (int i = 0; i < nnodes(); ++i)
         newVV[i] = origVV[perm[i]];
       for (int i = 0; i < nnodes(); ++i)
         origVV[i] = newVV[i];

       for (int i = 0; i < vo2n.size(); ++i)
         vo2n[i] = iperm[vo2n[i]];

       // hex
       for (int h = 0; h < hexs.size(); ++h)
         for (int n = 0; n < 8; ++n)
           hexs[h][n] = iperm[hexs[h][n]];

       // cell
       for (int h = 0; h < Bquads.size(); ++h)
         for (int n = 0; n < 4; ++n)
           Bquads[h][n] = iperm[Bquads[h][n]];


       check_mesh3d();


      n2n_set.clear();
       n2n_set.resize(nn);
       for (auto it : patches)
       {
         // list of vertices in patch
         const Hex &H = hex(it);

         std::set<int> nh; // vertices in patch
         assert(H.nchilds() == 8);
         for (int c = 0; c < H.nchilds(); ++c)
         {
           const Hex &C = hex(H.child(c));
           for (int n = 0; n < 8; ++n)
             nh.insert(C[n]);
         }
         assert(nh.size() == 27);
         for (auto i1 : nh)
           for (auto i2 : nh)
             n2n_set[i1].insert(i2);
       }

      for (int i = 0; i < n2n_set.size(); ++i)
         assert(n2n_set[i].size() > 1);

       std::ofstream matrix_sorted_log("sortedmatrix.txt");
       for (auto itn2n_set : n2n_set)
         matrix_sorted_log << itn2n_set << std::endl;
       matrix_sorted_log.close();
    }
  };
}

/*---------------------------------------------------*/

#endif
