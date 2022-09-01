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

#ifndef __hierarchicalmesh2d_h
#define __hierarchicalmesh2d_h

#include "boundaryfunction.h"
#include "boundaryline.h"
#include "hangcontainer2d.h"
#include "hierarchicalmesh.h"
#include "quad.h"
#include "quadlawandorder.h"
#include "vertex.h"

/*---------------------------------------------------*/

namespace Gascoigne {
class HierarchicalMesh2d : public HierarchicalMesh
{
protected:
  typedef std::vector<Vertex2d> VertexVec2d;
  typedef std::vector<Quad> QuadVec;
  typedef BoundaryFunction<2> BoundaryFunction2d;
  typedef BoundaryCell<2> BoundaryLine;
  typedef std::vector<BoundaryLine> BLineVec;
  typedef HangList<2> LineHangList;

  typedef triple<IndexType, IndexType, IndexType> tint;
  typedef std::map<std::array<IndexType, 2>, HierarchicalMesh2d::BoundaryLine>
    HangBLList;

  /*  Data  */

  CurvedShapes<2> _curvedshapes;

  VertexVec2d vertexs2d;

  QuadVec quads;
  BLineVec Blines;
  LineHangList LineHang;
  QuadLawAndOrder QuadLaO;
  std::map<IndexType, IndexType> quadofcurved;

  /*  Functionen  */

  void post_refine2d();

  void delete_vertexs2d(const IndexVector&);

  void new_edge_vertex2d(IndexType, const EdgeVector&);
  void new_face_vertex2d(IndexType, const FaceVector&);

  void check_mesh2d() const;
  void prepare2d(const IndexVector&, const IndexVector&, IndexSet&, IndexSet&);
  std::pair<bool, tint> check_inp(const std::string&);
  void ghost2d(HangContainer2d&, const IndexSet&, const IndexSet&);
  void ghostglobalcoarse(HangContainer2d&, const IndexSet&);
  void ghost_fill_neighbours2d();
  void basic_fill_neighbours2d();
  void new_vertexs2d(HangContainer2d&, const IndexVector&, const IndexSet&);
  void new_quads(const HangContainer2d&,
                 const IndexVector&,
                 const IndexVector&,
                 IndexType,
                 const IndexSet&);

  void change_hangs2d(const IndexVector&, const IndexVector&);
  void change_vertexs2d(const IndexVector&);
  void change_quads2d(const IndexVector&, const IndexVector&);
  void boundary_prepare2d(IndexSet&,
                          IndexSet&,
                          IndexSet&,
                          const HangContainer2d&);
  void new_boundary2d(IndexSet&, IndexSet&, IndexSet&);

  void basic_refine2d(HangContainer2d&, const IndexSet&, const IndexSet&);

  void init_line(BoundaryLine&);
  void new_lines(const IndexVector&, const IndexVector&, const IndexSet&);
  void boundary_newton2d();
  void inner_vertex_newton2d(const IndexVector&, const IndexSet&);
  void update_boundary_data2d(const IndexSet&);

  IndexType regular_grid2d_one(IndexSet&,
                               IndexVector&,
                               IndexSet&,
                               IndexSet&) const;
  IndexType regular_grid2d_two(IndexSet&, IndexSet&) const;
  IndexType regular_grid2d_three(IndexSet&, IndexSet&) const;
  IndexType regular_grid2d_three_refine(IndexSet&) const;
  IndexType regular_grid2d_three_coarse(IndexSet&, IndexSet&) const;

  void GetMinMaxLevels(IndexVector& maxi,
                       IndexVector& mini,
                       const IndexSet& CellRef) const;
  void init_edges2d();

  void LoadFathers(IndexVector& v) const;

  void _refine2d(IndexSet&, IndexSet&, const IndexVector&, const IndexVector&);
  void InitQuadOfCurved();
  IndexType FindPatchDepth() const;
  void FillVertexLevels(IndexVector& dst) const;
  void RefineCoarseNodes(IndexSet& dst,
                         const IndexVector& refnodes,
                         const IndexVector& vertexlevel) const;
  void VertexToCells(IndexVector& dst,
                     const IndexSet& src,
                     const IndexVector& vertexlevel) const;
  void VertexToCellsCoarsening(IndexVector& dst,
                               const IndexSet& src,
                               const IndexVector& vertexlevel) const;
  void recursive_childs(IndexType q, IndexVector& ref, IndexType d) const;

  struct sort_pred
  {
    bool operator()(const std::pair<Vertex3d, double>& left,
                    const std::pair<Vertex3d, double>& right)
    {
      return left.second < right.second;
    }
  };

public:
  HierarchicalMesh2d();
  HierarchicalMesh2d(const HierarchicalMesh2d& H);
  HierarchicalMesh2d& operator=(const HierarchicalMesh2d& H);
  HierarchicalMesh2d(const ParamFile& paramfile);
  ~HierarchicalMesh2d() { GetCurvedShapes().clear(); }

  std::string GetName() const { return "HierarchicalMesh2d"; }

  /*  Zugriff  */

  IndexType dimension() const { return 2; }
  IndexType nnodes() const { return vertexs2d.size(); }
  IndexType ncells() const { return quads.size(); }
  IndexType nblines() const { return Blines.size(); }
  IndexType nodes_per_cell(IndexType i) const { return 4; }
  IndexType VtkType(IndexType i) const { return 9; }

  const CurvedShapes<2>& GetCurvedShapes() const { return _curvedshapes; }
  CurvedShapes<2>& GetCurvedShapes() { return _curvedshapes; }

  const VertexVec2d& GetVertexVector() const { return vertexs2d; }
  VertexVec2d& GetVertexVector() { return vertexs2d; }

  const Vertex2d& vertex2d(IndexType i) const { return vertexs2d[i]; }

  const Quad& quad(IndexType i) const { return quads[i]; }
  const BoundaryLine& bline(IndexType i) const { return Blines[i]; }
  std::pair<IndexType, IndexType> GetBoundaryInformation(IndexType i) const;

  IndexType vertex_of_cell(IndexType i, IndexType ii) const
  {
    return quads[i].vertex(ii);
  }
  IndexType vertex_of_bline(IndexType i, IndexType ii) const
  {
    return Blines[i].vertex(ii);
  }
  IndexType edge_of_quad(IndexType i, IndexType ii) const
  {
    return quads[i].edge(ii);
  }
  IndexType level(IndexType i) const { return quads[i].level(); }
  bool sleep(IndexType i) const { return quads[i].sleep(); }

  IndexType child(IndexType i, IndexType ii) const
  {
    return quads[i].child(ii);
  }
  IndexType nchilds(IndexType i) const { return quads[i].nchilds(); }

  IndexType QuadNeighbour(const Quad&, IndexType) const;

  const QuadLawAndOrder& QuadLawOrder() const { return QuadLaO; }
  const LineHangList& linehanglist() const { return LineHang; }
  const BoundaryFunction2d* line_shape(IndexType i) const;

  const std::vector<BoundaryLine>& line_list() const { return Blines; }

  const VertexVec2d& vertex2d() const { return vertexs2d; }
  const QuadVec& quad() const { return quads; }
  const BLineVec& bline() const { return Blines; }
  const LineHangList& linehang() const { return LineHang; }
  const std::map<IndexType, IndexType>& GetQuadOfCurved() const
  {
    return quadofcurved;
  }

  /*  Functionen  */

  IndexType Vater(const IndexType i) const;
  IndexVector Nachkommen(const IndexType i) const;
  IndexVector Geschwister(const IndexType i) const;
  IndexVector Kinder(const IndexType i) const;
  IndexType nactivedescendants(IndexType i) const;
  IndexVector GetVertices(IndexType c) const;

  void write(const std::string&) const;
  void write_gup(const std::string&) const;
  void write_gip(const std::string&) const;

  void WriteAll(const std::string&) const;

  void write_inp(const std::string&) const;
  void read_inp(const std::string&);
  void write_vtk(const std::string&) const;

  void read_gup(const std::string&);
  void read_gip(const std::string&);

  void global_coarse();

  void refine(const IndexVector&, const IndexVector&);
  void patch_refine(IndexVector&, IndexVector&);
  IndexType smooth_edges();
  void FillAllBoundaryLines();

  pint EdgeNeighbour(IndexType i, IndexType e) const;
  void VertexNeighbours2d(std::set<IndexType>&, IndexType i) const;

  IndexType NodeOnEdge(IndexType e) const;
  std::array<IndexType, 2> ChildrenOfEdge(IndexType e) const;

  void GetVertexesOfEdge(std::array<IndexType, 3>&, IndexType) const;
  void GetVertexesOfEdge(std::array<IndexType, 2>&, IndexType) const;
  void GetAwakePatchs(std::set<IndexType>&) const;
  void GetAwakeCells(std::set<IndexType>&) const;
  void ConstructQ2PatchMesh(IndexVector& pm) const;
  IndexVector ConstructQ4Patch(IndexType c) const;
  std::set<IndexType> GetColors() const;
  int GetBoundaryCellOfCurved(IndexType iq) const
  {
    std::map<IndexType, IndexType>::const_iterator p = quadofcurved.find(iq);
    if (p != quadofcurved.end())
      return p->second;
    return -1;
  }

  std::set<IndexType> CellNeighbours(IndexType i) const;

  IndexType neighbour(IndexType c, IndexType le) const;
  void FillVolumes(DoubleVector& vol) const;
  void writeq2(const IndexVector& a,
               const std::vector<IndexType>& b,
               IndexType np) const;

  void AddShape(IndexType col, BoundaryFunction<2>* f)
  {
    GetCurvedShapes().AddShape(col, f);
  }
  void Sort()
  {
    abort();

    //   // connectivity - patch based
    //    std::set<IndexType> patches;
    //    GetAwakePatchs(patches);
    //    idx_t nn = nnodes();
    //    std::vector<std::set<IndexType>> n2n_set(nn);

    //    for (auto it : patches)
    //    {
    //      // list of vertices in patch
    //      const Quad &Q = quad(it);

    //      std::set<IndexType> nh; // vertices in patch
    //      assert(Q.nchilds() == 4);
    //      for (IndexType c = 0; c < Q.nchilds(); ++c)
    //      {
    //        const Quad &C = quad(Q.child(c));
    //        for (IndexType n = 0; n < 4; ++n)
    //          nh.insert(C[n]);
    //      }
    //      assert(nh.size() == 9);
    //      for (auto i1 : nh)
    //        for (auto i2 : nh)
    //          n2n_set[i1].insert(i2);
    //    }

    //    for (IndexType i = 0; i < n2n_set.size(); ++i)
    //      assert(n2n_set[i].size() > 1);

    //    std::ofstream matrix_log("unsortedmatrix.txt");
    //    for (auto itn2n_set : n2n_set)
    //      matrix_log << itn2n_set << std::endl;
    //    matrix_log.close();

    //   // copy coupling structure to METIS format
    //    std::vector<idx_t> adj(nn + 1, 0);
    //    std::vector<idx_t> adjncy;

    //    IndexType count = 0;
    //    adj[0] = 0;
    //    for (IndexType r = 0; r < n2n_set.size(); ++r) // loop over 'rows'
    //    {
    //      for (auto it : n2n_set[r]) // loop over 'cols'
    //      {
    //        if (it == r)
    //          continue;
    //        ++count;
    //        adjncy.push_back(it);
    //      }
    //      adj[r + 1] = count;
    //    }

    // //   /*
    // //       std::vector<idx_t> iperm(nn);
    // //       std::vector<idx_t> perm(nn);
    // //       for (IndexType i=0;i<nn;++i)
    // //         perm[i]=i;

    // //       //////////// METIS!!!!
    // //       idx_t options[METIS_NOPTIONS];
    // //       METIS_SetDefaultOptions(options);
    // //       options[METIS_OPTION_NUMBERING] = 0;
    // //       //////////// METIS!!!!

    // //       METIS_NodeND(&nn,&adj[0],&adjncy[0],NULL,
    // //          &options[0], &perm[0], &iperm[0]);
    // //       //////////// CMC
    // //       //abort();
    // //   */
    //    using namespace boost;
    //    using namespace std;
    //    typedef adjacency_list<vecS,
    //                           vecS,
    //                           undirectedS,
    //                           property<vertex_color_t,
    //                                    default_color_type,
    //                                    property<vertex_degree_t, IndexType>>>
    //        Graph;
    //    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    //    typedef graph_traits<Graph>::vertices_size_type size_type;

    //   Graph G(nn);
    //    for (IndexType i = 0; i < nn; i++)
    //    {
    //      for (auto j : n2n_set[i])
    //        if (j <= i)
    //          add_edge(i, j, G);
    //    }

    //    graph_traits<Graph>::vertex_iterator ui, ui_end;

    //   property_map<Graph, vertex_degree_t>::type deg = get(vertex_degree, G);
    //    for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
    //      deg[*ui] = degree(*ui, G);

    //    property_map<Graph, vertex_index_t>::type

    //        index_map = get(vertex_index, G);

    //    std::vector<Vertex> perm(num_vertices(G));
    //    std::vector<size_type> iperm(num_vertices(G));

    // //   //    Vertex s = vertex(6, G);
    // //   // reverse cuthill_mckee_ordering
    // //   // cuthill_mckee_ordering(G, s, perm.rbegin(), get(vertex_color, G),
    // //   //                       get(vertex_degree, G));
    //    cuthill_mckee_ordering(
    //        G, perm.rbegin(), get(vertex_color, G), get(vertex_degree, G));
    //    for (size_type c = 0; c != perm.size(); ++c)
    //     iperm[index_map[perm[c]]] = c;

    //    assert(nnodes() == perm.size());
    //    // vertices
    //    VertexVec2d newVV(nnodes());
    //    VertexVec2d &origVV = GetVertexVector();
    //    for (IndexType i = 0; i < nnodes(); ++i)
    //      newVV[i] = origVV[perm[i]];
    //    for (IndexType i = 0; i < nnodes(); ++i)
    //      origVV[i] = newVV[i];

    //    for (IndexType i = 0; i < vo2n.size(); ++i)
    //      	vo2n[i] = iperm[vo2n[i]];

    //    // hex
    //    for (IndexType q = 0; q < quads.size(); ++q)
    //      for (IndexType n = 0; n < 4; ++n)
    //       	quads[q][n] = iperm[quads[q][n]];

    //    // cell
    //    for (IndexType h = 0; h < Blines.size(); ++h)
    //      for (IndexType n = 0; n < 2; ++n)
    //        Blines[h][n] = iperm[Blines[h][n]];

    //    check_mesh2d();

    //   n2n_set.clear();
    //    n2n_set.resize(nn);
    //    for (auto it : patches)
    //    {
    //      // list of vertices in patch
    //      const Quad &Q = quad(it);

    //      std::set<IndexType> nh; // vertices in patch
    //      assert(Q.nchilds() == 4);
    //      for (IndexType c = 0; c < Q.nchilds(); ++c)
    //      {
    //        const Quad &C = quad(Q.child(c));
    //        for (IndexType n = 0; n < 4; ++n)
    //          nh.insert(C[n]);
    //      }
    //      assert(nh.size() == 9);
    //      for (auto i1 : nh)
    //        for (auto i2 : nh)
    //          n2n_set[i1].insert(i2);
    //    }

    //   for (IndexType i = 0; i < n2n_set.size(); ++i)
    //      assert(n2n_set[i].size() > 1);

    //    std::ofstream matrix_sorted_log("sortedmatrix.txt");
    //    for (auto itn2n_set : n2n_set)
    //      matrix_sorted_log << itn2n_set << std::endl;
    //    matrix_sorted_log.close();
  }
};
} // namespace Gascoigne

/*---------------------------------------------------*/

#endif
