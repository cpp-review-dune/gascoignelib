#ifndef __hierarchicalmesh3d_h
#define __hierarchicalmesh3d_h

#include  "vertex.h" 
#include  "boundaryquad.h" 
#include  "hex.h" 
#include  "hexlawandorder.h" 
#include  "boundaryfunction.h"
#include  "hangcontainer3d.h" 
#include  "hierarchicalmesh.h"

/*---------------------------------------------------*/

class HierarchicalMesh3d : public HierarchicalMesh
{
  protected :
  
  /*  typedef  */

  typedef  std::vector<Vertex3d>       VertexVec3d;
  typedef  BoundaryCell<4>          BoundaryQuad;
  typedef  std::vector<Hex>            HexVec;
  typedef  std::vector<BoundaryQuad>   BQuadVec;
  typedef  HangList<2>              LineHangList;
  typedef  HangList<4>              QuadHangList;
  typedef  BoundaryFunction<3>      BoundaryFunction3d;
  typedef  std::map<int,fixarray<8,int> >  HexChilds;

  /*  Data  */

  CurvedShapes<3>    _curvedshapes;

  VertexVec3d        vertexs3d; 

  /* info fuer interpolation auf neues gitter */
  HexChilds          hexchildsofdeleted;
  HexVec             hexs;
  BQuadVec           Bquads;
  LineHangList       LineHang;
  QuadHangList       QuadHang;
  HexLawAndOrder     HexLaO;
  std::map<int,int>       hexofcurved;

  /*  Functionen  */
  int    Vater(const int i) const;
  Gascoigne::IntVector Nachkommen(const int i) const;
  Gascoigne::IntVector Geschwister(const int i) const;
  Gascoigne::IntVector Kinder     (const int i) const;
  
  void post_refine3d();

  void delete_vertexs3d(const Gascoigne::IntVector&);

  void new_edge_vertex3d(int, const EdgeVector&);
  void new_face_vertex3d(int, const FaceVector&);
  void new_vertex3d     (int, const fixarray<6,int>&);

  void check_mesh3d() const;

  std::pair<bool,tint>  check_inp(const std::string&);
  std::pair<int,int>  GetBoundaryInformation(int i) const;

  void init_quad           (BoundaryQuad&);

  void  build_neighbours() const;

  void prepare3d  (const Gascoigne::IntVector&, const Gascoigne::IntVector&, Gascoigne::IntSet&, Gascoigne::IntSet&);
  void new_hexs   (const HangContainer3d&, const Gascoigne::IntVector&, 
		   const Gascoigne::IntVector&, int, const Gascoigne::IntSet&);
  void ghost_fill_neighbours2d();
  void ghost_fill_neighbours3d();
  void UpdateHangs(HangContainer3d& hangset,
		   const Gascoigne::IntSet& cellref, 
		   const Gascoigne::IntSet& cellcoarse);
  void FaceCoarse(HangContainer3d&, const Gascoigne::IntSet&) const;
  void FaceRefine(HangContainer3d&, const Gascoigne::IntSet&) const;
  void UpdateHangingEdges(HangContainer3d& hangset,
			  const Gascoigne::IntSet& cellref, 
			  const Gascoigne::IntSet& cellcoarse) const;
  void boundary_prepare3d(Gascoigne::IntSet&, Gascoigne::IntSet&, Gascoigne::IntSet&, const HangContainer3d&);
  void new_boundary3d    (Gascoigne::IntSet&, Gascoigne::IntSet&,Gascoigne::IntSet&);
  void new_vertexs3d     (HangContainer3d&, const Gascoigne::IntVector&, const Gascoigne::IntSet&);
  void basic_refine3d    (HangContainer3d&, const Gascoigne::IntSet&, const Gascoigne::IntSet&);
  void basic_fill_neighbours3d();
  void boundary_newton3d      (Gascoigne::IntSet&);
  void inner_vertex_newton3d (const Gascoigne::IntVector&, const Gascoigne::IntSet&, const Gascoigne::IntSet&);
  void update_boundary_data3d(const Gascoigne::IntSet&);
  void new_bquads            (const Gascoigne::IntVector&, const Gascoigne::IntVector&, const Gascoigne::IntSet&);
  void new_middle_vertex3d   (int,int);

  int  regular_grid3d_one  (Gascoigne::IntSet&, Gascoigne::IntVector&, const Gascoigne::IntSet&, const Gascoigne::IntSet&);
  int  regular_grid3d_one  (Gascoigne::IntVector&, Gascoigne::IntVector&, const Gascoigne::IntSet&, const Gascoigne::IntSet&);
  int  regular_grid3d_two  (Gascoigne::IntVector&, const Gascoigne::IntSet&);
  int  regular_grid3d_three_refine(Gascoigne::IntSet&) const;
  int  regular_grid3d_three_coarse(Gascoigne::IntSet&, Gascoigne::IntSet&) const;

  void GetMinMaxLevels(nvector<int>& maxi, nvector<int>& mini,
		       const Gascoigne::IntSet& CellRef) const;

  void init_edges3d();
  void LoadFathers3d(Gascoigne::IntVector& v) const;

  void    _refine3d (Gascoigne::IntSet&, Gascoigne::IntSet&, const Gascoigne::IntVector&, const Gascoigne::IntVector&);
  void FillNeighbourFaces(const Hex& father, const FaceVector& Face,
			  int rneigh);
  void FillNeighbourFaces(int M, int S, const FaceVector& Face);
  void   InitHexOfCurved();
  int   FindPatchDepth() const;
  void  FillVertexLevels(Gascoigne::IntVector& dst) const;
  void  RefineCoarseNodes(Gascoigne::IntSet& dst, const Gascoigne::IntVector& refnodes,
        		  const Gascoigne::IntVector& vertexlevel) const;
  void  VertexToCells(Gascoigne::IntVector& dst, const Gascoigne::IntSet& src, 
		      const Gascoigne::IntVector& vertexlevel) const;
  void VertexToCellsCoarsening(Gascoigne::IntVector& dst, const Gascoigne::IntSet& src, 
			       const Gascoigne::IntVector& vertexlevel) const;
  void recursive_childs(int q, Gascoigne::IntVector& ref, int d) const;

  const CurvedShapes<3>& GetCurvedShapes() const { return _curvedshapes;}
  CurvedShapes<3>& GetCurvedShapes() { return _curvedshapes;}

  public:

  HierarchicalMesh3d();
  HierarchicalMesh3d(const HierarchicalMesh3d& H);
  HierarchicalMesh3d& operator=(const HierarchicalMesh3d& H);
  HierarchicalMesh3d(const Gascoigne::ParamFile* paramfile);
  ~HierarchicalMesh3d() {  GetCurvedShapes().clear();}

  std::string GetName() const {return "HierarchicalMesh3d";}

  /*  Zugriff  */

  int  dimension()            const { return 3;}

  int  nnodes   ()            const { return vertexs3d.size();}
  int  ncells   ()            const { return hexs.size();}
  int  nbquads  ()            const { return Bquads.size();}

  int  nodes_per_cell(int i)  const { return 8;}
  int  VtkType(int i) const { return 12;}

  const Vertex3d& vertex3d(int i)         const { return vertexs3d[i];}

  const Hex&   hex (int i)                const { return hexs[i];}
  const BoundaryQuad&  bquad(int i)       const { return Bquads[i];}

  int  vertex_of_cell (int i, int ii)      const { return hexs[i].vertex(ii);}
  int  vertex_of_bquad(int i, int ii)      const { return Bquads[i].vertex(ii);}
  int  face_of_hex    (int i, int ii)      const { return hexs[i].edge(ii); }
  int  level(int i)                        const { return hexs[i].level();}
  bool sleep(int i)                        const { return hexs[i].sleep();}

  const HexLawAndOrder&     HexLawOrder()     const { return HexLaO; }
  const LineHangList&       linehang()        const { return LineHang;}
  const QuadHangList&       quadhanglist()    const { return QuadHang; }
  const BoundaryFunction3d* quad_shape(int i) const;

  const std::vector<BoundaryQuad>& quad_list() const { return Bquads; }

  const VertexVec3d&        vertex3d() const {return vertexs3d;}
  const HexVec&             hex     () const {return hexs;}
  const BQuadVec&           bquad   () const {return Bquads;}
  const QuadHangList&       quadhang() const {return QuadHang;}
  const std::map<int,int>&       GetHexOfCurved() const {return hexofcurved;}

  /*  Functionen  */

  void   write    (const std::string&) const;
  void   write_gup(const std::string&) const;

  void   WriteAll(const std::string&) const;

  void   write_inp(const std::string&) const;
  void   read_inp (const std::string&);
  void   read_gup (const std::string&);

  void   global_coarse3d();

  void   refine      (const Gascoigne::IntVector&, const Gascoigne::IntVector&);
  void   patch_refine(Gascoigne::IntVector&, Gascoigne::IntVector&);
  //  int    smooth_edges();
  void   FillAllBoundaryLines();

  pint  EdgeNeighbour(int i, int e) const;

  int  NodeOnFace(int e) const;
  fixarray<4,int> ChildrenOfFace(int e) const;

  void  GetVertexesOfFace(fixarray<4,int>&, int) const;
  void  GetVertexesOfFace(fixarray<5,int>&, int) const;
  void GetAwakePatchs(std::set<int>&) const;
  void GetAwakeCells(std::set<int>&) const;
  void ConstructQ2PatchMesh(Gascoigne::IntVector& pm) const;
  std::set<int> GetColors() const;

  
  int nactivedescendants(int i)      const;
  Gascoigne::IntVector GetVertices(int c) const;

  int GetBoundaryCellOfCurved(int iq) const
    {
      std::map<int,int>::const_iterator p = hexofcurved.find(iq);
      if( p!=hexofcurved.end() ) return p->second;
      return -1;
    }
  void Testing();
  int neighbour(int c, int le) const;

  void AddShape(int col, BoundaryFunction<3>* f) {
    GetCurvedShapes().AddShape(col,f);
  }
};

/*---------------------------------------------------*/

#endif
