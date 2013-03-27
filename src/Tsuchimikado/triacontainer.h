/*----------------------------   triacontainer.h     ---------------------------*/
#ifndef __triacontainer_H
#define __triacontainer_H
/*----------------------------   triacontainer.h     ---------------------------*/

#include "tsuchimikado.h"
#include "element.h"
#include "meshvertex.h"
#include <tr1/unordered_set>
#include <tr1/unordered_map>

#include <map>
#include <set>
#include <list>
#include <utility>
#include "boundaryfunction.h"




namespace Tsuchimikado
{
    
  template<int DIM>
    class TriaContainer
    {

      
    public:

      /**
       * Elements of the mesh
       **/
      typedef Element<3>       HEX;
      typedef Element<2>       QUAD;
      typedef Element<1>       LINE;
      typedef MeshVertex<DIM>  VERTEX;
      typedef Element<DIM>     CELL;

      //      typedef __gnu_cxx::hash_map<int,int> HASH_MAP;
      typedef std::tr1::unordered_map<int,int> HASH_MAP;
      
    protected:

      ////////////////////////////////////////////////// DATA

      unsigned int __mesh_flags;
      
      /**
       * Data of the Mesh
       **/
      std::vector<VERTEX> __vertices;
      std::vector<HEX>    __hexes;
      std::vector<QUAD>   __quads;
      std::vector<LINE>   __lines;

      /**
       * Hanging nodes
       *    map node -> line
       *    map node -> quad
       **/
      HASH_MAP    __hanging_lines;
      HASH_MAP    __hanging_quads;
  
      /**
       * Boundary
       *    map line -> color
       *    map quad -> color
       **/
      HASH_MAP __boundary_lines;
      HASH_MAP __boundary_quads;

      /**
       * Curved boundary
       **/
      std::map<int,const BoundaryFunction<DIM>* > __curved_boundary;
      

      ////////////////////////////////////////////////// INTERNAL FUNCTIONS

      
      /**
       * reinits the structure for hanging nodes
       * has to be called after refining the mesh
       **/
      void reinit_hanging_nodes();

      /**
       * adjust the curved boundaries
       **/
      void reinit_curved();

      /**
       * adds new element to the mesh and return id
       **/
      int add(const VERTEX& v);
      int add(const HEX&  h);
      int add(const QUAD& q);
      int add(const LINE& l);

      /**
       * add lines and quads to the boundary with given color
       **/
      void AddBoundaryLine(int l,int color);
      void AddBoundaryQuad(int l,int color);
      
      /**
       * set master/slave info of line 'li'. if master unset set it to id
       * otherwise slave has to be unset, and will be set to id.
       **/
      void set_adjacent_quad(int li,int id);

      /**
       * set master/slave info of line 'li'. if master unset set it to id
       * otherwise slave has to be unset, and will be set to id.
       * if both are set, replace the one pointing to old.
       **/
      void set_adjacent_quad_and_replace(int li,int id, int old);
      
      /**
       * set master/slave info of quad 'qi'.
       *   if slave/master was _old, set it to _new
       *   otherwise slave/mast must be -1, set it to _new
       **/
      void set_adjacent_hex_and_replace(int qi,int _old,int _new);


      /**
       * check if flagging is ok.
       **/
      int resolve_flags();


      /**
       * refines a given element by recursively refining the boundary first
       * can only be called, if the flagging in the mesh is ok.
       * 
       * Standard argument for refinement type is isotropic.
       **/
      bool  refine_line(const int i,int type = 1);
      bool  refine_quad(const int i,int type = 3);
      bool  refine_hex (const int i,int type = 7);
      bool  refine_cell(const int i,int type = (DIM*(DIM-1)+1));

      int refine_boundary_quad(int hn, int qn, int type = 3);
      
      
      /**
       * splits on hex anisotropically in one direction
       **/
      int split_hex_one_way(const int hn, const mesharray<27,int>& cn, int type);
      int split_hex_isotropic(const int hn, const mesharray<27,int>& cn, int type);

      /**
       * sets the boundary quad information for a hex
       **/
      void set_boundary_quads(int hn,int q0,int q1,int q2,int q3,int q4,int q5);
      /**
       * tries to set refinement flag on elements.
       * returns 1 if flag could be set on element and on boundary of element
       **/
      int set_refine_flag_on_line         (int c, int type = 1);
      int set_refine_flag_on_quad         (int c, int type = 3);
      int set_refine_flag_on_hex          (int c, int type = 7);

      /**
       * sets the refine flag on the boundary quad of a hex h,
       * the type is as seen from the hex in standard orientation
       **/
      int set_refine_flag_on_boundary_quad(int h, int bq, int type);

      /**
       * creates a new hex given by the 8 nodes indicated by "nodes"
       * and indiced by cnode. returns the number of the new element
       * the father is set to hex_number
       **/
      int create_child_hex(const int hex_number,const mesharray<27,int>& cnode,
			   const mesharray<8,int>& nodes,
			   const mesharray<6,int>& quads);  




      // Find element ID in the mesh.
      const int find(const HEX&  H) const;
      const int find(const QUAD& H) const;
      const int find(const LINE& H) const;


      // find local line-index 'i' of quad Q, such that Q.line(i)== line
      const int find_local_line_index_of_quad(int quad_number, int line_number) const;
      
      // middle lines of isotropically refined quads
      const LINE& middle_line_of_quad_at_node(int quad_number,int node) const;

      /**
       *
       * returns the index of a line of quad qn with
       * nodes n1 and n2
       *
       **/
      const int line_of_quad_at_nodes(int qn,int n1,int n2) const;  
      const int child_of_line_at_node(int li,int ni) const;
      const int child_of_quad_at_node(int qi,int ni) const;
      const int child_of_quad_at_line(int qi,int li) const;

      /**
       * auxiliary functions for i/o
       **/
      void read_data(HASH_MAP& v,std::istream& s) const
      {
	size_t n;
	s >> n;
	for (int i=0;i<n;++i)
	  {
	    int fi,se;
	    s >> fi >> se;
	    v.insert(std::make_pair<int,int>(fi,se));
	  }
      }
      void write_data(const HASH_MAP& v,std::ostream& s) const
      {
	s << v.size() << std::endl;
	for (typename HASH_MAP::const_iterator it = v.begin();it!=v.end();++it)
	  {
	    s << it->first << std::endl;
	    s << it->second << std::endl;
	  }
      }

      
      /**
       * Check if mesh is ok (for debugging)
       *  If something is not ok it's too late
       **/
      void sanity_test_lines() const;
      void sanity_test_quads() const;
      void sanity_test_hexes() const;
      void sanity_test_hanging_nodes() const;

      
    public:

      TriaContainer(unsigned int flags=MESH_IS_ISOTROPIC);
      void SetAnisotropic()
      {
	__mesh_flags &= (!MESH_IS_ISOTROPIC);
	assert(!isotropic());
      }
      

      //////////////////////////////////////// PUBLIC FUNCTIONS
      bool isotropic() const;
      /**
       * clears the mesh data
       **/
      void reset();

      /**
       * Add Curved Boundary function
       **/
      void AddCurvedShape(int col, const BoundaryFunction<DIM>* BF);
      const std::map<int,const BoundaryFunction<DIM>* >& GetCurvedShapes() const;      
      
  
      /**
       * access to data of mesh
       **/
      const int dimension() const { return DIM; }
      const int nvertices() const { return __vertices.size(); }
      const int nhexes()    const { return __hexes.size(); }
      const int nquads()    const { return __quads.size(); }
      const int nlines()    const { return __lines.size(); }
      const int ncells()    const;
      
      const VERTEX& vertex(int i) const { assert(i>=0); assert(i<nvertices()); return __vertices[i]; }
      const HEX&      hex (int i) const { assert(i>=0); assert(i<nhexes());    return __hexes[i]; }
      const QUAD&     quad(int i) const { assert(i>=0); assert(i<nquads());    return __quads[i]; }
      const LINE&     line(int i) const { assert(i>=0); assert(i<nlines());    return __lines[i]; }
      const CELL&     cell(int i) const;

      VERTEX& vertex(int i) { assert(i>=0); assert(i<nvertices()); return __vertices[i]; }
      HEX&      hex (int i) { assert(i>=0); assert(i<nhexes());    return __hexes[i]; }
      QUAD&     quad(int i) { assert(i>=0); assert(i<nquads());    return __quads[i]; }
      LINE&     line(int i) { assert(i>=0); assert(i<nlines());    return __lines[i]; }
      CELL&     cell(int i);
      
      const HASH_MAP& GetBoundaryLines() const { return __boundary_lines; }
      const HASH_MAP& GetBoundaryQuads() const { return __boundary_quads; }

      const HASH_MAP& GetHangingLines()  const { return __hanging_lines; }
      const HASH_MAP& GetHangingQuads()  const { return __hanging_quads; }

      /**
       * Get Information about the Mesh and elements of the mesh
       **/
      /**
       * returns index of line number ni of hex hn as
       * indicated in setup.h
       **/
      const int line_of_hex(int hn,int ni) const;
      /**
       * middle nodes of hexes, quads and lines
       **/
      const int middle_node(const HEX&  H) const;
      const int middle_node(const QUAD& Q) const;
      const int middle_node(const LINE& L) const;

      bool line_at_boundary(int l) const
      {
	assert(DIM>=2);
	return (__boundary_lines.find(l)!=__boundary_lines.end());
      }
      bool quad_at_boundary(int l) const
      {
	assert(DIM==3);
	return (__boundary_quads.find(l)!=__boundary_quads.end());
      }
      int color_of_boundary_line(int l) const
      {
	assert(line_at_boundary(l));
	return __boundary_lines.find(l)->second;
      }
      int color_of_boundary_quad(int l) const
      {
	assert(quad_at_boundary(l));
	return __boundary_quads.find(l)->second;
      }
      bool is_hanging(const LINE& L) const;
      bool is_hanging(const QUAD& Q) const;
      
      bool node_is_hanging(int node) const
      {
	if (DIM==2)
	  return (__hanging_lines.find(node)!=__hanging_lines.end());
	if (DIM==3)
	  {
	    bool hl =  (__hanging_lines.find(node)!=__hanging_lines.end());
	    if (hl) return true;
	    else return (__hanging_quads.find(node)!=__hanging_quads.end()); 
	  }
	abort();
      }

  
      /**
       *
       * returns rotation of boundary quad with respect
       * to the standard orientation.
       * The first return value specifies the place of quad.node(0)
       * in the standard numeration of quad.node(3) if in reverse direction
       * The second value indicates the direction of the numeration       
       *
       * i.e., quad 1: standard: 1 5 6 2
       *
       *  (0,1)    1 5 6 2   (0,-1) 2 6 5 1
       *  (1,1)    5 6 2 1   (1,-1) 1 2 6 5
       *
       * This function is necessary to apply the correct flagging
       * to the boundary of an element.
       * e.g. if hex is flagges for type REF_X, quads 0 2 3 4 have
       * to be flagged. Depending on the rotation, different flagging
       * has to be applied to this quads.
       **/
      std::pair<int,int> rotation_of_boundary_quad(const int hex_number,
						   const int boundary_quad_number) const;

      /**
       * returns the location of node 0 of boundary quad
       * in the standard ordering.
       * i.e. quad 1, standard H.node(1) H.node(5) H.node(6) H.node(2)
       * Q.node(0) = H.node(1) return 0
       * Q.node(0) = H.node(6) return 2
       * ...
       **/
      int node_0_of_boundary_quad(const int hex_number,const int boundary_quad_number) const;
  
      // refine, new vertices
      MeshVertex<DIM> new_middle_vertex(const LINE& L) const;
      MeshVertex<DIM> new_middle_vertex(const QUAD& Q) const;
      MeshVertex<DIM> new_middle_vertex (const HEX& H)  const;
  
      //////////////////////////////////////// PUBLIC FUNCTIONS FOR REFINEMENT
      /**
       * clears all the refine flags in mesh
       **/
      void clear_refine_flags();
      /**
       * Flag element c in direction type.
       * returns 1 if setting of new flag is possible
       **/
      int set_refine_flag        (int c, int type = (DIM*(DIM-1)+1));
      /**
       * preforms a global refinement of the mesh
       **/
      void global_refine();
      /**
       * preforms 'i' global refinements of the mesh
       **/
      void global_refine(int i){ for (int j=0;j<i;++j) global_refine(); }
      /**
       * refines all flagged element
       * first adds refinement flags to garuantee a valid mesh,
       * then refines everything.
       **/
      int  refine_cells();
      /**
       * takes care of jobs to be done after refinenent
       *  - curved boundaries
       **/
      void post_refine();

  
      //////////////////////////////////////// FUNCTIONS FOR I/O

      /**
       * simple gnuplot output of all active lines in mesh
       **/
      void print_gnuplot(std::string fn,int extra=0) const;
      
      /**
       * reads inp-file
       **/
      void read_inp(const std::string& fn);

      /**
       * write all mesh data to file
       **/
      void write_gup(const std::string& fn) const;
      void write_gip(const std::string& fn) const;
      /**
       * read all mesh data from file
       **/
      void read_gup (const std::string& fn);
      void read_gip (const std::string& fn);
  


      //////////////////////////////////////// OTHER STUFF
      
      /**
       * Check if mesh is ok (for debugging)
       *  If something is not ok it's too late
       **/
      void sanity_test() const;
      
  
  
    };

}

/*----------------------------   triacontainer.h     ---------------------------*/
/* end of #ifndef __triacontainer_H */
#endif
/*----------------------------   triacontainer.h     ---------------------------*/
