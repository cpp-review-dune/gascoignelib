/*----------------------------   meshlevel.h     ---------------------------*/
/*      $Id: meshlevel.h,v 1.7 2009/11/08 21:22:50 richter Exp $                 */
#ifndef __meshlevel_H
#define __meshlevel_H
/*----------------------------   meshlevel.h     ---------------------------*/

#include "triacontainer.h"

namespace Tsuchimikado
{

  /**
   * vector of elements from a triacontainer to be one mesh level
   * this always has to be a valid mesh:
   *   no hanging node on hanging node
   *
   * public functions:
   *   MeshLevel<DIM>(TC): sets reference to TriaContainer
   *   init_active(): initializes Mesh-level as Mesh fine-mesh of TC
   *   init_coarse(ML): initializes Mesh-Level as coarsening of ML.
   *
   * protected:
   *   init_hanging(): called after the inits
   *   mesh_ok(), repair_mesh(): look and resolve problems in mesh
   **/
  
  template<int DIM> 
    class MeshLevel : public std::vector<int>
    {
    protected:
      //      typedef HASH_SET     HASH_SET;
      typedef std::tr1::unordered_map<int,int> HASH_MAP;
      typedef std::tr1::unordered_set<int>     HASH_SET;
      //      typedef __gnu_cxx::hash_set<int>     HASH_SET;
      //      typedef __gnu_cxx::hash_map<int,int> HASH_MAP;
      
      typedef std::vector<int>::const_iterator IT;
      typedef Element<DIM>             CELL;
      typedef Element<1>               LINE;
      typedef Element<2>               QUAD;
      
      ////////////////////////////////////////////////// DATA
      const TriaContainer<DIM>* __TC;

      HASH_MAP __Cg2l;
      
      HASH_SET     __hanging_lines;
      HASH_SET     __hanging_quads;
      HASH_MAP __hanging_nodes_on_lines;
      HASH_MAP __hanging_nodes_on_quads;

      // nodes that hang but are needed to resolve another hanging node
      // used for generating a valid mesh. This set has to be empty for
      // a valid mesh.
      HASH_MAP     __double_hanging_nodes;
      
      HASH_SET     __changed;

      void repair_mesh();
      bool mesh_ok();
      void init_hanging();

      
    public:
      ~MeshLevel() {}
      MeshLevel(const TriaContainer<DIM>& TC) : __TC(&TC) {}
      MeshLevel() : __TC(0) { abort(); }
      
      int Cellg2l(int g) const
	{
	  assert(__Cg2l.find(g)!=__Cg2l.end());
	  return __Cg2l.find(g)->second;
	}
      int Celll2g(int g) const
	{
	  assert(g<this->size());
	  return (*this)[g];
	}
      bool CellInMesh(int g) const
	{
	  return(__Cg2l.find(g)!=__Cg2l.end());
	}

      ////////////////////////////////////////////////// Hanging
      bool line_is_hanging(int l) const
	{ return __hanging_lines.find(l)!=__hanging_lines.end(); }
      bool quad_is_hanging(int q) const
	{ return __hanging_quads.find(q)!=__hanging_quads.end(); }

      const HASH_SET& GetHangingLines() const { return __hanging_lines; }      
      const HASH_SET& GetHangingQuads() const { return __hanging_quads; }      

      
      // init with all active cells
      void init_active();
      // init by coarsening
      void init_from_meshlevel(const MeshLevel<DIM>& ML, const HASH_SET& dont_coarse);
      void init_from_meshlevel(const MeshLevel<DIM>& ML)
      {
	HASH_SET nix;
	init_from_meshlevel(ML,nix);
      }

      void print_gnuplot(const std::string& fname) const;
      
      //just for testing .... gnuplot with lines
      void print_gnuplot2(std::string fn) const;


    };
}



/*----------------------------   meshlevel.h     ---------------------------*/
/* end of #ifndef __meshlevel_H */
#endif
/*----------------------------   meshlevel.h     ---------------------------*/
