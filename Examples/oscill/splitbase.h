/*----------------------------   splitbase.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __splitbase_H
#define __splitbase_H
/*----------------------------   splitbase.h     ---------------------------*/


/*-------------------------------------------------------------*/

#include "meshinterface.h"

#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#define HASHMAP   std::tr1::unordered_map
#define HASHSET   std::tr1::unordered_set
#else
#include  <ext/hash_map>
#include  <ext/hash_set>
#define HASHMAP  __gnu_cxx::hash_map
#define HASHSET  __gnu_cxx::hash_set
#endif

namespace Gascoigne
{

  
  class SplitBase
  {
  public:

    virtual const MatrixInterface* GetFluid_NS()  const {abort();}
    virtual const MatrixInterface* GetFluid_EXT() const {abort();}
    virtual const MatrixInterface* GetFluid_ALE() const {abort();}
    virtual const MatrixInterface* GetSolid11()     const {abort();}
    virtual const MatrixInterface* GetSolid12()     const {abort();}
    virtual const MatrixInterface* GetSolid21()     const {abort();}
    virtual const MatrixInterface* GetSolid22()     const {abort();}
    

    const MeshInterface* __mesh;
    
    HASHSET<int>    * __interface_nodes;
    HASHMAP<int,int>* __fluid_g2l_map; 
    HASHMAP<int,int>* __solid_g2l_map; 
    std::vector<int>* __fluid_l2g; 
    std::vector<int>* __solid_l2g; 

    std::vector<int> __fluid_g2l; 
    std::vector<int> __solid_g2l; 
    
    HASHSET<int> __fluid_nodes,__solid_nodes;

    int __nnodes;
    
    
  public:
    SplitBase()
      { __mesh = 0; __nnodes = 0; __interface_nodes=0;
	__fluid_g2l_map=0; __solid_g2l_map=0;
	__fluid_l2g=0; __solid_l2g=0; }

    void SetInterface(int nnodes,
		      HASHSET<int>    & interface_nodes,
		      std::vector<int>& fluid_l2g, 
		      std::vector<int>& solid_l2g,
		      HASHMAP<int,int>& fluid_g2l,
		      HASHMAP<int,int>& solid_g2l,
		      const MeshInterface* mesh)
    {
      __mesh = mesh;
      
      __nnodes = nnodes;
      __interface_nodes = &interface_nodes;
      __fluid_g2l_map = &fluid_g2l;
      __solid_g2l_map = &solid_g2l;
      __fluid_l2g = &fluid_l2g;
      __solid_l2g = &solid_l2g;

      __fluid_g2l.clear(); __fluid_g2l.resize(nnodes,-1);
      __solid_g2l.clear(); __solid_g2l.resize(nnodes,-1);
      
      for (HASHMAP<int,int>::const_iterator it = fluid_g2l.begin();
	   it!=fluid_g2l.end();++it) __fluid_g2l[it->first] = it->second;
      for (HASHMAP<int,int>::const_iterator it = solid_g2l.begin();
	   it!=solid_g2l.end();++it) __solid_g2l[it->first] = it->second;
    }

    int n_fluid_nodes() const  { return __fluid_l2g->size(); }
    int n_solid_nodes() const  { return __solid_l2g->size(); }
    
    
    bool fluid_node(int n) const { return __fluid_g2l[n]!=-1; }
    bool solid_node(int n) const { return __solid_g2l[n]!=-1; }
    bool inter_node(int n) const
    { return (__interface_nodes->find(n)!=__interface_nodes->end()); }

    int Fg2l(int i) const {assert(i<__nnodes); return __fluid_g2l[i]; }
    int Sg2l(int i) const {assert(i<__nnodes); return __solid_g2l[i]; }
    int Fl2g(int i) const {assert(i<__fluid_l2g->size());  return (*__fluid_l2g)[i]; }
    int Sl2g(int i) const {assert(i<__solid_l2g->size());  return (*__solid_l2g)[i]; }
    


    
  };
}




/*----------------------------   splitbase.h     ---------------------------*/
/* end of #ifndef __splitbase_H */
#endif
/*----------------------------   splitbase.h     ---------------------------*/
