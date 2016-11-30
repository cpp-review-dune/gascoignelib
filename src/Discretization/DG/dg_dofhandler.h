#ifndef __dg_dofhandler__h
#define __dg_dofhandler__h


#include "gascoignemesh.h"
#include "matrixinterface.h"
#include "dgedge.h"

namespace Gascoigne
{

  class DofHandlerBase
  {
  protected:

    // base mesh
    const MeshInterface* _MP;

    // Number of total dofs in the discretization
    int _ndofs;
    // Number of dofs per element 
    int _ndofs_per_element;
    // Number of elements in discretization
    int _nelements;

  public:
  DofHandlerBase() : _MP(NULL), _ndofs(-1), _ndofs_per_element(-1), _nelements(-1)
      {}
    virtual ~DofHandlerBase(){}

    ////////// Access
    virtual const MeshInterface* GetMesh() const { assert(_MP); return _MP;}
    
    virtual int ndofs_per_element() const { return _ndofs_per_element; }
    virtual int nelements()         const { return _nelements; } 
    virtual int ndofs()             const { return _ndofs; }

    
    ////////// Initialize
    virtual void ReInit(const MeshInterface* MP) 
    {
      assert(_ndofs_per_element>0);
      _MP=MP;
      _nelements         = MP->ncells();
    }
    // Das muss automatisch aus der FE-Basis kommen
    virtual void SetDofsPerElement(const int dpe)
    {
      _ndofs_per_element = dpe;
    }
    
    //////////
    ////////// Dof-Handling
    //////////

    // init the size of a local vector used for integration 
    virtual void InitLocalVector(LocalVector& U, int ncomp) const
    {
      U.ncomp()=ncomp; U.resize(ndofs_per_element());
    }
    // init the size of a global dof-vector
    virtual void InitGlobalVector(GlobalVector& u, int ncomp) const
    {
      u.ncomp()=ncomp; u.resize(ndofs());
    }
    // get the dof-indices for one element
    virtual void GetIndices(int element, std::vector<int>& indices) const {assert(0);}

    // local to global and global to local functions
    virtual void GlobalToLocal(LocalVector& U, const GlobalVector& u, int el) const
    {
      assert(U.ncomp()==u.ncomp());
      assert(U.n()==this->ndofs_per_element());
      assert(el<this->nelements());
      
      nvector<int> in(this->ndofs_per_element());
      GetIndices(el,in);
      
      for (int i=0;i<this->ndofs_per_element();++i)
	U.equ_node(i, in[i],  u);
    }
    virtual void LocalToGlobal(GlobalVector& u, const LocalVector& U, int el, double d) const
    {
      assert(U.ncomp()==u.ncomp());
      assert(U.n()==this->ndofs_per_element());
      assert(el<this->nelements());
      
      nvector<int> in(this->ndofs_per_element());
      GetIndices(el,in);
      
      for (int i=0;i<this->ndofs_per_element();++i)
	u.add_node(in[i], d, i, U);
    }
    virtual void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int el, double s) const
    {
      std::vector<int> indices(this->ndofs_per_element());
      this->GetIndices(el,indices);
      
      //HN->CondenseHanging(E,indices);
      std::vector<int>::const_iterator  start = indices.begin();
      std::vector<int>::const_iterator  stop  = indices.end();
      A.entry(start,stop,E,s);
    }
    
    
  };


  // Basic dofhandler with vertices
  template<int DIM>
    class DofHandler : virtual public DofHandlerBase
    {
    protected:

      // coordinates for all dofs
      std::vector<Vertex<DIM> > _vertices;
      
    public:

      // Constructure/Destructor
      DofHandler<DIM>() : DofHandlerBase()
	{}
      ~DofHandler<DIM>() {}

      ////////// Access
      // get vertex of the discretization
      const Vertex<DIM>& vertex(int n) const
      {
	assert(n<_vertices.size());
	return _vertices[n];
      }
      // get vertex of the mesh
      const Vertex<DIM>& GetMeshVertex(int n) const
      {
	assert(this->_MP);
	assert(dynamic_cast<const GascoigneMesh<DIM> *> (this->_MP));
	return dynamic_cast<const GascoigneMesh<DIM> *> (this->_MP)->vertex(n);
      }
    };


  // Dof-Handler for DG-Discretizations
  // - dof's per element
  // - element is patch or cell? !!!!!!!!!! CELL-WISE
  // - number of dof's per element. template or parameter?
  // - how to set connectivity?
  //
  template<int DIM>
    class DGDofHandler : virtual public DofHandler<DIM>
    {
    protected:

      // the dg-dof handler has additional info's on the edges of the triangulation
      std::vector<DGEdge> _edges;
      
    public:
      
      DGDofHandler<DIM>() : DofHandler<DIM>() {}
      ~DGDofHandler<DIM>() {}


      
      
      
      ////////// initialize
      virtual void ReInit(const MeshInterface* MP) 
      {
	DofHandler<DIM>::ReInit(MP);

	// DG, set dofs
	this->_ndofs     = this->ndofs_per_element() * this->_nelements;
	      
	this->_vertices.resize(0); // we never need the coordinates of the vertices... I think...
      }

      ////////// Access
      // do not store vertices in DG Q1 - nevertheless, we should not never need them...
      const Vertex<DIM>& vertex(int n) const { abort(); }
      
      int           nedges()    const { return _edges.size(); }
      int  n_dofs_per_edge()    const { return 2 * this->ndofs_per_element(); }
      const DGEdge& edge(int e) const { assert(e<nedges()); return _edges[e]; }

      
      //////////
      ////////// Dof-Handling
      //////////
      void GetIndices(int element, std::vector<int>& indices) const
      {
	assert(indices.size()==this->ndofs_per_element());
	assert(element>=0);
	assert(element<this->nelements());
	for (int q=0;q<this->ndofs_per_element();++q)
	  indices[q]=element*this->ndofs_per_element()+q;
      }
      // get the dof-indices for one edge
      virtual void GetIndicesEdge(int edge, std::vector<int>& indices) const{};
    
      

    };







  // Dof-Handler for CG-Discretizations
  // - dof's per element
  // - element is patch or cell? !!!!!!!!!! CELL-WISE
  // - number of dof's per element. template or parameter?
  // - how to set connectivity?
  //
  template<int DIM>
    class CGDofHandler : virtual public DofHandler<DIM>
    {
    protected:

    public:
      
      CGDofHandler<DIM>() : DofHandler<DIM>() {}
      ~CGDofHandler<DIM>() {}

      
      ////////// initialize
      virtual void ReInit(const MeshInterface* MP) 
      {
	DofHandler<DIM>::ReInit(MP);

	this->_ndofs = this->GetMesh()->nnodes();
	this->_vertices.resize(0); // vertices are not stored in CG-Q1 disc.
      }
      
      
      ////////// Access
      const Vertex<DIM>& vertex(int n) const
      {
	abort();
	return this->GetMeshVertex(n);
      }

  
      //////////
      ////////// dof-Handling
      //////////
      void GetIndices(int element, std::vector<int>& indices) const
      {
	assert(indices.size()==this->ndofs_per_element());
	assert(element>=0);
	assert(element<this->nelements());
	indices = this->GetMesh()->IndicesOfCell(element);
      }
      
    };



}


#endif
