#ifndef __dg_disc__h
#define __dg_disc__h

#include "dg_dofhandler.h"
#include "discretizationinterface.h"
#include "galerkinintegrator.h"
#include "baseq1.h"
#include  "finiteelement.h"
#include  "transformation.h"
#include  "problemdescriptorinterface.h"
#include  <string>
#include  "edgeintegrator.h"

namespace Gascoigne
{

  
  template<int DIM>
  class Disc : public DiscretizationInterface
  {
  protected:

    DofHandlerBase*       _dof_handler;
    GalerkinIntegrator<DIM> _integrator;

    typedef Gascoigne::Transformation<DIM, BaseQ1<DIM> >     TransType;
    typedef FiniteElement<DIM,DIM-1,TransType,BaseQ1<DIM> >  FiniteElementType;

  public:

    Disc<DIM> () : _dof_handler(NULL)
      {
	_integrator.BasicInit();
      }
    ~Disc<DIM>()
      {
	if (_dof_handler!=0) delete _dof_handler;
	_dof_handler=NULL;
      }
    
    
    void ReInit   (const MeshInterface* M)
    {
      assert(_dof_handler);
      _dof_handler->ReInit(M);
    }
    
    
    void SetDofHandler(DofHandlerBase* DH)  { _dof_handler = DH;  abort(); }
    
    ////////// Access
    const DofHandlerBase* GetDofHandler() const { return _dof_handler; }
    const GascoigneMesh<DIM>* GetMesh() const
    {
      assert(_dof_handler);
      assert(dynamic_cast<const GascoigneMesh<DIM>* >(_dof_handler->GetMesh()));
      return dynamic_cast<const GascoigneMesh<DIM>* >(_dof_handler->GetMesh());
    }

    int ndofs() const
    { assert(_dof_handler); return GetDofHandler()->ndofs(); }
    int nelements() const
    { assert(_dof_handler); return GetDofHandler()->nelements(); }


    ////////// Handling of Dofs and Mesh
    void Transformation(FemInterface::Matrix& T, int iq) const;
    void InitColoring() 
    { std::cerr << "Disc:: InitColering not written!" << std::endl; }
    void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT)
    { std::cerr << "\"Disc::ConstructInterpolator \" not written!" << std::endl; }
    void InitFilter(DoubleVector&) const
    { std::cerr << "\"Disc::InitFilter\" not written!" << std::endl; }

    
    void Structure(SparseStructureInterface* S) const;

    // Integration-Stuff
    void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
    void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
    void Matrix(MatrixInterface& A, const GlobalVector& u, const ProblemDescriptorInterface* PD, double) const;

    // I/O
    void WriteVtk(std::string name, const GlobalVector& u) const;
  };



  template<int DIM>
    class DGDisc : public Disc<DIM>
    {
    protected:
      
      typedef Gascoigne::Transformation<DIM, BaseQ1<DIM> >     TransType;
      typedef FiniteElement<DIM,DIM-1,TransType,BaseQ1<DIM> >  FiniteElementType;

      EdgeIntegrator<DIM> _edge_integrator;

    public:
      
      DGDisc<DIM> () : Disc<DIM>()
      {
      }
      ~DGDisc<DIM>()
	{
	}
      
      ////////////////////////////////////////////////// Init
      void BasicInit(const ParamFile* pf)
      {
	assert(this->_dof_handler==NULL);
	this->_dof_handler = new DGDofHandler<DIM>();
	this->_dof_handler->SetDofsPerElement(4);
      }

      //////////////////////// Access
      const DGDofHandler<DIM>* GetDGDofHandler() const
      {
	assert(dynamic_cast<const DGDofHandler<DIM>* >(this->_dof_handler));
	return dynamic_cast<const DGDofHandler<DIM>* >(this->_dof_handler);
      }

      //////////////// INtegration DG-specific
      void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
      
    // I/O
      void WriteVtk(std::string name, const GlobalVector& u) const;
    };
  




  template<int DIM>
    class CGDisc : public Disc<DIM>
  {
  protected:

  public:

    CGDisc<DIM> () : Disc<DIM>()
      {
      }
    ~CGDisc<DIM>()
      {
      }

    ////////////////////////////////////////////////// Init
    void BasicInit(const ParamFile* pf)
    {
      assert(this->_dof_handler==NULL);
      this->_dof_handler = new CGDofHandler<DIM>();
      this->_dof_handler->SetDofsPerElement(4);
    }
  };



}


#endif


