#ifndef  __IntegratorInterface_h
#define  __IntegratorInterface_h


#include  "gascoigne.h"
#include  "equation.h"
#include  "exactsolution.h"
#include  "feminterface.h"
#include  "entrymatrix.h"
#include  "domainrighthandside.h"
#include  "domainfunctional.h"
#include  "boundaryfunctional.h"
#include  "boundaryrighthandside.h"
#include  "boundaryequation.h"
#include  "diracrighthandside.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments IntegratorInterface

  ///
  ///
  /////////////////////////////////////////////

  class IntegratorInterface
  {
    private:

    protected:

    public:
      IntegratorInterface() {}
      virtual ~IntegratorInterface() {}
  
      virtual std::string GetName() const=0;

      virtual void BasicInit() { assert(0);};
      
      virtual void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, 
          const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::Rhs\" not written!" << std::endl;
						assert(0);
      }
      virtual void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
          const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::Form\" not written!" << std::endl;
						assert(0);
      }
      virtual void AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FEM, 
          const LocalVector& U, const LocalNodeData& Q, const LocalCellData& QC) const { 
        std::cerr << "\"IntegratorInterface::AdjointForm\" not written!" << std::endl;
						assert(0);
      }

      virtual void BoundaryRhs(const BoundaryRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, int ile, 
          int col, const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::BoundaryRhs\" not written!" << std::endl;
						assert(0);
      }
      virtual void BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
          int ile, int col, const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::BoundaryForm\" not written!" << std::endl;
						assert(0);
      }
      virtual void BoundaryMatrix (const BoundaryEquation& BE, EntryMatrix& E, const FemInterface& FEM, 
          const LocalVector& U, int ile, int col, const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::BoundaryMatrix\" not written!" << std::endl;
						assert(0);
      }
      virtual void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
          const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::Matrix\" not written!" << std::endl;
						assert(0);
      }
      virtual double MassMatrix(EntryMatrix& E, const FemInterface& FEM) const {
        std::cerr << "\"IntegratorInterface::MassMatrix\" not written!" << std::endl;
        assert(0);
      }

      virtual void MassForm(const TimePattern& TP, LocalVector& F, const FemInterface& FEM, const LocalVector& U) const {
        std::cerr << "\"IntegratorInterface::MassForm\" not written!" << std::endl;
        abort();
      }

      virtual double ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, 
          const LocalVector& U, const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::ComputeDomainFunctional\" not written!" << std::endl;
						assert(0);
      }
      virtual double ComputeBoundaryFunctional(const BoundaryFunctional& F, const FemInterface& FEM, int ile,
          int col, const LocalVector& U, const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::ComputeBoundaryFunctional\" not written!" << std::endl;
						assert(0);
      }

        
      virtual void EvaluateCellRightHandSide(LocalVector &F, const DomainRightHandSide& CF,const FemInterface& FEM, 
          const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::EvaluateCellRightHandSide\" not written!" << std::endl;
						assert(0); }

      virtual void DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex2d& p, const DiracRightHandSide& DRHS, 
          int i, const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::DiracRhsPoint\" not written!" << std::endl;
						assert(0);
      }
      virtual void DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex3d& p, const DiracRightHandSide& DRHS, 
          int i, const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::DiracRhsPoint\" not written!" << std::endl;
						assert(0);
      }
      virtual double ComputePointValue(const FemInterface& E, const Vertex2d& p, const LocalVector& U, int comp) const {
        std::cerr << "\"IntegratorInterface::ComputePointValue\" not written!" << std::endl;
        assert(0);
      }
      virtual double ComputePointValue(const FemInterface& E, const Vertex3d& p, const LocalVector& U, int comp) const {
        std::cerr << "\"IntegratorInterface::ComputePointValue\" not written!" << std::endl;
        assert(0);
      }
      virtual void ErrorsByExactSolution(LocalVector& dst, const FemInterface& FE, const ExactSolution& ES, 
          const LocalVector& U, const LocalNodeData& Q, const LocalCellData& QC) const {
        std::cerr << "\"IntegratorInterface::ErrorsByExactSolution\" not written!" << std::endl;
        assert(0);
      }
      virtual void IntegrateMassDiag(DoubleVector& F, const FemInterface& FEM) const {
        std::cerr << "\"IntegratorInterface::IntegrateMassDiag\" not written!" << std::endl;
        assert(0);
      }

      virtual void IntegrateBoundaryMassDiag(DoubleVector& F, const FemInterface& FEM, int ile, int col) const {
        std::cerr << "\"IntegratorInterface::IntegrateBoundaryMassDiag\" not written!" << std::endl;
        assert(0);
      }
  };
}

#endif
