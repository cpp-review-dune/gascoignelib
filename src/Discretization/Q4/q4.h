#ifndef __Q4_h
#define __Q4_h

#include "patchdiscretization.h"
#include "hnstructureq1.h"

namespace Gascoigne
{

/**********************************************************/

  class Q4 : public PatchDiscretization
  {
    protected:
      HNStructureQ1 *HN;

      void GlobalToLocalCell(LocalVector& U, const GlobalVector& u, int iq) const;
      void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;
      void Transformation(FemInterface::Matrix& T, int iq) const;

      nvector<int> GetLocalIndices(int iq) const
      {
        return *GetPatchMesh()->IndicesOfQ4Patch(iq);
      }

    public:
      Q4() : PatchDiscretization(), HN(NULL) { }
      virtual ~Q4() { }

      void ReInit(const MeshInterface* MP);
      void Structure(SparseStructureInterface* S) const;

      void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
      void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
      void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
      void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const;
      void BoundaryMatrix(MatrixInterface& A, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
      void MassMatrix(MatrixInterface& A) const;
      void MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const;

      void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const;

      void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
      void BoundaryRhs(GlobalVector& f, const IntSet& Colors,  const BoundaryRightHandSide& NRHS, double s) const;

      void InitFilter(nvector<double>&) const;

      double ComputeBoundaryFunctional(const GlobalVector& u, const IntSet& Colors, const BoundaryFunctional& BF) const;
      double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const;

      void HNAverage   (GlobalVector& x) const { HN->Average(x); }
      void HNDistribute(GlobalVector& x) const { HN->Distribute(x); }
      void HNZero      (GlobalVector& x) const { HN->Zero(x); }
      bool HNZeroCheck (const GlobalVector& x) const { return HN->ZeroCheck(x); }

      int n() const                { return GetMesh()->nnodes(); }
      int nc() const               { return GetMesh()->ncells(); }
      int n_withouthanging() const { return GetMesh()->nnodes()-HN->nhnodes(); }

      std::string GetName() const {return "Q4";}
  };

/**********************************************************/

}

#endif