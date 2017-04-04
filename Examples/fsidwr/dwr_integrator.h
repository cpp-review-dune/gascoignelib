/*----------------------------   dwr_integrator.h     ---------------------------*/
/*      $Id: dwr_integrator.h,v 1.2 2009/07/08 20:59:17 richter Exp $                 */
#ifndef __dwr_integrator_H
#define __dwr_integrator_H
/*----------------------------   dwr_integrator.h     ---------------------------*/


#include "basicintegrator.h"
#include "domainrighthandside.h"

namespace Gascoigne
{

  /*---------------------------------------------------*/
  
  template<int DIM>
    class DWR_Integrator : public BasicIntegrator
    {
    protected:
      double Volume2MeshSize(double vol) const { return pow(vol,1./float(DIM));}
      int PatchMeshNr2IntegratorNr(int in) const;
      int PatchMeshNr2IntegratorNrBoundary(int in, int ile) const;
    public:
      
      DWR_Integrator<DIM>() : BasicIntegrator() {}
      ~DWR_Integrator<DIM>() {}
      
      std::string GetName() const {return "DWR Integrator";}
      void BasicInit();
      
      void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FemAnsatz, const FemInterface& FemTrial, const LocalData& Q, const LocalData& QC) const; 
      void Form(const Equation& EQ, LocalVector& F, const FemInterface& FemH, 
		const FemInterface& FemL, const LocalVector& U, const LocalData& Q, const LocalData& QC) const;
      void AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FemH, 
		       const FemInterface& FemL, const LocalVector& U, const LocalData& Q, const LocalData& QC) const;
      void BoundaryRhs(const BoundaryRightHandSide& RHS, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, 
		       int ile, int col, const LocalData& Q, const LocalData& QC) const;
      void BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, 
			const LocalVector& U, int ile, int col, LocalData& Q, const LocalData& QC) const;
      void DiracRhsPoint(LocalVector& b, const FemInterface& FemH, const FemInterface& FemL, const Vertex<DIM>& p, 
			 const DiracRightHandSide& DRHS, int j, const LocalData& Q, const LocalData& QC) const;
      
      double MassMatrix(EntryMatrix& E, const FemInterface& FemH, const FemInterface& FemL) const;
      void MassForm(const TimePattern& TP, LocalVector& F, const FemInterface& FemH, const FemInterface& FemL, const LocalVector& U) const;
    };

}




/*----------------------------   dwr_integrator.h     ---------------------------*/
/* end of #ifndef __dwr_integrator_H */
#endif
/*----------------------------   dwr_integrator.h     ---------------------------*/
