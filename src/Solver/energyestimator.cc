#include  "energyestimator.h"
#include  "edgeinfocontainer.h"
#include  "q12d.h"
#include  "q13d.h"
#include  "stdsolver.h"

namespace Gascoigne{

/*--------------------------------------------------------*/

EnergyEstimator::EnergyEstimator(SolverInterface& SR) : S(SR) 
{ 
  primalproblem  = S.GetProblemDescriptor();
  discretization = S.GetMeshInterpretor();
}

/*--------------------------------------------------------*/

double EnergyEstimator::Estimator(nvector<double>& eta, BasicGhostVector& gu, 
				  const BasicGhostVector& gf)
{
  const GlobalVector& u = S.GetGV(gu);

  const Equation*          EQ  = primalproblem->GetEquation();
  const RightHandSideData* RHS = primalproblem->GetRightHandSideData();

  S.HNAverage(gu);
  S.HNAverageData();

  eta.reservesize(u.n());
  eta.zero();
  
  const StdSolver* SS = dynamic_cast<const StdSolver*>(&S);
  assert(SS);

  if (S.GetMesh()->dimension()==2)
    {
      EdgeInfoContainer<2> EIC;
      EIC.BasicInit(SS->GetHierarchicalMesh(),u.ncomp());

      const Q12d* DP = dynamic_cast<const Q12d*>(discretization);
      assert(DP);
      if(RHS)
      {
        const DomainRightHandSide *DRHS = dynamic_cast<const DomainRightHandSide *>(RHS);
        assert(DRHS);
        DP->EnergyEstimator(EIC,eta,u,*EQ,*DRHS);
      }
      else
      {
        DP->EnergyEstimatorZeroRhs(EIC,eta,u,*EQ);
      }
    }
  else if (S.GetMesh()->dimension()==3)
    {
      EdgeInfoContainer<3> EIC;
      EIC.BasicInit(SS->GetHierarchicalMesh(),u.ncomp());
      const Q13d* DP = dynamic_cast<const Q13d*>(discretization);
      assert(DP);
      if(RHS)
      {
        const DomainRightHandSide *DRHS = dynamic_cast<const DomainRightHandSide *>(RHS);
        assert(DRHS);
        DP->EnergyEstimator(EIC,eta,u,*EQ,*DRHS);
      }
      else
      {
        DP->EnergyEstimatorZeroRhs(EIC,eta,u,*EQ);
      }
    }
  S.HNZero(gu);
  S.HNZeroData();

  return eta.norm();
}

}
/*--------------------------------------------------------*/


