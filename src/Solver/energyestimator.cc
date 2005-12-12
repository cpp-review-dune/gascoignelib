#include  "energyestimator.h"
#include  "edgeinfocontainer.h"
#include  "q12d.h"
#include  "q13d.h"
#include  "stdsolver.h"

namespace Gascoigne{

/*--------------------------------------------------------*/

EnergyEstimator::EnergyEstimator(SolverInterface& SR) : S(SR) 
{ 
  {
    DataFormatHandler DFH; 
    DFH.insert("estimator", &_s_energytype ,"energy_laplace"); 
    FileScanner FS(DFH); 
    FS.NoComplain(); 
    FS.readfile(SR.GetParamfile(),"Loop"); 
  }
  {
    DataFormatHandler DFH; 
    DFH.insert("visc", &_d_visc ,1); 
    FileScanner FS(DFH); 
    FS.NoComplain(); 
    FS.readfile(SR.GetParamfile(),"Equation"); 
  }
  primalproblem  = S.GetProblemDescriptor();
  discretization = dynamic_cast<Q1*>(S.GetDiscretization());
  assert(discretization);
}

/*--------------------------------------------------------*/

double EnergyEstimator::Estimator(DoubleVector& eta, VectorInterface& gu, 
				  const VectorInterface& gf)
{
  const GlobalVector& u = S.GetGV(gu);

  const Equation*    EQ  = primalproblem->GetEquation();
  const Application* RHS = primalproblem->GetRightHandSide();

  S.HNAverage(gu);
  S.HNAverageData();

  eta.reservesize(u.n());
  eta.zero();
  
  const StdSolver* SS = dynamic_cast<const StdSolver*>(&S);
  assert(SS);

  const DomainRightHandSide* DRHS = NULL;
  if(RHS)
    {
      DRHS = dynamic_cast<const DomainRightHandSide *>(RHS);
      assert(DRHS);
    }
  EdgeInfoContainerInterface* EIC = NULL;
  if      (S.GetMesh()->dimension()==2)  EIC = new EdgeInfoContainer<2>;
  else if (S.GetMesh()->dimension()==3)  EIC = new EdgeInfoContainer<3>;

  EIC->BasicInit(SS->GetHierarchicalMesh(),u.ncomp());
  discretization->EnergyEstimator(*EIC,eta,u,*EQ,DRHS,_s_energytype,_d_visc);
  delete EIC;

  S.HNZero(gu);
  S.HNZeroData();

  return eta.norm();
}

}
/*--------------------------------------------------------*/


