#include  "dwrq1q2.h"
#include  "pi.h"
#include  "dwrfem.h"

namespace Gascoigne
{
/*--------------------------------------------------------*/

DwrQ1Q2::DwrQ1Q2(SolverInterface& SR) : S(SR) 
{ 
  primalproblem  = S.GetProblemDescriptor();
  discretization = S.GetDiscretization();
}

/*--------------------------------------------------------*/

double DwrQ1Q2::ScalarProduct(nvector<double>& eta, const GlobalVector& f, 
			     const GlobalVector& z) const
{
  for(int i=0; i<z.n(); i++)
    {
      for (int c=0; c<z.ncomp(); c++)
	{
	  eta[i] += fabs(f(i,c)*z(i,c));
	}
    } 
  return z * f;
}

/*--------------------------------------------------------*/

double DwrQ1Q2::ScalarProduct(nvector<double>& eta, const VectorInterface& gf, 
			     const VectorInterface& gz) const
{
  const GlobalVector& f = S.GetGV(gf);
  const GlobalVector& z = S.GetGV(gz);

  return ScalarProduct(eta,f,z);
}

/*--------------------------------------------------------*/

double DwrQ1Q2::ScalarProductWithFluctuations(nvector<double>& eta, const VectorInterface& gf, 
					     const VectorInterface& gz) const
{
  const GlobalVector& f = S.GetGV(gf);
  const GlobalVector& z = S.GetGV(gz);

  GlobalVector dz(f.ncomp(),f.n());

  dz.zero();
  Pi pi;
  pi.Init(S.GetMesh());
  pi.vmult(dz,z);

  return ScalarProduct(eta,f,dz);
}


/*--------------------------------------------------------*/

DiscretizationInterface* DwrQ1Q2::GetOtherDiscretization() const
{
  DiscretizationInterface* D;

  if (S.GetMesh()->dimension()==2) 
    {
      D = new DwrFemQ1Q22d;    
    }
  else
    {
      D = new DwrFemQ1Q23d;
    }
  
  D->BasicInit(S.GetParamfile());
  return D;
}

/*-------------------------------------------------------*/

void DwrQ1Q2::PrimalResidualsHigher(VectorInterface& gf, const VectorInterface& gu)
{
  GlobalVector& f = S.GetGV(gf);

  f.zero();

  // only necessary if z has additional Dirichlet bc compared to u
  //
  S.Rhs(gf,-0.5);
  S.Form(gf,gu,0.5);

  DiscretizationInterface* D = GetOtherDiscretization();

  S.SetDiscretization(*D,true);
      
  S.Rhs(gf,0.5);
  S.Form(gf,gu,-0.5);
//  S.SetBoundaryVectorZero(gf);
  
  S.SetDiscretization(*discretization);
  delete D;
}

/*--------------------------------------------------------*/

void DwrQ1Q2::DualResidualsHigher(VectorInterface& gf, 
				  const VectorInterface& gu, 
				  const VectorInterface& gz, 
				  const ProblemDescriptorInterface& PDI)
{
  S.GetGV(gf).zero();
  // dual problem
  S.SetProblem(PDI);
  S.AddNodeVector("u",gu);

  // standard residual
  //
  {
    S.Rhs     (gf, -0.5);
    S.AdjointForm(gf,gz,0.5);
//    S.SetBoundaryVectorZero(gf);
    S.HNDistribute(gf);
  }
  // residual respect Q2 test functions
  //
  {  
    DiscretizationInterface* D = GetOtherDiscretization();
    S.SetDiscretization(*D,true);

    S.Rhs     (gf,   0.5);
    S.AdjointForm(gf,gz,-0.5);
//    S.SetBoundaryVectorZero(gf);
    S.HNDistribute(gf);

    S.SetDiscretization(*discretization);
    delete D;
  }

  S.DeleteNodeVector("u");
  S.SetProblem(*primalproblem);
}

/*--------------------------------------------------------*/

double DwrQ1Q2::Estimator(nvector<double>& eta, VectorInterface& gf, 
			  const VectorInterface& gu, const VectorInterface& gz,
			  const ProblemDescriptorInterface& PDI)
{
  double rho=0, rhostern=0;

  DualResidualsHigher(gf,gu,gz,PDI);
  rhostern =  ScalarProductWithFluctuations(eta,gf,gu);
  
  PrimalResidualsHigher(gf,gu);
  rho      =  ScalarProductWithFluctuations(eta,gf,gz);
      
  return rho + rhostern;
}

/*--------------------------------------------------------*/
}
