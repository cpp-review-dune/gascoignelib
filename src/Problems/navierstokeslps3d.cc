#include  "navierstokeslps3d.h"
#include  "filescanner.h"

namespace Gascoigne
{

/*-----------------------------------------*/

NavierStokesLps3d::NavierStokesLps3d(const ParamFile* filename) 
  : NavierStokes3d()
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 0.01);
  DFH.insert("cut"  , &_cut  ,  1.e8);
  DFH.insert("alpha", &ST.alpha0, 0.25);
  DFH.insert("delta", &ST.delta0, 0.25);
  DFH.insert("xeta" , &ST.xeta0, 6.);
  DFH.insert("penalty",&_penalty, 0.);

  FileScanner FS(DFH,filename,"Equation");
}

/*-----------------------------------------*/

std::string NavierStokesLps3d::GetName() const
{
  return "NavierStokesLps3d";
}

/*-----------------------------------------*/

void NavierStokesLps3d::init(const nmatrix<double>& H, const FemFunction& U, const Vertex3d& v)const
{
  ST.ReInit(H(0,0),_visc,U[1].m(),U[2].m(),U[3].m());
}

/*-----------------------------------------*/

void NavierStokesLps3d::lpspoint(double h, const FemFunction& U, const Vertex3d& v)const
{
  ST.ReInit(h,_visc,U[1].m(),U[2].m(),U[3].m());
}
 
/*-----------------------------------------*/

void NavierStokesLps3d::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
{
  b[0] += ST.alpha() * Laplace(UP[0],N);
  double betaN  = Convection(U,N);
  double betaU1 = Convection(U,UP[1]);
  double betaU2 = Convection(U,UP[2]);
  double betaU3 = Convection(U,UP[3]);
  b[1] += ST.delta() * betaU1*betaN;
  b[2] += ST.delta() * betaU2*betaN;
  b[3] += ST.delta() * betaU3*betaN;
}
 
/*-----------------------------------------*/

void NavierStokesLps3d::StabilizationResidual(LocalVector& F, const FemFunction& U, const FemFunction& UP, 
				 const FemFunction& N, const FemFunction& NP) const
{
  double betaU1 = Convection(U,UP[1]);
  double betaU2 = Convection(U,UP[2]);
  double betaU3 = Convection(U,UP[3]);
  for (int i=0; i<NP.size(); i++)
    {
      VectorIterator b = F.start(i);
      b[0] += ST.alpha() * Laplace(UP[0],NP[i]);
      double betaN  = Convection(U,NP[i]);
      b[1] += ST.delta() * betaU1*betaN;
      b[2] += ST.delta() * betaU2*betaN;
      b[3] += ST.delta() * betaU3*betaN;
    }
}

/*-----------------------------------------*/

void NavierStokesLps3d::StabilizationMatrix(EntryMatrix& A, const FemFunction& U, const FemFunction& UP, const FemFunction& M, const FemFunction& MP, const FemFunction& N, const FemFunction& NP) const
{
  // MP=NP  !!!
  nvector<double> betaN(NP.size());
  for (int i=0; i<NP.size(); i++)
    {
      betaN[i] = Convection(U,NP[i]);
    }
  for (int j=0; j<NP.size(); j++)
    {
      //double betaM = ST.delta() * Convection(U,MP[j]);
      double betaM = ST.delta() * betaN[j];
      for (int i=0; i<NP.size(); i++)
	{
	  A.SetDofIndex(i,j);
	  double betabeta = betaM*betaN[i];
	  A(0,0) += ST.alpha() * Laplace(NP[j],NP[i]);
	  A(1,1) += betabeta;
	  A(2,2) += betabeta;
	  A(3,3) += betabeta;
	}
    }
}

/*-----------------------------------------*/

void NavierStokesLps3d::StabMatrix(EntryMatrix& A,  const FemFunction& U, 
 const TestFunction& Np, const DerivativeVector& Mp) const
{
  double laplace = Laplace(Mp,Np);
  double betaM = Convection(U,Mp);
  double betaN = Convection(U,Np);
  double betabeta = ST.delta() * betaM*betaN;

  ////////////// Continuity ////////////////////////////////////////////////

  A(0,0) += ST.alpha() * laplace;
  A(1,1) += betabeta;
  A(2,2) += betabeta;
  A(3,3) += betabeta;
}
}

/*-----------------------------------------*/
