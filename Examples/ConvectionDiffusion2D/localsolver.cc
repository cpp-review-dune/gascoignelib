#include  "localsolver.h"
#include  "localmeshinterpretor.h"


using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

void LocalSolver::NewMeshInterpretor(int dimension, const string& discname)
{
  GetMeshInterpretorPointer()  = new LocalMeshInterpretor;
}

/*-------------------------------------------------------*/

void LocalSolver::InterpolateSolution(GlobalVector& u, const GlobalVector& uold) const
{
  u.zero();
  GetMeshInterpretor()->InterpolateSolution(u, uold);
  // geht nur mit ghostvector!
  //  PressureFilterIntegrate(u);
}

/*-----------------------------------------*/

void LocalSolver::HNAverage(GlobalVector& x) const
{
  GetMeshInterpretor()->HNAverage(x);
}

/*-------------------------------------------------------*/

void LocalSolver::HNZero(GlobalVector& x) const
{
  GetMeshInterpretor()->HNZero(x);
}

/*-------------------------------------------------------*/

void LocalSolver::AddGlobalData(const GlobalVector& d)
{
  GetMeshInterpretor()->AddGlobalData(&d);
}
