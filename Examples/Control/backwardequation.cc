#include  "backwardequation.h"
#include  "filescanner.h"


using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

BackwardEquation::BackwardEquation(const ParamFile* paramfile) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("visc",&_visc,1.);
  FileScanner FS(DFH, paramfile, "Backward");
}

/* ----------------------------------------- */

void BackwardEquation::SetTimePattern(TimePattern& P) const
{
  P.reservesize(GetNcomp(),GetNcomp(),0.);
  for(int c=0;c<GetNcomp();c++)
  { 
    P(c,c) = 1.;
  }
}

/* ----------------------------------------- */

void BackwardEquation::SetFemData(FemData& Q) const
{
  assert(Q.count("U"));
  q = &Q["U"];
}

/* ----------------------------------------- */

void BackwardEquation::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  b[0] += _visc* (U[0].x()*N.x()+U[0].y()*N.y());
  b[0] += N.m() * 2.*(*q)[0].m()*U[0].m();
}

/* ----------------------------------------- */

void BackwardEquation::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
{
  A(0,0) += _visc* (M.x()*N.x()+M.y()*N.y());
  A(0,0) += 2.*(*q)[0].m()*M.m()*N.m();
}
