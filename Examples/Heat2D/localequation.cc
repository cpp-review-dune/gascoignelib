#include  "localequation.h"
#include  "filescanner.h"


using namespace std;

/* ----------------------------------------- */

LocalEquation::LocalEquation(const ParamFile* paramfile) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("visc",&_visc,1.);
  DFH.insert("h",&_h,0.4);
  DFH.insert("k",&_k,2.);
  DFH.insert("r",&_r,0.3);
  FileScanner FS(DFH,paramfile,"Equation");

  _us = (_r*_h)/(1.-_r);
  _vs = (1.-_us)*(_h+_us);


  cerr << "****** us, vs = **** " << _us << " " << _vs << endl; 
}

/* ----------------------------------------- */

void LocalEquation::SetTimePattern(TimePattern& P) const
{
  P.reservesize(ncomp(),ncomp(),0.);
  P(0,0) = 1.;
  P(1,1) = 1.;
}

/* ----------------------------------------- */

void LocalEquation::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  b[0] += _visc* (U[0].x()*N.x()+U[0].y()*N.y());
  b[1] += _visc* (U[1].x()*N.x()+U[1].y()*N.y());

  double u = U[0].m();
  double v = U[1].m();
  double s = u/(u+_h);

  b[0] += N.m() * (-u*(1.-u) + s * v);
  b[1] += N.m() * (_k*_r* v - _k* s * v);
}

/* ----------------------------------------- */

void LocalEquation::Residual(Vector& b, const FemFunction& U, const DerivativeVector& N) const
{
  b[0] += _visc* (U[0].x()*N.x()+U[0].y()*N.y());
  b[1] += _visc* (U[1].x()*N.x()+U[1].y()*N.y());

  double u = U[0].m();
  double v = U[1].m();
  double s = u/(u+_h);

  b[0] += N.m() * (-u*(1.-u) + s * v);
  b[1] += N.m() * (_k*_r* v - _k* s * v);
}

/* ----------------------------------------- */

void LocalEquation::Matrix(EntryMatrix& A, const FemFunction& U, const DerivativeVector& M, const DerivativeVector& N) const
{
  A(0,0) += _visc* (M.x()*N.x()+M.y()*N.y());
  A(1,1) += _visc* (M.x()*N.x()+M.y()*N.y());

  double u = U[0].m();
  double v = U[1].m();
  double MM = M.m()*N.m();

  double s = u/(u+_h);
  double t = _h/((u+_h)*(u+_h));

  A(0,0) += ( (-1.+2.*u) + t*v ) * MM;
  A(0,1) += s * MM;
  A(1,0) += -_k * t*v * MM;
  A(1,1) += ( _k*_r - _k*s ) * MM;
}
