#include "convectiondiffusionequation2d.h"
#include "filescanner.h"

namespace Gascoigne{

/*---------------------------------------------------*/

ConvectionDiffusionEquation2d::ConvectionDiffusionEquation2d(const ParamFile* paramfile) 
{
  DataFormatHandler DFH;
  DFH.insert("visc" , &_visc , 0.01);
  DFH.insert("betax", &_bx   , 1.);
  DFH.insert("betay", &_by   , 1.);
  
  FileScanner FS(DFH,paramfile,"Equation");
};

/*---------------------------------------------------*/

void ConvectionDiffusionEquation2d::SetTimePattern(TimePattern& P) const 
{
  P.resize(GetNcomp(),GetNcomp());
  P(0,0) = 1.;
}

/*-----------------------------------------*/

void ConvectionDiffusionEquation2d::OperatorStrong(DoubleVector& b, const FemFunction& U)const
{
  b[0] += Convection(U[0]) - _visc*U[0].D();
}
 
/*---------------------------------------------------*/

double ConvectionDiffusionEquation2d::Convection(const TestFunction& N) const
{
  return betax()*N.x() + betay()*N.y();
}

/*---------------------------------------------------*/

void ConvectionDiffusionEquation2d::Form(VectorIterator b, const FemFunction& U, 
					 const TestFunction& N) const
{
  b[0] += _visc* (U[0].x()*N.x()+U[0].y()*N.y());
  b[0] += betax()* U[0].x()*N.m();
  b[0] += betay()* U[0].y()*N.m();
}

/*---------------------------------------------------*/

void ConvectionDiffusionEquation2d::Matrix(EntryMatrix& A, const FemFunction& U, 
					   const TestFunction& M, const TestFunction& N) const
{
  A(0,0) += _visc* (M.x()*N.x()+M.y()*N.y());
  A(0,0) += betax()* M.x()*N.m();
  A(0,0) += betay()* M.y()*N.m();
}
}
/*---------------------------------------------------*/
