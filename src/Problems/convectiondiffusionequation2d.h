#ifndef __ConvectionDiffusionEquation_h
#define __ConvectionDiffusionEquation_h

#include "equation.h"
#include "paramfile.h"

/*---------------------------------------------------*/

namespace Gascoigne{

class ConvectionDiffusionEquation2d : public virtual Equation
{
protected:

  double _visc, _bx, _by;

  double betax() const {return _bx;}
  double betay() const {return _by;}
  double Convection(const TestFunction& N) const;

public:

  ConvectionDiffusionEquation2d(const ParamFile* paramfile);
  ~ConvectionDiffusionEquation2d() {}
  void SetTimePattern(TimePattern& P) const;

  std::string GetName() const {return "ConvectionDiffusionEquation";}

  int  GetNcomp()    const { return 1;}

  void OperatorStrong(DoubleVector& b, const FemFunction& U) const;

  void point(double h, const FemFunction& U, FemData& Q, 
	     const Vertex2d& v) const {}

  void Form(VectorIterator b, const FemFunction& U, 
	    const TestFunction& N) const;

  void Matrix(EntryMatrix& A, const FemFunction& U, 
	      const TestFunction& M, const TestFunction& N) const;
};
}
/*---------------------------------------------------*/

#endif
