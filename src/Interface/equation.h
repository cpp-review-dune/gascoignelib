#ifndef __Equation_h
#define __Equation_h

#include  <string>

#include  "gascoigne.h"
#include  "vertex.h"
#include  "entrymatrix.h"
#include  "timepattern.h"
#include  "exactsolution.h"
#include  "paramfile.h"

//////////////////////////////////////////////
///
///@brief
/// Interface class for Equation

///
///
//////////////////////////////////////////////

/*-------------------------------------------------------------------------*/

class Equation
{
 private:
  
  mutable double _time, _dt;
  
 protected:
  double GetTime() const {return _time;}
  double GetTimeStep() const {return _dt;}

 public:
  
  //
  // Constructors
  //

  Equation() : _time(0.) {}
  Equation(const Equation& E) {assert(0);}
  virtual ~Equation() {}

  virtual std::string GetName() const=0;

  virtual int  ncomp() const=0;

  virtual void OperatorStrong(Gascoigne::DoubleVector& b, const Gascoigne::FemFunction& U) const { assert(0);} 

  virtual bool MatrixIsTransposed() const {return 0;}

  virtual void SetDeltaT(double k){assert(0);}
  virtual void SetTimePattern(TimePattern& P) const{assert(0);}
  virtual void SetTime(double time, double dt) const {_time = time; _dt = dt;}

  //
  // --------------------------------------
  //
  
  virtual void point(double h, const Gascoigne::FemFunction& U, const Vertex2d& v) const {}
  virtual void point(double h, const Gascoigne::FemFunction& U, const Vertex3d& v) const {}
  
  virtual void pointmatrix(double h, const Gascoigne::FemFunction& U, const Vertex2d& v) const {
    point(h,U,v);
  }
  virtual void pointmatrix(double h, const Gascoigne::FemFunction& U, const Vertex3d& v) const {
    point(h,U,v);
  }
  
  virtual void point(double h, const Gascoigne::FemFunction& U, Gascoigne::FemData& QH, const Vertex2d& v) const {
    point(h,U,v);
  }
  virtual void point(double h, const Gascoigne::FemFunction& U, Gascoigne::FemData& QH, const Vertex3d& v) const {
    point(h,U,v);
  }
  
  virtual void pointmatrix(double h, const Gascoigne::FemFunction& U, Gascoigne::FemData& QH, const Vertex2d& v) const {
    point(h,U,QH,v);
  }
  virtual void pointmatrix(double h, const Gascoigne::FemFunction& U, Gascoigne::FemData& QH, const Vertex3d& v) const {
    point(h,U,QH,v);
  }

  virtual void pointboundary
    (double h, const Gascoigne::FemFunction& U, const Vertex2d& v, 
     const Vertex2d& n) const {assert(0);}

  virtual void pointboundary 
    (double h, const Gascoigne::FemFunction& U, const Vertex3d& v, 
     const Vertex3d& n) const {assert(0);}

  virtual Gascoigne::IntSet GetBoundaryColors() const { return Gascoigne::IntSet();}

  virtual void SetParameterData(Gascoigne::LocalParameterData& q) const { }

  //
  // ---------------------------------------------
  //

  virtual void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const {assert(0);}

  virtual void Matrix(EntryMatrix& A, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const {assert(0);}

  virtual void BoundaryResidual(int col, Gascoigne::DoubleVector& b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const {assert(0);}

  virtual void BoundaryMatrix(int col, EntryMatrix& D, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const {assert(0);}
};


#endif
