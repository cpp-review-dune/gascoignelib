#ifndef __Equation_h
#define __Equation_h

#include  <vector>
#include  <set>
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

  typedef  nvector<double>        Vector;
  typedef  std::set<int>          IntSet;
  typedef Gascoigne::FemFunction  FemFunction;
  typedef Gascoigne::FemData      FemData;
  typedef Gascoigne::TestFunction TestFunction;

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

  virtual void OperatorStrong(Vector& b, const FemFunction& U) const { assert(0);} 

  virtual bool MatrixIsTransposed() const {return 0;}

  virtual void SetDeltaT(double k){assert(0);}
  virtual void SetTimePattern(TimePattern& P) const{assert(0);}
  virtual void SetTime(double time, double dt) const {_time = time; _dt = dt;}

  //
  // --------------------------------------
  //
  
  virtual void point(double h, const FemFunction& U, const Vertex2d& v) const {}
  virtual void point(double h, const FemFunction& U, const Vertex3d& v) const {}
  
  virtual void pointmatrix(double h, const FemFunction& U, const Vertex2d& v) const {
    point(h,U,v);
  }
  virtual void pointmatrix(double h, const FemFunction& U, const Vertex3d& v) const {
    point(h,U,v);
  }
  
  virtual void point(double h, const FemFunction& U, Gascoigne::FemData& QH, const Vertex2d& v) const {
    point(h,U,v);
  }
  virtual void point(double h, const FemFunction& U, Gascoigne::FemData& QH, const Vertex3d& v) const {
    point(h,U,v);
  }
  
  virtual void pointmatrix(double h, const FemFunction& U, Gascoigne::FemData& QH, const Vertex2d& v) const {
    point(h,U,QH,v);
  }
  virtual void pointmatrix(double h, const FemFunction& U, Gascoigne::FemData& QH, const Vertex3d& v) const {
    point(h,U,QH,v);
  }

  virtual void pointboundary
    (double h, const FemFunction& U, const Vertex2d& v, 
     const Vertex2d& n) const {assert(0);}

  virtual void pointboundary 
    (double h, const FemFunction& U, const Vertex3d& v, 
     const Vertex3d& n) const {assert(0);}

  virtual IntSet GetBoundaryColors() const { return IntSet();}

  virtual void SetParameterData(Gascoigne::LocalData& q) const { }

  //
  // ---------------------------------------------
  //

  virtual void Form(Gascoigne::VectorIterator b, const FemFunction& U, const TestFunction& N) const {assert(0);}

  virtual void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const {assert(0);}

  virtual void BoundaryResidual 
    (int col, Vector& b, const FemFunction& U, 
     const TestFunction& N) const {assert(0);}

  virtual void BoundaryMatrix
    (int col, EntryMatrix& D, const FemFunction& U, 
     const TestFunction& M, const TestFunction& N) const {assert(0);}
};


#endif
