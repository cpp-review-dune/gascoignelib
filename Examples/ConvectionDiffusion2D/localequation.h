#ifndef  __LocalEquation_h
#define  __LocalEquation_h


/////////////////////////////////////////////
////
////@brief
////  ... comments LocalEquation

////
////
/////////////////////////////////////////////

#include  "glsequation.h"
#include  "glsstabilization.h"

class LocalEquation : public Gascoigne::GlsEquation
{
private:

  mutable Gascoigne::GlsStabilization ST;
  mutable const Gascoigne::FemFunction* q;

protected:

  mutable double visc, sigma;

  double betax() const {return (*q)[1].m();}
  double betay() const {return (*q)[2].m();}

public:


//
////  Con(De)structor 
//

  LocalEquation(const Gascoigne::ParamFile* paramfile);
  ~LocalEquation() {}

  std::string GetName() const {return "Local";}

  int  ncomp() const {return 1;}

  void point(double h, const Gascoigne::FemFunction& U, Gascoigne::FemData& Q, const Gascoigne::Vertex2d& v) const;
  void glspoint(double h, const Gascoigne::FemFunction& U, Gascoigne::FemData& Q, const Gascoigne::Vertex2d& v) const {
    LocalEquation::point(h,U,Q,v);
  }

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;
  
  void Matrix(Gascoigne::EntryMatrix& D, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const;


  void L(Gascoigne::DoubleVector& dst, const Gascoigne::FemFunction& U) const;
  void S(Gascoigne::nmatrix<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;
  void LMatrix(Gascoigne::nmatrix<double>& dst, const Gascoigne::FemFunction& U, const Gascoigne::TestFunction& N) const;
  
};


#endif
