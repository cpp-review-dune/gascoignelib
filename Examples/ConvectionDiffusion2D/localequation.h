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

class LocalEquation : public GlsEquation
{
private:

  mutable GlsStabilization ST;
  mutable const FemFunction* q;

protected:

  mutable double visc, sigma;

  double betax() const {return (*q)[1].m();}
  double betay() const {return (*q)[2].m();}

public:


//
////  Con(De)structor 
//

  LocalEquation(const std::string& filename);
  ~LocalEquation() {}

  std::string GetName() const {return "Local";}

  int  ncomp() const {return 1;}

  void point(double h, const FemFunction& U, FemData& Q, const Vertex2d& v) const;
  void glspoint(double h, const FemFunction& U, FemData& Q, const Vertex2d& v) const {
    LocalEquation::point(h,U,Q,v);
  }

  void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
  
  void Matrix(EntryMatrix& D, const FemFunction& U, const DerivativeVector& M, const TestFunction& N) const;


  void L(nvector<double>& dst, const FemFunction& U) const;
  void S(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;
  void LMatrix(nmatrix<double>& dst, const FemFunction& U, const TestFunction& N) const;
  
};


#endif
