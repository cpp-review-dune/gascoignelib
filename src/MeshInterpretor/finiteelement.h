#ifndef __FiniteElement_h
#define __FiniteElement_h

#include  "feminterface.h"


/////////////////////////////////////////////
///
///@brief
///  FE based on Transformation (TRAFO) and Rerefernzelement (BASE)

///
///
/////////////////////////////////////////////


/*-----------------------------------------------------*/

template<int DIM, int BDIM, class TRAFO, class BASE>
class FiniteElement : public FemInterface
{
 protected:

  TRAFO                   T;
  BASE                    B;
  mutable std::vector<Vertex<DIM> >    grad;
  mutable double                  det;

  virtual void ComputeGrad() const;

 public:

  FiniteElement();

  std::string GetName() const {return "FiniteElement";}

  int    n()          const { return B.n(); }
  double N   (int i)  const { return B.phi(i); }
  double N_x (int i)  const { return grad[i].x(); }
  double N_y (int i)  const { return grad[i].y(); }
  double N_z (int i)  const { return grad[i].z(); }
  double J()          const { return det;}
  double G()          const { return T.G();}

  void x     (Vertex<DIM>& v) const { v = T.x();}
  void normal(Vertex<DIM>& v) const { v = T.normal();};

  void point(const Vertex<DIM>&) const;
  void point_boundary(int ie, const Vertex<BDIM>& s1) const;
  /// depreciated
  void init (const Matrix& M) { T.init(M); }
  void ReInit(const Matrix& M) const { T.ReInit(M); }

  void SetNx(int n) {T.SetNx(n);}

  void  init_test_functions(Gascoigne::TestFunction& Phi, double w, int i) const;
};

/*-----------------------------------------------------*/

#endif
