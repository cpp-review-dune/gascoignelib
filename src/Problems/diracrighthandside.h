#ifndef  __DiracRightHandSide_h
#define  __DiracRightHandSide_h

#include  "righthandsidedata.h"

/*-----------------------------------------*/


namespace Gascoigne
{
class DiracRightHandSide : public RightHandSideData
{
protected:

  std::vector<Vertex2d>  v0;
  DoubleVector   w;
  int mycomp, ncomp;
  double eps;

public:

  DiracRightHandSide(int n, int m, const Vertex2d& v)
    : ncomp(n), mycomp(m)
    {
      eps = 1e-4;
      v0.resize(1,v);
    };
  DiracRightHandSide(int n, int m, const std::vector<Vertex2d>& v,
		     const DoubleVector& w0)
    : ncomp(n), mycomp(m), v0(v), w(w0)
    {
      eps = 1e-4;
    };
  ~DiracRightHandSide(){};
  
  std::string GetName() const {return "DiracRightHandSide";}
  int size() const {return ncomp;}

  double operator()(int c, const Vertex2d& x) const
    {
      if (c==mycomp)
	{
	  for (int i=0; i<v0.size(); i++)
	    {
	      Vertex2d w(v0[i]);
	      w -= x;  
	      if (w.norm()<eps) { return 1.;}
	    }
	}
      return 0.;
    }
  
  int GetNcomp() const { return mycomp; }
  
  const DoubleVector& GetWeights()    const { return w; }
  const std::vector<Vertex2d>& GetPoints2d() const 
    { return v0; }
};
}

#endif
