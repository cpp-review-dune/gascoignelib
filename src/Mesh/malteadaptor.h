#ifndef  __MalteAdaptor_h
#define  __MalteAdaptor_h

#include  "gascoigne.h"
#include  "adaptordata.h"
#include  "paramfile.h"

//
/// Minimizes E*L by global search,
/// where E = extrapolated error estimator
///       L = extrapolated costs
/// f(x) = [theta(1)+gamma*theta(x)] * [1+p*x]
/// f(x)  --> min
/// p = additional cells for each refined cell
/// x = fraction of cells to be refined
/// gamma = 2^(-alpha) -1   
/// alpha = local convergence rate (h)
/// theta(x) = int_0^x eta(t)dt
//

namespace Gascoigne
{
class MalteAdaptor
{
protected:

  const  DoubleVector&   eta;
  int    ppp, coarsening, refining, maxnodes, N;
  double etasum, gamma, alpha, beta, yfactor;

  double Expectation(double theta, double x) const;
  double Expectation(double thetax, double thetay, double x, double y) const;
  double ExpectationCoarsening(double theta, double x) const;
  void   refine_and_coarse(IntVector& ref, IntVector& coarse) const;

public:

  MalteAdaptor(const ParamFile* pf, const DoubleVector& eta);
  void coarse(IntVector& coarse) const;
  void refine(IntVector& ref) const;
  void refine(IntVector& ref, IntVector& coarse) const;
};
}

#endif
