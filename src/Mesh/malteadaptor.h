#ifndef  __MalteAdaptor_h
#define  __MalteAdaptor_h

#include  "nvector.h"
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

class MalteAdaptor
{
protected:

  typedef nvector<double>  dvector;
  typedef nvector<int>     ivector;

  const  nvector<double>&   eta;
  int    ppp, coarsening, refining, maxnodes, N;
  double etasum, gamma, alpha, beta, yfactor;

  double Expectation(double theta, double x) const;
  double Expectation(double thetax, double thetay, double x, double y) const;
  double ExpectationCoarsening(double theta, double x) const;
  void   refine_and_coarse(ivector& ref, ivector& coarse) const;

public:

  MalteAdaptor(const Gascoigne::ParamFile* pf, const dvector& eta);
  void coarse(nvector<int>& coarse) const;
  void refine(nvector<int>& ref) const;
  void refine(ivector& ref, ivector& coarse) const;
};

#endif
