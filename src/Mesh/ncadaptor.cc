#include  "ncadaptor.h"
#include  "compareclass.h"
#include  "giota.h"
#include  "filescanner.h"
#include  <climits>


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
NCAdaptor::NCAdaptor(const ParamFile* paramfile, const DoubleVector& _eta) 
  :  eta(_eta)
{
  DataFormatHandler DH;

#ifdef __OLDCOMPILER__
  DH.insert("n"   ,& _n,INT_MAX);
#else
  DH.insert("n"   ,& _n,numeric_limits<int>::max());
#endif

  DH.insert("p"  ,& _p,0.9);
  FileScanner FS(DH, paramfile, "Adaptor");

  etasum = accumulate(eta.begin(),eta.end(),0.);
}

/*-----------------------------------------*/

void NCAdaptor::refine(IntVector& ref, IntVector& coars) const
{
  int n = eta.size();
  IntVector C(n); 
  iota(C.begin(),C.end(),0);
  typedef CompareObjectBigToSmall<DoubleVector >  CoC;
  sort(C.begin(),C.end(),CoC(eta));

  double eta_limit = _p*etasum;

  int    i    = 0;
  double done = 0.;

  ref.clear();
  while((done<eta_limit) && (i<_n))
    {
      done += fabs(eta[C[i]]);
      ref.push_back(C[i]);
      i++;
    }
}
}
