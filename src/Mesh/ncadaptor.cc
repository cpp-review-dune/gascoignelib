#include  "ncadaptor.h"
#include  "compareclass.h"
#include  "giota.h"
#include  "filescanner.h"
#include  <climits>

/*-----------------------------------------*/

NCAdaptor::NCAdaptor(const std::string& filename, const nvector<double>& _eta) 
  :  eta(_eta)
{
  DataFormatHandler DH;

#ifdef __OLDCOMPILER__
  DH.insert("n"   ,& _n,INT_MAX);
#else
  DH.insert("n"   ,& _n,std::numeric_limits<int>::max());
#endif

  DH.insert("p"  ,& _p,0.9);
  FileScanner FS(DH,filename,"Adaptor");

  etasum = accumulate(eta.begin(),eta.end(),0.);
}

/*-----------------------------------------*/

void NCAdaptor::refine(nvector<int>& ref, nvector<int>& coars) const
{
  int n = eta.size();
  nvector<int> C(n); 
  iota(C.begin(),C.end(),0);
  typedef CompareObjectBigToSmall<nvector<double> >  CoC;
  std::sort(C.begin(),C.end(),CoC(eta));

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
