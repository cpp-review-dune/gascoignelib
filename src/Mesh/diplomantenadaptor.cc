#include  "diplomantenadaptor.h"
#include  "compareclass.h"
#include  "giota.h"

#include  <fstream>

using namespace std;

/*-----------------------------------------*/

DiplomantenAdaptor::DiplomantenAdaptor(AdaptorData& _info, const dvector& _eta) :
  info(_info), eta(_eta)
{
  if (info.dim()==2) ppp = 3;
  else ppp=7;

  assert(info.rfactor()<=1);
  assert(info.rfactor()>0);
}

/*-----------------------------------------*/

void DiplomantenAdaptor::analyse() const
{
  double s = accumulate(eta.begin(),eta.end(),0.);
  //double s = 1.;
  double reduction = 1.-pow(0.5,info.local_conv());

  nvector<int> C(eta.size()); 
  iota(C.begin(),C.end(),0);

  typedef CompareObjectBigToSmall<nvector<double> >  CoC;

  sort(C.begin(),C.end(),CoC(eta));

  double theta = 0;
  ofstream file("hyundai");
  double dx = 1./eta.size();
  for (int i=0; i<eta.size(); i++)
    {
      double x = float(i)*dx;
      theta += eta[C[i]]*dx;
      double f = theta + eta[C[i]] * (x+1./ppp) + s/(reduction*eta.size());
      //theta += dx*dx;
      //double f = theta + dx * (x+1./ppp) + s/(reduction*eta.size());
      file << x << " " << f << std::endl;
      std::cout << x << " " << f << std::endl;
    }
  file.close();
}

/*-----------------------------------------*/

void DiplomantenAdaptor::MalteRefine(nvector<int>& ref) const
{
  if (eta.size()==0) return;

  double alpha = info.local_conv();  // lokale konvergenzordnung in h
  double beta  = info.global_conv(); // globale konvergenzordnung
  
  nvector<int> C(eta.size()); 
  iota(C.begin(),C.end(),0);

  typedef CompareObjectBigToSmall<nvector<double> >  CoC;

  sort(C.begin(),C.end(),CoC(eta));
  
  int n;
  for (n=0; n<eta.size(); n++)
    {
      if (eta[C[n]]==0.) break;
    }
  double etasum = accumulate(eta.begin(),eta.end(),0.);
  double mingoal = etasum*info.rfactor();
  double etaaktuell = etasum;
  double reduction = 1.-pow(0.5,alpha);

  int i;
  for (i=0; (i<n) && (etaaktuell>mingoal); i++)
    {
      etaaktuell -= eta[C[i]] * reduction;
    }
  int nelem = n + i*ppp;
  i++;

  double xlimit = etasum*pow(nelem,beta);
  for( ; i<n; i++)
    {
      etaaktuell -= eta[C[i]] * reduction;
      nelem  += ppp;
      double x = etaaktuell * pow(nelem,beta);
      if (x>xlimit) break;
    }
  ref.resize(i);
  for(int j=0; j<i; j++)
    {
      ref[j] = C[j];
    }
}

/*-----------------------------------------*/

void DiplomantenAdaptor::refine(nvector<int>& ref)
{
  int n = eta.size();

  if (n==0) return;


  double alpha = info.local_conv(); // konvergenzordnung
  double t     = 1.-pow(0.5,alpha);
  
  nvector<int> C(eta.size()); 
  iota(C.begin(),C.end(),0);

  typedef CompareObjectBigToSmall<nvector<double> >  CoC;

  std::sort(C.begin(),C.end(),CoC(eta));
  
  int i;
  for(i=0; i<n; i++)
    {
      if(eta[C[i]]==0.) break;
    }
  n = i;
  
  double etasum = accumulate(eta.begin(),eta.end(),0.);
  
  int nelem = n;

  double value = etasum*pow(nelem,info.global_conv());

  double goal = etasum*(info.rfactor()+(1.-info.rfactor())*pow(0.5,alpha));
  
  for(i=0; i<n; i++)
    {
      etasum -= t*eta[C[i]];
      nelem  += ppp;
      double newvalue = etasum*pow(nelem,info.global_conv());
      if ((newvalue>value)&&(etasum<goal)) break;

      value = newvalue;
    }
  if(i==0) i = n/10;  // 10% verfeinern

  ref.resize(i);
  for(int j=0; j<i; j++)
    {
      ref[j] = C[j];
    }
}

/*-----------------------------------------*/
