#include  "optadaptor.h"
#include  "compareclass.h"

#include  <fstream>
#include  "giota.h"

/*********************************************************************/

OptAdaptor::OptAdaptor
(AdaptorData& inf, nvector<double>& e,const nvector<double>& v) 
: info(inf), eta(e), vol(v)
{
  p = info.dim();
  d = 2;
  if(d==2)
    {
      p2 = 4;
      p4 = 16;
    }
  else
    {
      p2 = 8;
      p4 = 64;
    }
  dd = d;
  pp = p;
  factor = 1.;
  
  prepare();
}

/*********************************************************************/

void OptAdaptor::prepare()
{
  n_aimed = GascoigneMath::max_int(1,int(info.rfactor()*info.ncells()));
  n_aimed = GascoigneMath::min_int(info.maxnodes(),n_aimed);

  info.reset();

  int alpha = int(pp+dd);

  info.eta() = eta.norm_l1();

  double integral = 0.;
  for (int i=0; i<eta.size(); i++)
    {
      if (vol[i]>0)
	{
	  double h  = pow(vol[i],1./dd);
	  eta[i]   /= pow(h,alpha);
	  integral += vol[i] * pow(eta[i],dd/alpha);
	}
    }
  //info.cs() = integral;
  co = pow(integral/n_aimed,1./dd);

  double min = 1.e8, max = 0.;
  for (int i=0; i<eta.size(); i++)
    {
      if (vol[i]>0.)
	{
	  double h       = pow(vol[i],1./dd);
	  double h_aimed = co * pow(eta[i],-1./alpha);
	  double fac = h/h_aimed;
	  eta[i] = fac;
	  min = GascoigneMath::min(min,fac);
	  max = GascoigneMath::max(max,fac);
	}
    }
  info.minf() = min;
  info.maxf() = max;
  marge = n_aimed - info.ncells();
}

/*********************************************************************/

void OptAdaptor::coarse(nvector<int>& coarselist)
{
  //  eta = h/h_aimed

  if (info.cfactor()<=0.) return;

  for (int i=0; i<eta.size(); i++)
    {
      if (vol[i]>0.)
	{
	  if (info.cfactor() * eta[i] < 1.)
	    {
	      coarselist.push_back(i);
	    }
	}
    }
  return;

  /*************************** old *******************/

  for (int i=0; i<eta.size(); i++)
    {
      int level = 0; // ???????????
      if ((vol[i]>0) && (level>1))
	{
	  int test = 0;
	  int nchilds = 0;  //  ??????????????
	  for (int ch=0; ch<nchilds; ch++)
	    {
	      int child = 0;//  ??????????????
	      if (info.cfactor() * eta[child] < 1)
		{
		  test++;
		}
	    }
	  if (test==nchilds)
	    {
	      info.nc() += test;
	      marge     += test;
	      for (int ch=0; ch<nchilds; ch++)
		{
		  int child = 0;//  ??????????????
		  coarselist.push_back(child);
		}	      
	    }
	}
    }
}

/*********************************************************************/

void OptAdaptor::refine(nvector<int>& reflist)
{
  reflist.resize(0);

  nvector<int> C(eta.size()); iota(C.begin(),C.end(),0);
  std::sort(C.begin(),C.end(),CompareObjectBigToSmall<nvector<double> > (eta)); 
  
  int i = 0; used = 0;

  while((marge>used) && (i<eta.size()))
    {
      reflist.push_back(C[i]);
      i++; 
      used += 4;
    }

  info.nr() += i;
  marge     -= used;
}

/*********************************************************************/


/*********************************************************************/

void OptAdaptor::RefineGnuplot(nvector<int>& reflist)
{
  reflist.resize(0);

  nvector<int> C(eta.size()); 
  iota(C.begin(),C.end(),0);

  typedef CompareObjectBigToSmall<nvector<double> >  CoC;

  std::sort(C.begin(),C.end(),CoC(eta)); 
  
  int i = 0; used = 0;

  std::ofstream cmdfile("cmd.gpl");
  cmdfile << "plot  \"eta.dat\" using 1:2 title \"eta\" with lines lw 0" << std::endl;
  cmdfile << " pause -1" << std::endl;
  cmdfile << "plot  \"eta.dat\" using 1:3 title \"delta\" with lines lw 0"<<std::endl;
  cmdfile << " pause -1" << std::endl;
  cmdfile.close();

  std::ofstream file("eta.dat");

  for(int ii=0;ii<eta.size();ii++)
    {
      double e = eta[C[ii]];      
      if( (ii>=1) && (e>0) )
	{
 	  double delta = eta[C[ii-1]]-e;
	  file << ii << "\t" << e <<  "\t"<<  delta << std::endl;
	}
    }

  file.close();
  system("gnuplot cmd.gpl");

  while((marge>used) && (i<eta.size()))
    {
     reflist.push_back(C[i]);
      i++; 
      used += 4;
    }

  info.nr() += i;
  marge     -= used;
}

/*********************************************************************/

