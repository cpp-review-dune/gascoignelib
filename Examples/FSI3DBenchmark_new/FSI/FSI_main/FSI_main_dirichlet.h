
#ifndef  __FSI_main_dirichlet_h
#define  __FSI_main_dirichlet_h

#include  "dirichletdata.h"
#include  "paramfile.h"
#include  "filescanner.h"


using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

extern double __TIME;


class FSI_main_MyDD : public DirichletData
{
protected:
  double vmean;
  
public:
  FSI_main_MyDD(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("vmean" ,    &vmean , 0.0);
    FileScanner FS(DFH, pf, "Equation");
  }
  
  std::string GetName() const {return "MyDD";}
  
  void operator()(DoubleVector& b, const Vertex2d& v, int color) const
  {
    b.zero();
   
    // double teff = __TIME;
    // int Ti=0;
    // while (teff>5.0)
    //   {
    // 	++Ti;
    // 	teff -=5.0;
    //   }
    
    // // transition by 0.1
    // double vnew = 1.0+0.1*Ti;
    // double vold = 1.0+0.1*(Ti-1.0);
    // if (vold<1.0) vold = 1.0;
    // double sc = 1.0;
    // if (teff<0.1)
    //   sc = 0.5*(1.0-cos(M_PI*teff/0.1));

    // //    sc = 1.0;


    // //    ////// Std-Werte v=2
    // vnew = vold = 1.15; // re = 130
    


    // double veff = vold + sc*(vnew-vold);
    // //    cout << veff << endl;
    
    // if (color==0)
    //   b[1] += v.y()*(0.41-v.y())/0.205/0.205 * veff * 1.5;
    

    // return;
    



        double sc = 1.0;
    // // transition from zero to vmean at time t=0 (2 sec)
    if (__TIME<2.0) 
       {
     		sc = 0.5*(1.0-cos(M_PI*__TIME/2.0));
       	if (color==0)
     	  b[1] += v.y()*(0.41-v.y())/0.205/0.205 * vmean * 1.5 * sc;	
       }
    // else if (__TIME<10.0)
    //   {
    // 	sc = 1.0;
    // 	if (color==0)
    // 	  b[1] += v.y()*(0.41-v.y())/0.205/0.205 * vmean * 1.5 * sc;
    //   }
    // else if (__TIME<11.0) // transition from 2.0 to vmean in t in [10,11]
    //   {
    // 	sc = 2.0 + 0.5*(1.0-cos(M_PI*(__TIME-10.0)))*(vmean-2.0);
	
    // 	if (color==0)
    // 	  b[1] += v.y()*(0.41-v.y())/0.205/0.205 *  1.5 * sc;
	
    //   }
     else
      {
     		sc = 1.0;
     		if (color==0)
     	  	b[1] += v.y()*(0.41-v.y())/0.205/0.205 * vmean * 1.5 * sc;
      }
    
    
    
    

    // // transition from 2.0 to vmean at time t=10 (1 sec)
    // if (__TIME>10.0)
    //   {	
    // 	sc = 0.5*(1.0-cos(M_PI*__TIME/2.0));
	
    //   }
    
    
    
  }
};

class FSI_main_MyDD3d : public DirichletData
{
protected:
  double vmean;
  
public:
  FSI_main_MyDD3d(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("vmean" ,    &vmean , 0.0);
    FileScanner FS(DFH, pf, "Equation");
  }
  
  std::string GetName() const {return "MyDD";}
  
  void operator()(DoubleVector& b, const Vertex3d& v, int color) const
  {
    double sc = 1.0;
    if (__TIME<2.0) sc = 0.5*(1.0-cos(M_PI*__TIME/2.0));
    
    if (color==0)
      {
      	b[1] = 9.0/4.0*v.y()*(0.41-v.y())/0.205/0.205 * v.z()*(0.41-v.z())/0.205/0.205 * vmean * sc ;
      	//b[1] = v.y()*(0.4-v.y())/0.4/0.4*4*1.5* vmean * sc ;
      }
  }
};

#endif
