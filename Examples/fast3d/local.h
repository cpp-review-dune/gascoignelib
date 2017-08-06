
#ifndef  __local_h
#define  __local_h

#include  "fsi.h"
#include  "dirichletdata.h"
#include  "problemdescriptorbase.h"
#include  "domainrighthandside.h"
#include  "componentinformationbase.h"
#include <stdio.h>
#include <stdlib.h>


using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

extern double __TIME;


template<int DIM>
class FSI_CI : public ComponentInformationBase
{
 public:
  void BasicInit(const ParamFile* pf) {}

  std::string GetName() const { return "FSI CI"; }
  
  const int GetNScalars() const  {return 2*DIM+1;}

  void      GetScalarName   (int i, std::string& s_name) const
  {
    if (i==0) s_name="p";
    if (DIM==2)
      {
	if (i==1) s_name = "vx";
	if (i==2) s_name = "vy";
	if (i==3) s_name = "ux";
	if (i==4) s_name = "uy";	
      }
    if (DIM==3)
      {
	if (i==1) s_name = "vx";
	if (i==2) s_name = "vy";
	if (i==3) s_name = "vz";
	if (i==4) s_name = "ux";
	if (i==5) s_name = "uy";
	if (i==6) s_name = "uz";	
      }
  }
  
  const int GetNVectors     () const {return 2;}
  void      GetVectorName   (int i, std::string& s_name) const
  {
    if (i==0) s_name = "V";
    else if (i==1) s_name = "U";
    else abort();
  }
  void      GetVectorIndices(int i, fixarray<3,int>& fa_vectorindices) const
  {
    if (i==0) 
      {
	fa_vectorindices[0]=1;
	fa_vectorindices[1]=2;
	fa_vectorindices[2]=-1;
	if (DIM==3) fa_vectorindices[2]=3;
      }
    else
      {
	fa_vectorindices[0]=DIM+1;
	fa_vectorindices[1]=DIM+2;
	fa_vectorindices[2]=-1;
	if (DIM==3) fa_vectorindices[2]=DIM+3;
      }
  }
  
};


class MyDD : public DirichletData
{
 protected:
  double vmean;
  
 public:
  MyDD(const ParamFile* pf)
    {
      DataFormatHandler DFH;
      DFH.insert("vmean" ,    &vmean , 0.0);
      FileScanner FS(DFH, pf, "Equation");
    }
  
  std::string GetName() const {return "MyDD";}
  
  void operator()(DoubleVector& b, const Vertex2d& v, int color) const
  {
    b.zero();

    double teff = __TIME;
    int Ti=0;
    while (teff>5.0)
      {
	++Ti;
	teff -=5.0;
      }
    
    // transition by 0.1
    double vnew = 1.0+0.1*Ti;
    double vold = 1.0+0.1*(Ti-1.0);
    if (vold<1.0) vold = 1.0;
    double sc = 1.0;
    if (teff<0.1)
      sc = 0.5*(1.0-cos(M_PI*teff/0.1));

    //    sc = 1.0;


    //    ////// Std-Werte v=2
    vnew = vold = 1.15; // re = 130
    


    double veff = vold + sc*(vnew-vold);
    //    cout << veff << endl;
    
    if (color==0)
      b[1] += v.y()*(0.41-v.y())/0.205/0.205 * veff * 1.5;
    

    return;
    



    //    double sc = 1.0;
    // transition from zero to vmean at time t=0 (2 sec)
    if (__TIME<2.0) 
      {
	sc = 0.5*(1.0-cos(M_PI*__TIME/2.0));
	if (color==0)
	  b[1] += v.y()*(0.41-v.y())/0.205/0.205 * vmean * 1.5 * sc;	
      }
    else if (__TIME<10.0)
      {
	sc = 1.0;
	if (color==0)
	  b[1] += v.y()*(0.41-v.y())/0.205/0.205 * vmean * 1.5 * sc;
      }
    else if (__TIME<11.0) // transition from 2.0 to vmean in t in [10,11]
      {
	sc = 2.0 + 0.5*(1.0-cos(M_PI*(__TIME-10.0)))*(vmean-2.0);
	
	if (color==0)
	  b[1] += v.y()*(0.41-v.y())/0.205/0.205 *  1.5 * sc;
	
      }
    else
      {
	sc = 1.0;
	if (color==0)
	  b[1] += v.y()*(0.41-v.y())/0.205/0.205 * vmean * 1.5 * sc;
      }
    
    
    
    

    // transition from 2.0 to vmean at time t=10 (1 sec)
    if (__TIME>10.0)
      {	
	sc = 0.5*(1.0-cos(M_PI*__TIME/2.0));
	
      }
    
    
    
  }
};

class MyDD3d : public DirichletData
{
 protected:
  double vmean;
  
 public:
  MyDD3d(const ParamFile* pf)
    {
      DataFormatHandler DFH;
      DFH.insert("vmean" ,    &vmean , 0.0);
      FileScanner FS(DFH, pf, "Equation");
    }
  
  std::string GetName() const {return "MyDD";}
  
  void operator()(DoubleVector& b, const Vertex3d& v, int color) const
  {
    b.zero();




    double teff = __TIME;
    
    double vvmax = 2.25;
    
    double sc = vvmax;
    if ((teff>=0.0) && (teff<1.0) )  sc =0.0 +  vvmax/2.0*(1.0-cos(M_PI*teff));
    
    if (color==0)
      b[1] += v.y()*(0.41-v.y())/0.205/0.205 * v.z()*(0.41-v.z())/0.205/0.205 * vmean * sc * 9.0/4.0;
    
  }
};

// -----------------------------------------

class ProblemDescriptor2d : public ProblemDescriptorBase
{
 public:
    
  std::string GetName() const {return "fsi";}
  void BasicInit(const ParamFile* pf)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new FSI<2>(GetParamFile());
    //    GetBoundaryEquationPointer() = new FSI<2>(GetParamFile());
    GetDirichletDataPointer() = new MyDD(GetParamFile());
    
    ProblemDescriptorBase::BasicInit(pf);

    GetComponentInformationPointer() = new FSI_CI<2>;
    
  }
};


class ProblemDescriptor3d : public ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "fsi";}
  void BasicInit(const ParamFile* pf)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new FSI<3>(GetParamFile());
    //    GetBoundaryEquationPointer() = new FSI<3>(GetParamFile());
    GetDirichletDataPointer() = new MyDD3d(GetParamFile());
    
    ProblemDescriptorBase::BasicInit(pf);
    
    GetComponentInformationPointer() = new FSI_CI<3>;
    
  }
};






#endif
