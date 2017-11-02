
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

class MyVPS : public DirichletData
{
  protected:
    double vmean;

  public:
    MyVPS(const ParamFile* pf)
      {
        DataFormatHandler DFH;
        DFH.insert("vmean" ,  &vmean , 0.0);
        FileScanner FS(DFH, pf, "Equation");
      }

    std::string GetName() const {return "MyVPS";}

    void operator()(DoubleVector& b, const Vertex2d& v, int color) const
    {
      b.zero();
      if (color==8)
        b[1] = 8.0 * (1.0+ tanh( 8.0 * (__TIME - 0.5))) * v.x()*v.x()*(1.0-v.x())*(1.0-v.x());

      //stress
      b[3]=1.0;
      b[4]=0.0;
      b[5]=1.0;

    }

};

// -----------------------------------------
class VPSProblem : public ProblemDescriptorBase
{
  public:

    std::string GetName() const {return "vps";}
    void BasicInit(const ParamFile* pf)
    {
      GetParamFilePointer() = pf;
      GetEquationPointer()      = new VPSEQ(GetParamFile());
      GetDirichletDataPointer() = new MyVPS(GetParamFile());

      ProblemDescriptorBase::BasicInit(pf);
    }
};


#endif
