#ifndef  __ProblemDescriptorInterface_h
#define  __ProblemDescriptorInterface_h

#include  "gostream.h"
#include  "stringutil.h"
#include  "filescanner.h"

#include  "boundarymanager.h"
#include  "equation.h"
#include  "dirichletdata.h"
#include  "boundaryrighthandside.h"
#include  "boundaryequation.h"
#include  "exactsolution.h"
#include  "boundarymanager.h"
#include  "componentinformation.h"


namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments ProblemDescriptorInterface

  ///
  ///
  /////////////////////////////////////////////

  class ProblemDescriptorInterface
  {
    private:
      
    protected:

    public:
      ProblemDescriptorInterface() {}
      virtual ~ProblemDescriptorInterface() {}
  
      virtual void BasicInit(const ParamFile* pf) {}

      virtual std::string GetName() const=0;
      virtual std::ostream& OutputSettings(std::ostream& os) const=0;
      virtual void SetTime(double time, double dt) const=0;

      virtual const ParamFile* GetParamFile() const=0;

      virtual const Application*            GetRightHandSide        () const=0;
      virtual const BoundaryRightHandSide*  GetBoundaryRightHandSide() const=0;
      virtual const Equation*               GetEquation             () const=0;
      virtual const BoundaryEquation*       GetBoundaryEquation     () const=0;
      virtual const DirichletData*          GetDirichletData        () const=0;
      virtual const Application*            GetInitialCondition     () const=0;
      virtual const ExactSolution*          GetExactSolution        () const=0;
      virtual const BoundaryManager*        GetBoundaryManager      () const=0;
      virtual const ComponentInformation*   GetComponentInformation () const=0;
  };
}

#endif
