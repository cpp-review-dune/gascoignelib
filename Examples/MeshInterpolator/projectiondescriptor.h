#ifndef  __ProjectionDescriptor_h
#define  __ProjectionDescriptor_h

#include  "problemdescriptorbase.h"

/*---------------------------------------------------*/

class ProjectionEquation : public Gascoigne::Equation
{
 public:

  ProjectionEquation() : Gascoigne::Equation() { }
  ~ProjectionEquation() { }

  int GetNcomp() const { return 1; }
  std::string GetName() const { return "ProjectionEquation"; }

  void Form(Gascoigne::VectorIterator b, const Gascoigne::FemFunction& U, 
	    const Gascoigne::TestFunction& N) const
  {
    b[0] += U[0].m() * N.m();
  }
  void Matrix(Gascoigne::EntryMatrix& A, const Gascoigne::FemFunction& U, 
	      const Gascoigne::TestFunction& M, const Gascoigne::TestFunction& N) const
  {
    A(0,0) += M.m() * N.m();
  }
};

/*---------------------------------------------------*/

class ProjectionProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "Projection";}
  void BasicInit(const Gascoigne::ParamFile* pf) 
  {
    GetEquationPointer() = new ProjectionEquation();
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

#endif
