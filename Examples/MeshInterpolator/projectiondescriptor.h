#ifndef  __ProjectionDescriptor_h
#define  __ProjectionDescriptor_h

#include  "problemdescriptorbase.h"
#include  "dirichletdatabycolor.h"
#include  "domainrighthandside.h"

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

class ProjectionRightHandSide : public Gascoigne::DomainRightHandSide
{
 protected:
   mutable Gascoigne::FemFunction _U;
   int _ncomp;

 public:

  ProjectionRightHandSide(const Gascoigne::Equation* EQ) : Gascoigne::DomainRightHandSide() 
    { 
      _ncomp = EQ->GetNcomp(); 
    }
  ~ProjectionRightHandSide() { }
  
  int GetNcomp() const { return _ncomp; }
  std::string GetName() const { return "ProjectionRightHandSide"; }

  void SetFemData(Gascoigne::FemData& q) const
  {
    assert(q.count("U")==1);
    _U = q["U"];
  }

  void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, 
		  const Gascoigne::Vertex2d& v) const 
  {
    b[0] += _U[0].m() * N.m();
  }
};

/*---------------------------------------------------*/

class ProjectionProblemDescriptor : public Gascoigne::ProblemDescriptorBase
{
 public:
  
  std::string GetName() const {return "Projection";}
  void BasicInit(const Gascoigne::ParamFile* pf) 
  {
    GetEquationPointer()      = new ProjectionEquation();
    GetRightHandSidePointer() = new ProjectionRightHandSide(GetEquation());
    
    ProblemDescriptorBase::BasicInit(pf);
  }
};

/*---------------------------------------------------*/

#endif
