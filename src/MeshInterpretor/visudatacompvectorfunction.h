#ifndef __visudatacompvectorfunction_h
#define __visudatacompvectorfunction_h

#include  "visudatacompvector.h"

/*-------------------------------------------------------------------------*/

class SolutionFunction
{
 public:
  virtual ~SolutionFunction(){}
  virtual double operator()(int i, const CompVector<double>&) const =0;
};

/*-------------------------------------------------------------------------*/

class MySolutionFunction : public SolutionFunction
{
 public:

  double operator()(int i, const CompVector<double>& U) const;
};

/***************************************************************/

// mit noch weitere Funktionen die man ausrechnen moechte

class VisuDataCompVectorFunction : public VisuDataCompVector
{
 protected:

  std::vector<const SolutionFunction*>    outcome;
  std::pair<int,int> VisuDataCompVectorFunction::GetIndex(int c) const;

 public:

  VisuDataCompVectorFunction();
  VisuDataCompVectorFunction(const CompVector<double>& v);

  std::vector<const SolutionFunction*> Outcome() { return outcome;}

  int    visucomp()     const;
  double visudata(int i,int c) const;
};

/***************************************************************/

#endif
