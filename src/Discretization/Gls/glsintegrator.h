#ifndef  __GlsIntegrator_h
#define  __GlsIntegrator_h


#include  "basicintegrator.h"
#include  "integrationformula.h"

/*-----------------------------------------*/

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments GlsIntegrator

////
////
/////////////////////////////////////////////

template<int DIM>
class GlsIntegrator : public BasicIntegrator
{
private:

protected:

  IntegrationFormulaInterface* IF;

  const IntegrationFormulaInterface& FormFormula() const { return *IF;}
  double Volume2MeshSize(double vol) const { return pow(vol,1./float(DIM));}

 public:


  GlsIntegrator<DIM>();
  ~GlsIntegrator<DIM>() {}

  std::string GetName() const {return "Gls";}

  void Rhs(const DomainRightHandSide& RHS, LocalVector& F, const FemInterface& FEM, const LocalNodeData& Q) const {};
  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector&U, const LocalNodeData& Q) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const;
};
}

/*-----------------------------------------*/


#endif
