#ifndef  __GlsIntegrator_h
#define  __GlsIntegrator_h


/////////////////////////////////////////////
////
////@brief
////  ... comments GlsIntegrator

////
////
/////////////////////////////////////////////

#include  "basicintegrator.h"
#include  "integrationformula.h"

using namespace Gascoigne;

/*-----------------------------------------*/

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

  void Rhs(const RightHandSideData& RHS, LocalVector& F, const FemInterface& FEM, const LocalData& Q) const {};
  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector&U, const LocalData& Q) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalData& Q) const;
  void MassMatrix(EntryMatrix& E, const FemInterface& FEM) const { assert(0); }
};

/*-----------------------------------------*/


#endif
