#ifndef  __LpsIntegrator_h
#define  __LpsIntegrator_h


/////////////////////////////////////////////
////
////@brief
////  ... comments LpsIntegratorQ1

////
////
/////////////////////////////////////////////

#include  "basicintegrator.h"
#include  "integrationformula.h"

namespace Gascoigne
{

template<int DIM>
class LpsIntegrator : public BasicIntegrator
{
protected:

  mutable Gascoigne::FemFunction   UHP;
  mutable Gascoigne::FemFunction   NLPS, MLPS;
  void Projection(const FemInterface& FEM) const;
  double  CellWeight;

  IntegrationFormulaInterface* _IF;

  const IntegrationFormulaInterface& FormFormula() const { return *_IF;}
  double Volume2MeshSize(double vol) const { return pow(vol,1./float(DIM));}

public:

//
////  Con(De)structor 
//

  LpsIntegrator<DIM>() : CellWeight(1.) {};
  ~LpsIntegrator<DIM>() {}

  std::string GetName() const {return "Lps";}

  void Rhs(const RightHandSideData& RHS, LocalVector& F, const FemInterface& FEM, const LocalNodeData& Q) const {};
  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector&U, const LocalNodeData& Q) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const;
};

/*-----------------------------------------*/

template<int DIM>
class LpsIntegratorQ1 : public LpsIntegrator<DIM>
{
protected:

  std::string GetName() const {return "LpsQ1";}
  void ResetMatrix(EntryMatrix& E, int n, int ncomp) const;
  void VectorReinit(Gascoigne::LocalVector& F, int n, int ncomp) const;

public:

  LpsIntegratorQ1<DIM>();
  ~LpsIntegratorQ1<DIM>() {}

  void Form(const Equation& EQ, Gascoigne::LocalVector& F, const FemInterface& FEM, const LocalVector&U, const LocalNodeData& Q) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const;
};

/*-----------------------------------------*/

template<int DIM>
class LpsIntegratorQ2 : public LpsIntegrator<DIM>
{
protected:

  std::string GetName() const {return "LpsQ2";}

public:

  LpsIntegratorQ2<DIM>();
  ~LpsIntegratorQ2<DIM>() {}
};

/*-----------------------------------------*/

}

#endif
