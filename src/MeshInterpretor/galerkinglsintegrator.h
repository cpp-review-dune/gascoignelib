#ifndef  __GalerkinGlsIntegrator_h
#define  __GalerkinGlsIntegrator_h

/////////////////////////////////////////////
////
////@brief
////  ... comments GlsIntegrator

////
////
/////////////////////////////////////////////

#include  "galerkinintegrator.h"
#include  "glsintegrator.h"

using namespace Gascoigne;

/*-----------------------------------------*/

template<int DIM>
class GalerkinGlsIntegrator : public GalerkinIntegrator<DIM>
{
 protected:

  GlsIntegrator<DIM>  Gls;

 public:


  GalerkinGlsIntegrator<DIM>() : GalerkinIntegrator<DIM>() {};
  ~GalerkinGlsIntegrator<DIM>() {}

  std::string GetName() const {return "GalerkinGlsIntegrator";}

  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector&U, const LocalData& Q) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalData& Q) const;
};

/*-----------------------------------------*/


#endif
