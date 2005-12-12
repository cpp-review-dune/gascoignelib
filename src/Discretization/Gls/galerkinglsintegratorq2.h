#ifndef  __GalerkinGlsIntegratorQ2_h
#define  __GalerkinGlsIntegratorQ2_h


#include  "galerkinintegrator.h"
#include  "glsintegratorq2.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments GalerkinGlsIntegratorQ2

////
////
/////////////////////////////////////////////

template<int DIM>
class GalerkinGlsIntegratorQ2 : public GalerkinIntegrator<DIM>
{
protected:

  GlsIntegratorQ2<DIM>  Gls;

public:


//
////  Con(De)structor 
//

  GalerkinGlsIntegratorQ2<DIM>() : GalerkinIntegrator<DIM>() {}
  ~GalerkinGlsIntegratorQ2<DIM>() {}

  std::string GetName() const {return "GalerkinGlsIntegratorQ2";}

  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector&U, 
      const LocalNodeData& Q, const LocalCellData& QC) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
      const LocalNodeData& Q, const LocalCellData& QC) const;

};
}

#endif
