#ifndef  __GalerkinLpsIntegratorQ2_h
#define  __GalerkinLpsIntegratorQ2_h


/////////////////////////////////////////////
////
////@brief
////  ... comments GalerkinLpsIntegratorQ2

////
////
/////////////////////////////////////////////

#include  "galerkinintegrator.h"
#include  "lpsintegrator.h"

namespace Gascoigne
{
template<int DIM>
class GalerkinLpsIntegratorQ2 : virtual public GalerkinIntegrator<DIM>
{
protected:

  LpsIntegratorQ2<DIM>  Lps;

public:

//
////  Con(De)structor 
//

  GalerkinLpsIntegratorQ2<DIM>() {
    if (DIM==2)
      {
	FormFormulaPointer() = new QuadGauss9;
	ErrorFormulaPointer() = new QuadGauss16;
	BoundaryFormulaPointer() = new LineGauss3;
      }
    else if (DIM==3)
      {
	FormFormulaPointer() = new HexGauss27;
	ErrorFormulaPointer() = new HexGauss64;
	BoundaryFormulaPointer() = new QuadGauss9;
      }
    assert(FormFormulaPointer());
    assert(ErrorFormulaPointer());
    assert(BoundaryFormulaPointer());
  }
  
  ~GalerkinLpsIntegratorQ2<DIM>() {}

  std::string GetName() const {return "GalerkinLpsIntegratorQ2";}

  void Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector&U, const LocalNodeData& Q) const;
  void Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, const LocalNodeData& Q) const;
};
}

#endif
