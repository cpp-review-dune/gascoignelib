/*----------------------------   faceintegrator.h     ---------------------------*/
/*      $Id$                 */
#ifndef __faceintegrator_H
#define __faceintegrator_H
/*----------------------------   faceintegrator.h     ---------------------------*/

#include "faceequation.h"
#include "feminterface.h"
#include "integrationformula.h"

namespace Gascoigne
{
  class FaceIntegratorInterface
    {

    protected:
      IntegrationFormulaInterface*  __IFF;
      
      mutable TestFunction          __NN;
      
      IntegrationFormulaInterface*& FaceFormulaPointer() { return __IFF;}
      const IntegrationFormulaInterface* FaceFormula() const { assert(__IFF); return __IFF; }

      // --------------------------------------------------

      virtual void universal_point_face(const FemInterface& FEM1,const FemInterface& FEM2, FemFunction& U1, FemFunction& U2, const LocalVector& U) const=0;


    public:

      FaceIntegratorInterface() : __IFF(0) {}
      virtual ~FaceIntegratorInterface() { }
      
      
      virtual void FaceForm(const FaceEquation& EQ, LocalVector& F, const FemInterface& FEM1, const FemInterface& FEM2, const LocalVector& U) const=0;
      virtual void FaceMatrix(const FaceEquation& EQ, EntryMatrix& E, const FemInterface& FEM1, const FemInterface& FEM2, const LocalVector& U) const=0;

    };

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  template<int DIM, int Q>
    class FaceIntegrator : public FaceIntegratorInterface
    {
    protected:
      void universal_point_face(const FemInterface& FEM1,const FemInterface& FEM2, FemFunction& U1, FemFunction& U2, const LocalVector& U) const;
      
    public:
      FaceIntegrator();
      ~FaceIntegrator()
	{
	  if (__IFF) { delete __IFF; __IFF=0; }
	}
      
      void FaceForm(const FaceEquation& EQ, LocalVector& F, const FemInterface& FEM1, const FemInterface& FEM2, const LocalVector& U) const;
      void FaceMatrix(const FaceEquation& EQ, EntryMatrix& E, const FemInterface& FEM1, const FemInterface& FEM2, const LocalVector& U) const;
    };
  
  
}




/*----------------------------   faceintegrator.h     ---------------------------*/
/* end of #ifndef __faceintegrator_H */
#endif
/*----------------------------   faceintegrator.h     ---------------------------*/
