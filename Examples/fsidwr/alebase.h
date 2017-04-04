/*----------------------------   alebase.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __alebase_H
#define __alebase_H
/*----------------------------   alebase.h     ---------------------------*/

#include "nvector.h"
#include "nmatrix.h"
#include "gascoigne.h"


namespace Gascoigne
{
  
  template<int DIM>
    class AleBase
    {
    protected:
      mutable nmatrix<double> __F,__Ftilde,__nablaV_Ftilde,__E;
      mutable double          __J,__trace_E;
      
      
    public:

      AleBase();

      //
      // domain: -1 fluid, 1 solid
      void compute_transformation(const FemFunction& U, int domain) const;


      double DU_trace_E(int d, const FemFunction& U, const TestFunction& M) const ;
     
       double DU_Ftilde(int i, int j, int d,
		       const FemFunction& U, const TestFunction& M) const;

       double DU_nablaV_Ftilde(int i, int j, int d,
			       const FemFunction& U, const TestFunction& M) const;
       double DV_nablaV_Ftilde(int i, int j, int d,
			       const FemFunction& U, const TestFunction& M) const;
	     
      double DU_J(int k, const FemFunction& U, const TestFunction& M) const;
      double DU_F(int i, int j, int d, const FemFunction& U, const TestFunction& M) const;
      double DU_E(int i, int j, int d, const FemFunction& U, const TestFunction& M) const;
 
    };

}




/*----------------------------   alebase.h     ---------------------------*/
/* end of #ifndef __alebase_H */
#endif
/*----------------------------   alebase.h     ---------------------------*/
