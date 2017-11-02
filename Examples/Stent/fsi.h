/*----------------------------   fsi.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __fsi_H
#define __fsi_H
/*----------------------------   fsi.h     ---------------------------*/



#include  "equation.h"
#include  "boundaryequation.h"
#include  "paramfile.h"
#include  "lpsequation.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  class FSI : public BoundaryEquation, public LpsEquation
    {
      
    protected:
      double rho_f, nu_f;
      double P2min,P2max,P2period,P4,P7;
      double lps0;

      mutable FemFunction *OLD;
      mutable double P2;
      
      mutable double lps;
      mutable Vertex3d normal;
      
    public:

      void SetFemData(FemData& q) const
      {
	assert(q.find("OLD")!=q.end());
	OLD = &q["OLD"];
      }

    
      ~FSI() { }
      FSI() { abort(); }
      FSI(const ParamFile* pf);

      void point(double h, const FemFunction& U, const Vertex<3>& v) const;
      std::string GetName() const { return "FSI"; }
      
      int    GetNcomp  () const { return 4; }


      void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
      void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
      


      void Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const;
      void Matrix(EntryMatrix& E, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const;
      void pointboundary(double h, const FemFunction& U, const Vertex3d& v, const Vertex3d& n) const;



      void lpspoint(double h, const FemFunction& U, const Vertex<3>& v) const;
      void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
        
      void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const;
      

        
  };
  
  
}

/*----------------------------   fsi.h     ---------------------------*/
/* end of #ifndef __fsi_H */
#endif
/*----------------------------   fsi.h     ---------------------------*/
