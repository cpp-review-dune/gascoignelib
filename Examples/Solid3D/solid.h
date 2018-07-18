/*----------------------------   solid.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __solid_H
#define __solid_H
/*----------------------------   solid.h     ---------------------------*/



#include  "equation.h"
#include  "boundaryequation.h"
#include  "paramfile.h"
#include  "lpsequation.h"
#include  "eigen3/Eigen/Dense"

#include <array>

/*-----------------------------------------*/

namespace Gascoigne
{
  template<int DIM>
    class Solid : public Equation // , public BoundaryEquation
    {
      
    protected:
      typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
      typedef Eigen::Matrix<double, DIM, 1>   VECTOR;

      // problem parameter
      double rho_s, lambda_s, mu_s,kapa_s;
      std::string mat_law;
  
      // solid
      mutable VECTOR SIGMA_dF[DIM];
      mutable MATRIX SIGMA_dU[DIM];
      
      // stuff from point
      mutable double J;
      mutable MATRIX E, SIGMAs,F,C;
      

      //mutable FemFunction *OLD, *DEF, *DEFOLD;
      
      void SetFemData(FemData& q) const
      {
			/*	assert(q.find("OLD")!=q.end());
				OLD = &q["OLD"];

				assert(q.find("DEF")!=q.end());
				DEF = &q["DEF"];

				assert(q.find("DEFOLD")!=q.end());
				DEFOLD = &q["DEFOLD"];
				*/
      }
      
    public:
      ~Solid() { }
      Solid() { abort(); }
      Solid(const ParamFile* pf);
      
      
      std::string GetName() const { return "Solid"; }
      
      int    GetNcomp  () const { return DIM; }
      
      void point(double h, const FemFunction& U, const Vertex<DIM>& v) const;
      void point_M(int j, const FemFunction& U, const TestFunction& M) const;

      
      void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
            
      void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
      void MatrixBlock(EntryMatrix& A, const FemFunction& U,  const FemFunction& NNN) const;
      
     
    };
  
  
}

/*----------------------------   solid.h     ---------------------------*/
/* end of #ifndef __solid_H */
#endif
/*----------------------------   solid.h     ---------------------------*/
