/*----------------------------   def_solid.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __def_solid_H
#define __def_solid_H
/*----------------------------   def_solid.h     ---------------------------*/



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
    class Def_Solid : public LpsEquation // , public BoundaryEquation
    {
      
    protected:
      typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
      typedef Eigen::Matrix<double, DIM, 1>   VECTOR;

 
      mutable double      __h;
      mutable Vertex<DIM> __v;
      
      mutable int domain;
      
      mutable FemFunction *UOLD_Vec, *U_Vec;
      
      void SetFemData(FemData& q) const
      {
		assert(q.find("U_Vec")!=q.end());
		U_Vec = &q["U_Vec"];

		assert(q.find("UOLD_Vec")!=q.end());
		UOLD_Vec = &q["UOLD_Vec"];
      }
      
    public:
      ~Def_Solid() { }
      Def_Solid() { abort(); }
      Def_Solid(const ParamFile* pf);
      
      
      std::string GetName() const { return "Def_Solid"; }
      
      int  GetNcomp  () const { return DIM; }
      
      void point(double h, const FemFunction& U, const Vertex<DIM>& v) const;
      void point_M(int j, const FemFunction& U, const TestFunction& M) const;
      void point_cell(int material) const ;
      
      void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
            
      void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
      void MatrixBlock(EntryMatrix& A, const FemFunction& U,  const FemFunction& NNN) const;
      
      ////////////////////////////////////////////////// LPS
      
      void lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const;
      
      void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
      
      void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const;
       
    };
  
  
}

/*----------------------------   def_solid.h     ---------------------------*/
/* end of #ifndef __def_solid_H */
#endif
/*----------------------------   def_solid.h     ---------------------------*/
