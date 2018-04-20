/*----------------------------    boundaryfsi.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __BoundaryFSI_H
#define __BoundaryFSI_H
/*----------------------------    boundaryfsi.h     ---------------------------*/



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
    class BoundaryFSI :  public BoundaryEquation
    {
      
    protected:
		typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
   	typedef Eigen::Matrix<double, DIM, 1> VECTOR;

     mutable VECTOR g,g_OLD;
     mutable VECTOR __n,PHI;
     mutable  MATRIX  F,NV, NU,F_OLD,NV_OLD, NU_OLD,NPHI;
     mutable double J,J_OLD ;
     double __nu_f,__rho_f;
     double p_2,p_4;
     mutable FemFunction *OLD, *DEF, *DEFOLD;
     
     void SetFemData(FemData& q) const
      {
				assert(q.find("OLD")!=q.end());
				OLD = &q["OLD"];

				assert(q.find("DEF")!=q.end());
				DEF = &q["DEF"];

				assert(q.find("DEFOLD")!=q.end());
				DEFOLD = &q["DEFOLD"];
      } 
    public:
      ~BoundaryFSI() { }
      BoundaryFSI() {abort(); }
      BoundaryFSI(const ParamFile* pf);
      
      
      std::string GetName() const { return "BoundaryFSI"; }
      
      int    GetNcomp  () const { return DIM+1; }
      


      void Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const; 
      void Matrix(EntryMatrix& E, const FemFunction& U, const TestFunction& M, const TestFunction& N, int col) const; 

      void pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v, const Vertex<DIM>& n) const; 


       
    };
  
  
}

/*----------------------------   boundaryfsi.h     ---------------------------*/
/* end of #ifndef __ boundaryfsi_H */
#endif
/*----------------------------    boundaryfsi.h     ---------------------------*/
