/*----------------------------    boundarysolid.h ---------------------------*/
/*      $Id:$                 */
#ifndef __Boundarysolid_H
#define __Boundarysolid_H
/*----------------------------    boundarysolid.h ---------------------------*/

#include "boundaryequation.h"
#include "eigen3/Eigen/Dense"
#include "equation.h"
#include "lpsequation.h"
#include "paramfile.h"

#include <array>

/*-----------------------------------------*/

namespace Gascoigne {
template<int DIM>
class BoundarySolid : public BoundaryEquation
{

protected:
  typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
  typedef Eigen::Matrix<double, DIM, 1> VECTOR;

  mutable VECTOR __n;
  mutable MATRIX F;
  mutable double J, __pressure;

  // solid
  mutable VECTOR g_dU[DIM];

  // mutable FemFunction *OLD, *DEF, *DEFOLD;

  void SetFemData(FemData& q) const
  {
    // assert(q.find("OLD")!=q.end());
    // OLD = &q["OLD"];

    // assert(q.find("DEF")!=q.end());
    // DEF = &q["DEF"];

    // assert(q.find("DEFOLD")!=q.end());
    // DEFOLD = &q["DEFOLD"];
  }

public:
  ~BoundarySolid() {}
  BoundarySolid() { abort(); }
  BoundarySolid(const ParamFile* pf);

  std::string GetName() const { return "BoundarySolid"; }

  int GetNcomp() const { return DIM; }
  void Form(VectorIterator b,
            const FemFunction& U,
            const TestFunction& N,
            int col) const;
  void Matrix(EntryMatrix& E,
              const FemFunction& U,
              const TestFunction& M,
              const TestFunction& N,
              int col) const;
  void pointboundary(double h,
                     const FemFunction& U,
                     const Vertex<DIM>& v,
                     const Vertex<DIM>& n) const;
  void point_M(int j, const FemFunction& U, const TestFunction& M) const;
};

} // namespace Gascoigne

/*----------------------------   boundarysolid.h ---------------------------*/
/* end of #ifndef __ boundarysolid_H */
#endif
/*----------------------------    boundarysolid.h ---------------------------*/
