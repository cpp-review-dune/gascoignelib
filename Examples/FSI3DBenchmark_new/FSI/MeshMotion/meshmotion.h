/*----------------------------   meshmotion.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __meshmotion_H
#define __meshmotion_H
/*----------------------------   meshmotion.h     ---------------------------*/

#include "boundaryequation.h"
#include "eigen3/Eigen/Dense"
#include "equation.h"
#include "lpsequation.h"
#include "paramfile.h"

#include <array>

/*-----------------------------------------*/

namespace Gascoigne {
template <int DIM>
class MeshMotion : public LpsEquation // , public BoundaryEquation
{

protected:
  typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
  typedef Eigen::Matrix<double, DIM, 1> VECTOR;

  mutable double __h, ext;
  mutable Vertex<DIM> __v;

  mutable int domain;

  mutable FemFunction *UOLD_Vec, *U_Vec;

  void SetFemData(FemData &q) const {
    assert(q.find("U_Vec") != q.end());
    U_Vec = &q["U_Vec"];

    assert(q.find("UOLD_Vec") != q.end());
    UOLD_Vec = &q["UOLD_Vec"];
  }

public:
  ~MeshMotion() {}
  MeshMotion() { abort(); }
  MeshMotion(const ParamFile *pf);

  std::string GetName() const { return "MeshMotion"; }

  int GetNcomp() const { return DIM; }

  void point(double h, const FemFunction &U, const Vertex<DIM> &v) const;
  void point_M(int j, const FemFunction &U, const TestFunction &M) const;
  void point_cell(int material) const;

  void Form(VectorIterator b, const FemFunction &U,
            const TestFunction &N) const;

  void Matrix(EntryMatrix &A, const FemFunction &U, const TestFunction &M,
              const TestFunction &N) const;
  void MatrixBlock(EntryMatrix &A, const FemFunction &U,
                   const FemFunction &NNN) const;

  ////////////////////////////////////////////////// LPS

  void lpspoint(double h, const FemFunction &U, const Vertex<DIM> &v) const;

  void StabForm(VectorIterator b, const FemFunction &U, const FemFunction &UP,
                const TestFunction &N) const;

  void StabMatrix(EntryMatrix &A, const FemFunction &U, const TestFunction &Np,
                  const TestFunction &Mp) const;
};

} // namespace Gascoigne

/*----------------------------   meshmotion.h     ---------------------------*/
/* end of #ifndef __meshmotion_H */
#endif
/*----------------------------   meshmotion.h     ---------------------------*/
