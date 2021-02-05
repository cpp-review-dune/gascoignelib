/*----------------------------   EQ_MAIN.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __EQ_MAIN_H
#define __EQ_MAIN_H
/*----------------------------   EQ_MAIN.h     ---------------------------*/

#include "boundaryequation.h"
#include "eigen3/Eigen/Dense"
#include "equation.h"
#include "lpsequation.h"
#include "paramfile.h"

#include <array>

/*-----------------------------------------*/

namespace Gascoigne {
template <int DIM>
class FSI_main : public LpsEquation // , public BoundaryEquation
{

protected:
  typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
  typedef Eigen::Matrix<double, DIM, 1> VECTOR;

public:
  ~FSI_main() {}
  FSI_main() { abort(); }
  FSI_main(const ParamFile *pf) {}

  std::string GetName() const { return "FSI_main"; }

  int GetNcomp() const { return DIM + DIM + 1; }

  void point(double h, const FemFunction &U, const Vertex<DIM> &v) const {}
  void point_M(int j, const FemFunction &U, const TestFunction &M) const {}
  void point_cell(int material) const {}

  void Form(VectorIterator b, const FemFunction &U,
            const TestFunction &N) const {}

  void Matrix(EntryMatrix &A, const FemFunction &U, const TestFunction &M,
              const TestFunction &N) const {}
  void MatrixBlock(EntryMatrix &A, const FemFunction &U,
                   const FemFunction &NNN) const {}

  ////////////////////////////////////////////////// LPS

  void lpspoint(double h, const FemFunction &U, const Vertex<DIM> &v) const {}

  void StabForm(VectorIterator b, const FemFunction &U, const FemFunction &UP,
                const TestFunction &N) const {}

  void StabMatrix(EntryMatrix &A, const FemFunction &U, const TestFunction &Np,
                  const TestFunction &Mp) const {}
};

} // namespace Gascoigne

/*----------------------------   EQ_MAIN.h     ---------------------------*/
/* end of #ifndef __EQ_MAIN_H */
#endif
/*----------------------------   EQ_MAIN.h     ---------------------------*/
