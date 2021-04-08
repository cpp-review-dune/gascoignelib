/*----------------------------   fsi.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __fsi_H
#define __fsi_H
/*----------------------------   fsi.h     ---------------------------*/

#include "../eigen3/Eigen/Dense"
#include "boundaryequation.h"
#include "equation.h"
#include "paramfile.h"

/*-----------------------------------------*/

namespace Gascoigne {
template <int DIM> class EQ : public Equation {

protected:
  double lambda, mu;

public:
  ~EQ() {}
  EQ() { abort(); }
  EQ(const ParamFile *pf);

  std::string GetName() const { return "EQ"; }

  int GetNcomp() const { return DIM; }

  void point(double h, const FemFunction &U, const Vertex<DIM> &v) const {}
  void Form(VectorIterator b, const FemFunction &U,
            const TestFunction &N) const;

  void Matrix(EntryMatrix &A, const FemFunction &U, const TestFunction &M,
              const TestFunction &N) const;
};

} // namespace Gascoigne

/*----------------------------   fsi.h     ---------------------------*/
/* end of #ifndef __fsi_H */
#endif
/*----------------------------   fsi.h     ---------------------------*/
