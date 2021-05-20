/*----------------------------    boundaryfsi.h ---------------------------*/
/*      $Id:$                 */
#ifndef __boundaryfsi_H
#define __boundaryfsi_H
/*----------------------------    boundaryfsi.h ---------------------------*/

#include "boundaryequation.h"
#include "eigen3/Eigen/Dense"
#include "equation.h"
#include "lpsequation.h"
#include "paramfile.h"

/*-----------------------------------------*/

namespace Gascoigne {

class BoundaryFSI : public BoundaryEquation
{

protected:
  mutable Vertex3d __n3d;

public:
  ~BoundaryFSI() {}
  BoundaryFSI() {}
  BoundaryFSI(const ParamFile* pf) {}

  std::string GetName() const { return "BoundaryFSI"; }

  int GetNcomp() const { return 3 + 1; }

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
                     const Vertex3d& v,
                     const Vertex3d& n) const;
};

} // namespace Gascoigne

/*----------------------------   boundaryfsi.h     ---------------------------*/
/* end of #ifndef __ boundaryfsi_H */
#endif
/*----------------------------    boundaryfsi.h ---------------------------*/
