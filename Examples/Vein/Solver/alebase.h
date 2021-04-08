/*----------------------------   alebase.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __alebase_H
#define __alebase_H
/*----------------------------   alebase.h     ---------------------------*/

#include "gascoigne.h"
#include "nmatrix.h"
#include "nvector.h"
#include "vertex.h"

namespace Gascoigne {

class AleBase {
public:
  mutable nmatrix<double> __F, __Ftilde, __nablaV_Fi;
  mutable double __J;
  mutable Vertex3d __v;

public:
  AleBase();

  void compute_transformation(const Vertex3d &v, const FemFunction &U,
                              const FemFunction &H) const;

  double DV_nablaV_Fi(int i, int j, int d, const FemFunction &U,
                      const FemFunction &H, const TestFunction &M) const;
};

} // namespace Gascoigne

/*----------------------------   alebase.h     ---------------------------*/
/* end of #ifndef __alebase_H */
#endif
/*----------------------------   alebase.h     ---------------------------*/
