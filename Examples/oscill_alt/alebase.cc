#include "alebase.h"

using namespace std;

namespace Gascoigne {

void AleBase::compute_transformation(const Vertex3d&    v,
                                     const FemFunction& U,
                                     const FemFunction& H) const {
  __v= v;
  // F
  __F(2, 0)= H[0].x() * v.z();
  __F(2, 1)= H[0].y() * v.z();
  __F(2, 2)= H[0].m();
  //
  __J= H[0].m();
  assert(__J > 0);

  // J F^{-1}
  __Ftilde(0, 0)= H[0].m();
  __Ftilde(1, 1)= H[0].m();
  __Ftilde(2, 0)= -H[0].x() * v.z();
  __Ftilde(2, 1)= -H[0].y() * v.z();
  __Ftilde(2, 2)= 1.0;

  // nabla V F^{-1}
  __nablaV_Fi.zero();
  for (int i= 0; i < 3; ++i)
    for (int j= 0; j < 2; ++j)
      __nablaV_Fi(i, j)= U[i + 1][j + 1] - v.z() * H[0][j + 1] * U[i + 1].z() / __J;
  for (int i= 0; i < 3; ++i)
    __nablaV_Fi(i, 2)= U[i + 1].z() / __J;
}

//////////////////////////////////////////////////

AleBase::AleBase() {
  __F.resize(3, 3);
  __F.zero();
  __F(0, 0)= 1.0;
  __F(1, 1)= 1.0;

  __Ftilde.resize(3, 3);
  __Ftilde.zero();
  __Ftilde(2, 2)= 1.0;

  __nablaV_Fi.resize(3, 3);
}

// --------------------------------------------------

double AleBase::DV_nablaV_Fi(int                 i,
                             int                 j,
                             int                 d,
                             const FemFunction&  U,
                             const FemFunction&  H,
                             const TestFunction& M) const {
  if (i != d)
    return 0.0;
  if (i < 2)
    return M[j + 1] - __v.z() * H[0][j + 1] * M.z() / __J;
  else
    return M.z() / __J;
}

}  // namespace Gascoigne
