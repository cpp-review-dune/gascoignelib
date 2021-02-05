/*----------------------------   cgbase.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __cgbase_H
#define __cgbase_H
/*----------------------------   cgbase.h     ---------------------------*/

#include "numfixarray.h"
#include "vertex.h"
#include <array>
#include <cmath>
#include <iostream>

namespace Gascoigne {

/**
 *  Lagrange Finite Element in DIM
 *  dimensions and M unknowns in each
 *  direction.
 **/
template <int DIM, int M> class CGBase {
public:
  // access and info
  constexpr int n() const { return (DIM == 2) ? M * M : M * M * M; }

private:
protected:
#define NDOFS ((DIM == 2) ? (M * M) : (M * M * M))
  mutable numfixarray<NDOFS, double>
      N; // the values of the basis functions in a given point
  mutable std::array<Vertex<DIM>, NDOFS> DN; // first derivatives
#undef NDOFS

  std::array<std::array<double, M>, M> alpha; // the 1d coefficients

public:
  CGBase() { abort(); }
  double psi(int i, double x) const { abort(); }
  double psi_x(int i, double x) const { abort(); }

  // initialize
  void point(const Vertex<DIM> &v) const {
    if (DIM == 2) {
      int i = 0;
      for (int iy = 0; iy < M; ++iy) {
        // maybe optimization by first computing psi's?
        for (int ix = 0; ix < M; ++ix, ++i) {
          N[i] = psi(ix, v[0]) * psi(iy, v[1]);
          DN[i][0] = psi_x(ix, v[0]) * psi(iy, v[1]);
          DN[i][1] = psi(ix, v[0]) * psi_x(iy, v[1]);
        }
      }
    } else if (DIM == 3) {
      int i = 0;
      for (int iz = 0; iz < M; ++iz)
        for (int iy = 0; iy < M; ++iy) {
          // maybe optimization by first computing psi's?
          for (int ix = 0; ix < M; ++ix, ++i) {
            N[i] = psi(ix, v[0]) * psi(iy, v[1]) * psi(iz, v[2]);
            DN[i][0] = psi_x(ix, v[0]) * psi(iy, v[1]) * psi(iz, v[2]);
            DN[i][1] = psi(ix, v[0]) * psi_x(iy, v[1]) * psi(iz, v[2]);
            DN[i][2] = psi(ix, v[0]) * psi(iy, v[1]) * psi_x(iz, v[2]);
          }
        }
    } else
      abort();
  }

  void point_boundary(int ie, const Vertex<DIM - 1> &v) const { abort(); }
  const std::array<int, 2> *faces() const { abort(); }

  // access values
  double phi(int i) const {
    assert(i < n());
    return N[i];
  }
  double phi_x(int i) const {
    assert(i < n());
    return DN[i][0];
  }
  double phi_y(int i) const {
    assert(i < n());
    return DN[i][1];
  }
  double phi_z(int i) const {
    assert(i < n());
    return DN[i][2];
  }

  const Vertex<DIM> &phi_grad(int i) const {
    assert(i < n());
    return DN[i];
  }

  const Vertex2d *normal2d() const { abort(); }
  const Vertex2d *tangent2d() const { abort(); }
  const Vertex3d *normal3d() const { abort(); }
  const Vertex3d *tangent3d() const { abort(); }
};

////////////////////////////////////////////////// Constructor

////////////////////////////////////////////////// 2d
template <> inline CGBase<2, 2>::CGBase() {
  alpha[0][0] = -1.0;
  alpha[0][1] = 1.0;
  alpha[1][0] = 1.0;
  alpha[1][1] = 0.0;
}
template <> inline CGBase<2, 3>::CGBase() {
  alpha[0][0] = 2.0;
  alpha[0][1] = -3.0;
  alpha[0][2] = 1.0;
  alpha[1][0] = -4.0;
  alpha[1][1] = 4.0;
  alpha[1][2] = 0.0;
  alpha[2][0] = 2.0;
  alpha[2][1] = -1.0;
  alpha[2][2] = 0.0;
}
template <> inline CGBase<2, 5>::CGBase() {
  // L1:   32./3. *x^4  -80./3. *x^3  +70./3. *x^2 -25./3. *x +1.
  // L2: -128./3. *x^4  +96.    *x^3 -208./3. *x^2 +16.    *x +0.
  // L3:   64.    *x^4 -128.    *x^3  +76.    *x^2 -12.    *x +0.
  // L4: -128./3. *x^4 +224./3. *x^3 -112./3. *x^2 +16./3. *x +0.
  // L5:   32./3. *x^4  -16.    *x^3  +22./3. *x^2 -1.0    *x
  alpha[0][0] = 32. / 3.;
  alpha[0][1] = -80. / 3.;
  alpha[0][2] = +70. / 3.;
  alpha[0][3] = -25. / 3.;
  alpha[0][4] = +1.;
  alpha[1][0] = -128. / 3.;
  alpha[1][1] = +96.;
  alpha[1][2] = -208. / 3.;
  alpha[1][3] = +16.;
  alpha[1][4] = +0.;
  alpha[2][0] = 64.;
  alpha[2][1] = -128.;
  alpha[2][2] = +76.;
  alpha[2][3] = -12.;
  alpha[2][4] = +0.;
  alpha[3][0] = -128. / 3.;
  alpha[3][1] = +224. / 3.;
  alpha[3][2] = -112. / 3.;
  alpha[3][3] = +16. / 3.;
  alpha[3][4] = +0.;
  alpha[4][0] = 32. / 3.;
  alpha[4][1] = -16.;
  alpha[4][2] = +22. / 3.;
  alpha[4][3] = -1.0;
  alpha[4][4] = +0.;
}
// q1
template <> inline double CGBase<2, 2>::psi(int i, double x) const {
  return alpha[i][0] * x + alpha[i][1];
}
template <> inline double CGBase<2, 2>::psi_x(int i, double x) const {
  return alpha[i][0];
}
// q2
template <> inline double CGBase<2, 3>::psi(int i, double x) const {
  return (alpha[i][0] * x + alpha[i][1]) * x + alpha[i][2];
}
template <> inline double CGBase<2, 3>::psi_x(int i, double x) const {
  return 2.0 * alpha[i][0] * x + alpha[i][1];
}
// q4
template <> inline double CGBase<2, 5>::psi(int i, double x) const {
  return alpha[i][4] +
         x * (alpha[i][3] +
              x * (alpha[i][2] + x * (alpha[i][1] + x * alpha[i][0])));
}
template <> inline double CGBase<2, 5>::psi_x(int i, double x) const {
  return alpha[i][3] + x * (2.0 * alpha[i][2] +
                            x * (3.0 * alpha[i][1] + x * 4.0 * alpha[i][0]));
}

////////////////////////////////////////////////// 3d
template <> inline CGBase<3, 2>::CGBase() {
  alpha[0][0] = -1.0;
  alpha[0][1] = 1.0;
  alpha[1][0] = 1.0;
  alpha[1][1] = 0.0;
}
template <> inline CGBase<3, 3>::CGBase() {
  alpha[0][0] = 2.0;
  alpha[0][1] = -3.0;
  alpha[0][2] = 1.0;
  alpha[1][0] = -4.0;
  alpha[1][1] = 4.0;
  alpha[1][2] = 0.0;
  alpha[2][0] = 2.0;
  alpha[2][1] = -1.0;
  alpha[2][2] = 0.0;
}
template <> inline CGBase<3, 5>::CGBase() {
  // L1:   32./3. *x^4  -80./3. *x^3  +70./3. *x^2 -25./3. *x +1.
  // L2: -128./3. *x^4  +96.    *x^3 -208./3. *x^2 +16.    *x +0.
  // L3:   64.    *x^4 -128.    *x^3  +76.    *x^2 -12.    *x +0.
  // L4: -128./3. *x^4 +224./3. *x^3 -112./3. *x^2 +16./3. *x +0.
  // L5:   32./3. *x^4  -16.    *x^3  +22./3. *x^2 -1.0    *x
  alpha[0][0] = 32. / 3.;
  alpha[0][1] = -80. / 3.;
  alpha[0][2] = +70. / 3.;
  alpha[0][3] = -25. / 3.;
  alpha[0][4] = +1.;
  alpha[1][0] = -128. / 3.;
  alpha[1][1] = +96.;
  alpha[1][2] = -208. / 3.;
  alpha[1][3] = +16.;
  alpha[1][4] = +0.;
  alpha[2][0] = 64.;
  alpha[2][1] = -128.;
  alpha[2][2] = +76.;
  alpha[2][3] = -12.;
  alpha[2][4] = +0.;
  alpha[3][0] = -128. / 3.;
  alpha[3][1] = +224. / 3.;
  alpha[3][2] = -112. / 3.;
  alpha[3][3] = +16. / 3.;
  alpha[3][4] = +0.;
  alpha[4][0] = 32. / 3.;
  alpha[4][1] = -16.;
  alpha[4][2] = +22. / 3.;
  alpha[4][3] = -1.0;
  alpha[4][4] = +0.;
}
// q1
template <> inline double CGBase<3, 2>::psi(int i, double x) const {
  return alpha[i][0] * x + alpha[i][1];
}
template <> inline double CGBase<3, 2>::psi_x(int i, double x) const {
  return alpha[i][0];
}
// q2
template <> inline double CGBase<3, 3>::psi(int i, double x) const {
  return (alpha[i][0] * x + alpha[i][1]) * x + alpha[i][2];
}
template <> inline double CGBase<3, 3>::psi_x(int i, double x) const {
  return 2.0 * alpha[i][0] * x + alpha[i][1];
}
// q4
template <> inline double CGBase<3, 5>::psi(int i, double x) const {
  return alpha[i][4] +
         x * (alpha[i][3] +
              x * (alpha[i][2] + x * (alpha[i][1] + x * alpha[i][0])));
}
template <> inline double CGBase<3, 5>::psi_x(int i, double x) const {
  return alpha[i][3] + x * (2.0 * alpha[i][2] +
                            x * (3.0 * alpha[i][1] + x * 4.0 * alpha[i][0]));
}

} // namespace Gascoigne

/*----------------------------   cgbase.h     ---------------------------*/
/* end of #ifndef __cgbase_H */
#endif
/*----------------------------   cgbase.h     ---------------------------*/
