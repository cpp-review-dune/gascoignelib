/*----------------------------   triabase2d.h     ---------------------------*/
#ifndef __triabase2d_H
#define __triabase2d_H
/*----------------------------   triabase2d.h     ---------------------------*/

#include "numfixarray.h"
#include "vertex.h"
#include <array>
#include <cmath>
#include <iostream>

namespace Gascoigne {

// 2d Basis functions on the reference  triangle of size 1x1
/**
 *  M(M+1)/2-1
 *     +
 *  .  | \
 *  .  |   \
 *  M  |     \  2M-2
 *     +------+
 *     0 ... (M-1)
 *
 *  M is the number of DOF's along a line, e.g. M=2 P1, M=3 P2, M=5 P4
 *
 **/
template<int M>
class TriaBase
{
#define NDOFS (M * (M + 1) / 2)
protected:
  mutable numfixarray<NDOFS, double>
    N; // the values of the basis functions in a given point
  mutable std::array<Vertex2d, NDOFS> DN; // first derivatives

public:
  constexpr int n() const { return NDOFS; }

  double psi(int i, double x) const { abort(); }
  double psi_x(int i, double x) const { abort(); }

  // initialize
  void point(const Vertex2d& v) const { abort(); }

  void point_boundary(int ie, const Vertex<1>& v) const { abort(); }

  // access values
  double phi(int i) const
  {
    assert(i < n());
    return N[i];
  }
  double phi_x(int i) const
  {
    assert(i < n());
    return DN[i][0];
  }
  double phi_y(int i) const
  {
    assert(i < n());
    return DN[i][1];
  }
  double phi_z(int i) const
  {
    assert(i < n());
    return DN[i][2];
  }
  const Vertex2d& phi_grad(int i) const
  {
    assert(i < n());
    return DN[i];
  }

  const Vertex2d* normal2d() const { abort(); }
  const Vertex2d* tangent2d() const { abort(); }
  const Vertex3d* normal3d() const { abort(); }
  const Vertex3d* tangent3d() const { abort(); }
#undef NDOFS
};

// P1 - Element
template<>
inline void
TriaBase<2>::point(const Vertex2d& v) const
{
  // Besser mit den Abfragen, da ansonsten der Gradient nicht definiert ist.
  /* assert(v.x()>0); */
  /* assert(v.y()>0); */
  /* assert(v.x()<1); */
  /* assert(v.y()<1); */
  N[0] = 1.0 - v.x() - v.y();
  N[1] = v.x();
  N[2] = v.y();

  DN[0].x() = -1.0;
  DN[0].y() = -1.0;
  DN[1].x() = 1.0;
  DN[1].y() = 0.0;
  DN[2].x() = 0.0;
  DN[2].y() = 1.0;
}
// P2 - Element
template<>
inline void
TriaBase<3>::point(const Vertex2d& v) const
{
  N[0] = 2.0 * v.x() * v.x() + 4.0 * v.x() * v.y() + 2.0 * v.y() * v.y() -
         3.0 * v.x() - 3.0 * v.y() + 1.0;
  N[1] = -4.0 * v.x() * v.x() - 4.0 * v.x() * v.y() + 4.0 * v.x();
  N[2] = 2.0 * v.x() * v.x() - 1.0 * v.x();
  N[3] = -4.0 * v.x() * v.y() - 4.0 * v.y() * v.y() + 4.0 * v.y();
  N[4] = 4.0 * v.x() * v.y();
  N[5] = 2.0 * v.y() * v.y() - 1.0 * v.y();

  DN[0][0] = 4.0 * v.x() + 4.0 * v.y() - 3.0;
  DN[1][0] = -8.0 * v.x() - 4.0 * v.y() + 4.0;
  DN[2][0] = 4.0 * v.x() - 1.0;
  DN[3][0] = -4.0 * v.y();
  DN[4][0] = 4.0 * v.y();
  DN[5][0] = 0.0;

  DN[0][1] = 4.0 * v.x() + 4.0 * v.y() - 3.0;
  DN[1][1] = -4.0 * v.x();
  DN[2][1] = 0.0;
  DN[3][1] = -4.0 * v.x() - 8.0 * v.y() + 4.0;
  DN[4][1] = 4.0 * v.x();
  DN[5][1] = 4.0 * v.y() - 1.0;
}

// 2D base functions on a quad that is split into 2 triangles
/**
 *  Type0   Type1  Type2
 *  +---+   +---+  +---+
 *  |   |   | / |  | \ |
 *  +---+   +---+  +---+
 **/
template<int M>
class TriaQuadBase
{
#define NDOFS (M * M)
protected:
  mutable TriaBase<M> TB; // Basis function on one triangle
  mutable numfixarray<NDOFS, double>
    N; // the values of the basis functions in a given point
  mutable std::array<Vertex2d, NDOFS> DN; // first derivatives

  mutable Vertex2d bn, bt; // normal and tangential on the boundary

public:
  constexpr int n() const { return NDOFS; }
#undef NDOFS

  // initialize
  void point(const Vertex2d& w) const
  {
    Vertex2d v = w;

    int type = 1;

    assert((type == 0) || (type == 1) || (type == 2));
    for (auto& it : N)
      it = 0;
    for (auto& it : DN)
      it.zero();

    if (type == 1) {
      if (v.y() >=
          v.x()) // Eigentlich y>x damit Punkte auf der Linie nicht vorkommen
      {
        v.y() = 1.0 - v.y();
        TB.point(v);
        int i = 0;
        for (int iy = M - 1; iy >= 0; --iy)
          for (int ix = 0; ix <= iy; ++ix, ++i) {
            N[M * iy + ix] = TB.phi(i);
            DN[M * iy + ix].x() = TB.phi_x(i);
            DN[M * iy + ix].y() = -TB.phi_y(i);
          }
        v.y() = 1.0 - v.y();
      } else if (v.x() > v.y()) {
        v.x() = 1.0 - v.x();
        TB.point(v);
        int i = 0;
        for (int iy = 0; iy < M; ++iy)
          for (int ix = M - 1; ix >= iy; --ix, ++i) {
            N[M * iy + ix] = TB.phi(i);
            DN[M * iy + ix].x() = -TB.phi_x(i);
            DN[M * iy + ix].y() = TB.phi_y(i);
          }
        v.x() = 1.0 - v.x();
      } else {
        std::cerr << "TriaQuad::point(v). v on the line" << std::endl;
        abort();
      }
    } else if (type == 2) {
      if (v.x() + v.y() <= 1) // Besser x+y < 1
      {
        TB.point(v);
        int i = 0;
        for (int iy = 0; iy < M; ++iy)
          for (int ix = 0; ix <= (M - 1) - iy; ++ix, ++i) {
            N[M * iy + ix] = TB.phi(i);
            DN[M * iy + ix].x() = TB.phi_x(i);
            DN[M * iy + ix].y() = TB.phi_y(i);
          }
      } else if (v.x() + v.y() > 1) {
        v.x() = 1.0 - v.x();
        v.y() = 1.0 - v.y();
        TB.point(v);
        int i = 0;
        for (int iy = M - 1; iy >= 0; --iy)
          for (int ix = M - 1; ix >= (M - 1) - iy; --ix, ++i) {
            N[M * iy + ix] = TB.phi(i);
            DN[M * iy + ix].x() = -TB.phi_x(i);
            DN[M * iy + ix].y() = -TB.phi_y(i);
          }
        v.x() = 1.0 - v.x();
        v.y() = 1.0 - v.y();
      } else {
        std::cerr << "TriaQuad::point(v). v on the line" << std::endl;
        abort();
      }
    }
  }

  void point_boundary(int ie, const Vertex1d& s1) const
  {
    Vertex2d s;
    if (ie == 0) {
      s.x() = s1.x();
      s.y() = 0.;
      bn.x() = 0.;
      bn.y() = -1.;
      bt.x() = 1.;
      bt.y() = 0.;
    } else if (ie == 1) {
      s.x() = 1.;
      s.y() = s1.x();
      bn.x() = 1.;
      bn.y() = 0.;
      bt.x() = 0.;
      bt.y() = 1.;
    } else if (ie == 2) {
      s.x() = 1 - s1.x();
      s.y() = 1.;
      bn.x() = 0.;
      bn.y() = 1.;
      bt.x() = -1.;
      bt.y() = 0.;
    } else {
      s.x() = 0.;
      s.y() = 1 - s1.x();
      bn.x() = -1.;
      bn.y() = 0.;
      bt.x() = 0.;
      bt.y() = -1.;
    }
    point(s);
  }

  // access values
  double phi(int i) const
  {
    assert(i < n());
    return N[i];
  }
  double phi_x(int i) const
  {
    assert(i < n());
    return DN[i][0];
  }
  double phi_y(int i) const
  {
    assert(i < n());
    return DN[i][1];
  }

  const Vertex2d& phi_grad(int i) const
  {
    assert(i < n());
    return DN[i];
  }

  const Vertex2d* normal2d() const { return &bn; }
  const Vertex2d* tangent2d() const { return &bt; }
  const Vertex3d* normal3d() const { abort(); }
  const Vertex3d* tangent3d() const { abort(); }
};

}

/*----------------------------   triabase2d.h     ---------------------------*/
/* end of #ifndef __triabase2d_H */
#endif
/*----------------------------   triabase2d.h     ---------------------------*/
