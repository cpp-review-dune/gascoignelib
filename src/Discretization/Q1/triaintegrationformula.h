/*----------------------------   triaintegrationformulas.h
 * ---------------------------*/
/*      $Id:$                 */
#ifndef __triaintegrationformulas_H
#define __triaintegrationformulas_H
/*----------------------------   triaintegrationformulas.h
 * ---------------------------*/

#include "integrationformula.h"
#include "integrationformulainterface.h"

namespace Gascoigne {

class TriaOnepointFormula : public IntegrationFormulaInterface
{
public:
  int n() const { return 1; }
  double w(int k) const { return 0.5; }
  void xi(Vertex2d& v, int k) const
  {
    assert(k == 0);
    v.x() = 1.0 / 3;
    v.y() = 1.0 / 3;
  }
};

class TriaSixpointFormula : public IntegrationFormulaInterface
{
public:
  int n() const { return 6; }
  double w(int k) const { return 0.5 / 6.0; }
  void xi(Vertex2d& v, int k) const
  {
    if (k == 0) {
      v.x() = 0.659027622374092;
      v.y() = 0.231933368553031;
    } else if (k == 1) {
      v.x() = 0.659027622374092;
      v.y() = 0.109039009072877;
    } else if (k == 2) {
      v.x() = 0.231933368553031;
      v.y() = 0.659027622374092;
    } else if (k == 3) {
      v.x() = 0.231933368553031;
      v.y() = 0.109039009072877;
    } else if (k == 4) {
      v.x() = 0.109039009072877;
      v.y() = 0.659027622374092;
    } else if (k == 5) {
      v.x() = 0.109039009072877;
      v.y() = 0.231933368553031;
    } else
      abort();
  }
};

class TriaTwelvepointFormula : public IntegrationFormulaInterface
{
public:
  int n() const { return 12; }
  double w(int k) const
  {

    if (k == 0)
      return 0.116786275726379 / 2.0;
    else if (k == 1)
      return 0.116786275726379 / 2.0;
    else if (k == 2)
      return 0.116786275726379 / 2.0;
    else if (k == 3)
      return 0.050844906370207 / 2.0;
    else if (k == 4)
      return 0.050844906370207 / 2.0;
    else if (k == 5)
      return 0.050844906370207 / 2.0;
    else if (k == 6)
      return 0.082851075618374 / 2.0;
    else if (k == 7)
      return 0.082851075618374 / 2.0;
    else if (k == 8)
      return 0.082851075618374 / 2.0;
    else if (k == 9)
      return 0.082851075618374 / 2.0;
    else if (k == 10)
      return 0.082851075618374 / 2.0;
    else if (k == 11)
      return 0.082851075618374 / 2.0;

    else
      abort();
  }
  void xi(Vertex2d& v, int k) const
  {
    if (k == 0) {
      v.x() = 0.249286745170910;
      v.y() = 0.249286745170910;
    } else if (k == 1) {
      v.x() = 0.249286745170910;
      v.y() = 0.501426509658179;
    } else if (k == 2) {
      v.x() = 0.501426509658179;
      v.y() = 0.249286745170910;
    } else if (k == 3) {
      v.x() = 0.063089014491502;
      v.y() = 0.063089014491502;
    } else if (k == 4) {
      v.x() = 0.063089014491502;
      v.y() = 0.873821971016996;
    } else if (k == 5) {
      v.x() = 0.873821971016996;
      v.y() = 0.063089014491502;
    } else if (k == 6) {
      v.x() = 0.310352451033784;
      v.y() = 0.636502499121399;
    } else if (k == 7) {
      v.x() = 0.636502499121399;
      v.y() = 0.053145049844817;
    } else if (k == 8) {
      v.x() = 0.053145049844817;
      v.y() = 0.310352451033784;
    } else if (k == 9) {
      v.x() = 0.636502499121399;
      v.y() = 0.310352451033784;
    } else if (k == 10) {
      v.x() = 0.310352451033784;
      v.y() = 0.053145049844817;
    } else if (k == 11) {
      v.x() = 0.053145049844817;
      v.y() = 0.636502499121399;
    } else
      abort();
  }
};

// Integrationformula of 2 tria-formulas
template<class TRIAFORMULA>
class TriaQuadFormula : public IntegrationFormulaInterface
{
  TRIAFORMULA triaformula;

  /**
   *  TYPE0 TYPE1  TYPE2
   *  +---+ +---+  +---+
   *  |   | | / |  | \ |
   *  +---+ +---+  +---+
   *  TYPE0 wird wie TYPE2 verwendet
   **/

public:
  int n() const { return 2 * triaformula.n(); } //
  double w(int k) const { return triaformula.w(k % triaformula.n()); }

  void xi(Vertex1d& v, int k) const { assert(0); }
  void xi(Vertex3d& v, int k) const { assert(0); }

  void xi(Vertex2d& v, int k) const
  {
    int type = 1;

    type = type % 2; // damit TYPE 2 und TYPE 0 gleich behandelt werden
    assert(k < n());
    int intria = k / triaformula.n();
    int kintria = k % triaformula.n();
    triaformula.xi(v, kintria);

    if (type == 1) {
      if (intria == 0) // upper left
        v.y() = 1.0 - v.y();
      else // lower right
        v.x() = 1.0 - v.x();
    } else if (type == 0) // 0 = 2 wird gleich behandelt
    {
      if (intria == 1) // upper right
      {
        v.y() = 1.0 - v.y();
        v.x() = 1.0 - v.x();
      }
    } else
      abort();
  }
};

}

/*----------------------------   triaintegrationformulas.h
 * ---------------------------*/
/* end of #ifndef __triaintegrationformulas_H */
#endif
/*----------------------------   triaintegrationformulas.h
 * ---------------------------*/
