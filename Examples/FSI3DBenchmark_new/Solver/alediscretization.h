/*----------------------------   alediscretization.h
 * ---------------------------*/
/*      $Id: alediscretization.h,v 1.16 2010/09/02 09:14:44 richter Exp $ */
#ifndef __alediscretization_H
#define __alediscretization_H
/*----------------------------   alediscretization.h
 * ---------------------------*/

#include "alebasediscretization.h"
#include "baseq22d.h"
#include "dwrfem.h"
#include "galerkinintegratorq2.h"
#include "integratorq1q2.h"
#include "q2lps2d.h"

namespace Gascoigne {
class AleQ12d : public AleDiscretization<Q12d> {
public:
  std::string GetName() const { return "Q1 Ale 2d"; }
};

class AleQ1Lps2d : public AleDiscretization<Q1Lps2d> {
public:
  std::string GetName() const { return "Q1 Ale 2d Lps"; }
  AleQ1Lps2d() : Q1Lps2d() {}
};

class AleQ2Lps2d : public AleDiscretization<Q2Lps2d> {
public:
  std::string GetName() const { return "Q2 Ale 2d Lps"; }
  AleQ2Lps2d() : Q2Lps2d() {}
};
class AleQ2Lps3d : public AleDiscretization<Q2Lps3d> {
public:
  std::string GetName() const { return "Q2 Ale 3d Lps"; }
  AleQ2Lps3d() : Q2Lps3d() {}
};

class AleQ1Lps3d : public AleDiscretization<Q1Lps3d> {
public:
  std::string GetName() const { return "Q1 Ale 3d Lps"; }

  AleQ1Lps3d() : Q1Lps3d() {}
};

class AleDwrQ1Q22d : public AleDiscretization<DwrFemQ1Q22d> {
public:
  std::string GetName() const { return "DwrFem Q1 Q2 Ale 2d"; }

  void AdjointForm(GlobalVector &f, const GlobalVector &u, const Equation &EQ,
                   double d) const {
    abort();
  }

  /* ----------------------------------------- */

  void BoundaryForm(GlobalVector &f, const GlobalVector &u,
                    const IntSet &Colors, const BoundaryEquation &BE,
                    double d) const {
    nmatrix<double> TH, TL;

    GlobalToGlobalData();
    BE.SetParameterData(__QP);

    const IntegratorQ1Q2<2> *I =
        dynamic_cast<const IntegratorQ1Q2<2> *>(GetIntegrator());
    assert(I);

    const FemInterface &HighOrderFem(*GetFem());

    for (IntSet::const_iterator p = Colors.begin(); p != Colors.end(); p++) {
      int col = *p;

      HASHSET<int> habschon;

      const IntVector &q = *GetMesh()->PatchOnBoundary(col);
      const IntVector &l = *GetMesh()->LocalPatchOnBoundary(col);
      for (int i = 0; i < q.size(); i++) {
        int ip = q[i];
        int ile = l[i];

        Transformation(TH, ip);
        TransformationQ1(TL, ip);

        HighOrderFem.ReInit(TH);
        LowOrderFem.ReInit(TL);

        PatchDiscretization::GlobalToLocal(__U, u, ip);
        if (__ADJOINT) {
          if (FluidInterfaceCell(ip))
            for (int i = 0; i < __DEL_F.size(); ++i)
              DeleteTestFunctionsVector(__U, FluidInterfaceNodes(ip),
                                        __DEL_F[i]);
          if (SolidInterfaceCell(ip))
            for (int i = 0; i < __DEL_S.size(); ++i)
              DeleteTestFunctionsVector(__U, SolidInterfaceNodes(ip),
                                        __DEL_S[i]);
        }

        I->BoundaryForm(BE, __F, HighOrderFem, LowOrderFem, __U, ile, col, __QN,
                        __QC);

        if (!__ADJOINT) {
          if (FluidInterfaceCell(ip))
            for (int i = 0; i < __DEL_F.size(); ++i)
              DeleteTestFunctionsVector(__F, FluidInterfaceNodes(ip),
                                        __DEL_F[i]);
          if (SolidInterfaceCell(ip))
            for (int i = 0; i < __DEL_S.size(); ++i)
              DeleteTestFunctionsVector(__F, SolidInterfaceNodes(ip),
                                        __DEL_S[i]);
        }
        PatchDiscretization::LocalToGlobal(f, __F, ip, d);
      }
    }
  }

  void Form(GlobalVector &f, const GlobalVector &u, const Equation &EQ,
            double d) const {
    nmatrix<double> TH, TL;

    GlobalToGlobalData();
    EQ.SetParameterData(__QP);

    const IntegratorQ1Q2<2> *I =
        dynamic_cast<const IntegratorQ1Q2<2> *>(GetIntegrator());
    assert(I);

    const FemInterface &HighOrderFem(*GetFem());

    for (int iq = 0; iq < GetPatchMesh()->npatches(); ++iq) {
      Transformation(TH, iq);
      TransformationQ1(TL, iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem.ReInit(TL);

      PatchDiscretization::GlobalToLocal(__U, u, iq);
      if (__ADJOINT) {
        if (FluidInterfaceCell(iq))
          for (int i = 0; i < __DEL_F.size(); ++i)
            DeleteTestFunctionsVector(__U, FluidInterfaceNodes(iq), __DEL_F[i]);
        if (SolidInterfaceCell(iq))
          for (int i = 0; i < __DEL_S.size(); ++i)
            DeleteTestFunctionsVector(__U, SolidInterfaceNodes(iq), __DEL_S[i]);
      }

      I->Form(EQ, __F, HighOrderFem, LowOrderFem, __U, __QN, __QC);

      if (!__ADJOINT) {
        if (FluidInterfaceCell(iq))
          for (int i = 0; i < __DEL_F.size(); ++i)
            DeleteTestFunctionsVector(__F, FluidInterfaceNodes(iq), __DEL_F[i]);
        if (SolidInterfaceCell(iq))
          for (int i = 0; i < __DEL_S.size(); ++i)
            DeleteTestFunctionsVector(__F, SolidInterfaceNodes(iq), __DEL_S[i]);
      }
      PatchDiscretization::LocalToGlobal(f, __F, iq, d);
    }
  }
};

class AleDwrQ1Q23d : public AleDiscretization<DwrFemQ1Q23d> {
public:
  std::string GetName() const { return "DwrFem Q1 Q2 Ale 3d"; }

  void AdjointForm(GlobalVector &f, const GlobalVector &u, const Equation &EQ,
                   double d) const {
    abort();
  }

  /* ----------------------------------------- */

  void BoundaryForm(GlobalVector &f, const GlobalVector &u,
                    const IntSet &Colors, const BoundaryEquation &BE,
                    double d) const {
    nmatrix<double> TH, TL;

    GlobalToGlobalData();
    BE.SetParameterData(__QP);

    const IntegratorQ1Q2<3> *I =
        dynamic_cast<const IntegratorQ1Q2<3> *>(GetIntegrator());
    assert(I);

    const FemInterface &HighOrderFem(*GetFem());

    for (IntSet::const_iterator p = Colors.begin(); p != Colors.end(); p++) {
      int col = *p;

      HASHSET<int> habschon;

      const IntVector &q = *GetMesh()->PatchOnBoundary(col);
      const IntVector &l = *GetMesh()->LocalPatchOnBoundary(col);
      for (int i = 0; i < q.size(); i++) {
        int ip = q[i];
        int ile = l[i];

        Transformation(TH, ip);
        TransformationQ1(TL, ip);

        HighOrderFem.ReInit(TH);
        LowOrderFem.ReInit(TL);

        PatchDiscretization::GlobalToLocal(__U, u, ip);
        if (__ADJOINT) {
          if (FluidInterfaceCell(ip))
            for (int i = 0; i < __DEL_F.size(); ++i)
              DeleteTestFunctionsVector(__U, FluidInterfaceNodes(ip),
                                        __DEL_F[i]);
          if (SolidInterfaceCell(ip))
            for (int i = 0; i < __DEL_S.size(); ++i)
              DeleteTestFunctionsVector(__U, SolidInterfaceNodes(ip),
                                        __DEL_S[i]);
        }

        I->BoundaryForm(BE, __F, HighOrderFem, LowOrderFem, __U, ile, col, __QN,
                        __QC);

        if (!__ADJOINT) {
          if (FluidInterfaceCell(ip))
            for (int i = 0; i < __DEL_F.size(); ++i)
              DeleteTestFunctionsVector(__F, FluidInterfaceNodes(ip),
                                        __DEL_F[i]);
          if (SolidInterfaceCell(ip))
            for (int i = 0; i < __DEL_S.size(); ++i)
              DeleteTestFunctionsVector(__F, SolidInterfaceNodes(ip),
                                        __DEL_S[i]);
        }
        PatchDiscretization::LocalToGlobal(f, __F, ip, d);
      }
    }
  }

  void Form(GlobalVector &f, const GlobalVector &u, const Equation &EQ,
            double d) const {
    nmatrix<double> TH, TL;

    GlobalToGlobalData();
    EQ.SetParameterData(__QP);

    const IntegratorQ1Q2<3> *I =
        dynamic_cast<const IntegratorQ1Q2<3> *>(GetIntegrator());
    assert(I);

    const FemInterface &HighOrderFem(*GetFem());

    for (int iq = 0; iq < GetPatchMesh()->npatches(); ++iq) {
      Transformation(TH, iq);
      TransformationQ1(TL, iq);

      HighOrderFem.ReInit(TH);
      LowOrderFem.ReInit(TL);

      PatchDiscretization::GlobalToLocal(__U, u, iq);
      if (__ADJOINT) {
        if (FluidInterfaceCell(iq))
          for (int i = 0; i < __DEL_F.size(); ++i)
            DeleteTestFunctionsVector(__U, FluidInterfaceNodes(iq), __DEL_F[i]);
        if (SolidInterfaceCell(iq))
          for (int i = 0; i < __DEL_S.size(); ++i)
            DeleteTestFunctionsVector(__U, SolidInterfaceNodes(iq), __DEL_S[i]);
      }

      I->Form(EQ, __F, HighOrderFem, LowOrderFem, __U, __QN, __QC);

      if (!__ADJOINT) {
        if (FluidInterfaceCell(iq))
          for (int i = 0; i < __DEL_F.size(); ++i)
            DeleteTestFunctionsVector(__F, FluidInterfaceNodes(iq), __DEL_F[i]);
        if (SolidInterfaceCell(iq))
          for (int i = 0; i < __DEL_S.size(); ++i)
            DeleteTestFunctionsVector(__F, SolidInterfaceNodes(iq), __DEL_S[i]);
      }
      PatchDiscretization::LocalToGlobal(f, __F, iq, d);
    }
  }
};

} // namespace Gascoigne

/*----------------------------   alediscretization.h
 * ---------------------------*/
/* end of #ifndef __alediscretization_H */
#endif
/*----------------------------   alediscretization.h
 * ---------------------------*/
