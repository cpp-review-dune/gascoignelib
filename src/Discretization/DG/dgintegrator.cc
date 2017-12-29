#include "dgintegrator.h"


namespace Gascoigne
{
  template <int DIM>
  DGIntegrator<DIM>::DGIntegrator()
      : GalerkinIntegrator<DIM>()
  {
  }


  template <>
  void DGIntegrator<2>::EdgeForm(bool internaledge,
                                 const DGEquation &EQ,
                                 LocalVector &F1,
                                 LocalVector &F2,
                                 const FemInterface &FEMASTER,
                                 const FemInterface &FESLAVE,
                                 const int &masterli,
                                 const int &slaveli,
                                 const LocalVector &U1,
                                 const LocalVector &U2,
                                 const LocalData &Q,
                                 const LocalData &QC) const
  {
    F1.ReInit(EQ.GetNcomp(), FEMASTER.n());
    F1.zero();
    if (internaledge)
    {
      F2.ReInit(EQ.GetNcomp(), FESLAVE.n());
      F2.zero();
    }

    const LineGauss2 IF1, IF2;

    Vertex<2> x1, x2, n1, n2;
    Vertex<1> xi1, xi2;

    FemFunction UH1, UH2;

    for (int k = 0; k < IF1.n(); k++)
    {
      // initialize quadrature rule and element
      IF1.xi(xi1, k);
      FEMASTER.point_boundary(masterli, xi1);
      FEMASTER.x(x1);
      FEMASTER.normal(n1);
      universal_point(FEMASTER, UH1, U1);
      double h = FEMASTER.G();
      double weight = IF1.w(k) * h;

      // slave, check ordering
      if (internaledge)
      {
        IF2.xi(xi2, k);
        FESLAVE.point_boundary(slaveli, xi2);
        FESLAVE.x(x2);
        assert(fabs(x1[0] - x2[0]) < 1.e-12);
        assert(fabs(x1[1] - x2[1]) < 1.e-12);
        FESLAVE.normal(n2);
        assert(fabs(n1[0] + n2[0]) < 1.e-12);
        assert(fabs(n1[1] + n2[1]) < 1.e-12);

        universal_point(FESLAVE, UH2, U2);
        assert(fabs(h - FESLAVE.G()) < 1.e-12);
      }
      EQ.point_edge(internaledge, h, UH1, UH2, x1, n1);

      TestFunction N;
      for (int i1 = 0; i1 < FEMASTER.n(); ++i1)
      {
        N.zero();
        FEMASTER.init_test_functions(N, weight, i1);
        EQ.EdgeForm1(F1.start(i1), UH1, UH2, N);
      }
      if (internaledge)
        for (int i2 = 0; i2 < FESLAVE.n(); ++i2)
        {
          N.zero();
          FESLAVE.init_test_functions(N, weight, i2);
          EQ.EdgeForm2(F2.start(i2), UH1, UH2, N);
        }
    }
  }


  // Die Matrix-Funktion wird 4mal aufgerufen,
  // jeweils Kopplungen zwischen den einzelnen Test-Ansatzfunktionen
  template <>
  void DGIntegrator<2>::EdgeMatrix(bool internaledge,
                                   const DGEquation &EQ,
                                   EntryMatrix &E11,
                                   EntryMatrix &E12,
                                   EntryMatrix &E21,
                                   EntryMatrix &E22,
                                   const FemInterface &FEMASTER,
                                   const FemInterface &FESLAVE,
                                   const int &masterli,
                                   const int &slaveli,
                                   const LocalVector &U1,
                                   const LocalVector &U2,
                                   const LocalData &Q,
                                   const LocalData &QC) const
  {
    E11.SetDimensionDof(FEMASTER.n(), FEMASTER.n());
    E11.SetDimensionComp(U1.ncomp(), U1.ncomp());
    E11.resize();
    E11.zero();

    if (internaledge)
    {
      E12.SetDimensionDof(FEMASTER.n(), FESLAVE.n());
      E12.SetDimensionComp(U1.ncomp(), U2.ncomp());
      E12.resize();
      E12.zero();

      E21.SetDimensionDof(FESLAVE.n(), FEMASTER.n());
      E21.SetDimensionComp(U2.ncomp(), U1.ncomp());
      E21.resize();
      E21.zero();

      E22.SetDimensionDof(FESLAVE.n(), FESLAVE.n());
      E22.SetDimensionComp(U2.ncomp(), U2.ncomp());
      E22.resize();
      E22.zero();
    }

    const LineGauss2 IF1, IF2;

    Vertex<2> x1, x2, n1, n2;
    Vertex<1> xi1, xi2;

    FemFunction UH1, UH2;

    for (int k = 0; k < IF1.n(); k++)
    {
      // initialize quadrature rule and element
      IF1.xi(xi1, k);
      FEMASTER.point_boundary(masterli, xi1);
      FEMASTER.x(x1);
      FEMASTER.normal(n1);
      universal_point(FEMASTER, UH1, U1);
      double h = FEMASTER.G();
      double weight = IF1.w(k) * h;

      // slave, check ordering
      if (internaledge)
      {
        IF2.xi(xi2, k);
        FESLAVE.point_boundary(slaveli, xi2);
        FESLAVE.x(x2);
        assert(fabs(x1[0] - x2[0]) < 1.e-12);
        assert(fabs(x1[1] - x2[1]) < 1.e-12);
        FESLAVE.normal(n2);
        assert(fabs(n1[0] + n2[0]) < 1.e-12);
        assert(fabs(n1[1] + n2[1]) < 1.e-12);

        universal_point(FESLAVE, UH2, U2);
        assert(fabs(h - FESLAVE.G()) < 1.e-12);
      }
      EQ.point_edge(internaledge, h, UH1, UH2, x1, n1);


      // init test and trial
      std::vector<TestFunction> NNN1(FEMASTER.n());
      std::vector<TestFunction> NNN2(FESLAVE.n());
      for (int i1 = 0; i1 < FEMASTER.n(); ++i1)
      {
        NNN1[i1].zero();
        FEMASTER.init_test_functions(NNN1[i1], sqrt(weight), i1);
      }
      
      for (int j = 0; j < FEMASTER.n(); j++)
	for (int i = 0; i < FEMASTER.n(); i++)
          {
            E11.SetDofIndex(i, j);
            EQ.EdgeMatrix11(E11, UH1, UH2, NNN1[j], NNN1[i]);
          }

      if (internaledge)
      {
        for (int i2 = 0; i2 < FESLAVE.n(); ++i2)
        {
          NNN2[i2].zero();
          FESLAVE.init_test_functions(NNN2[i2], sqrt(weight), i2);
        }
        for (int j = 0; j < FESLAVE.n(); j++)
          for (int i = 0; i < FEMASTER.n(); i++)
          {
            E12.SetDofIndex(i, j);
            EQ.EdgeMatrix12(E12, UH1, UH2, NNN2[j], NNN1[i]);
          }
        for (int j = 0; j < FEMASTER.n(); j++)
          for (int i = 0; i < FESLAVE.n(); i++)
          {
            E21.SetDofIndex(i, j);
            EQ.EdgeMatrix21(E21, UH1, UH2, NNN1[j], NNN2[i]);
          }
        for (int j = 0; j < FESLAVE.n(); j++)
          for (int i = 0; i < FESLAVE.n(); i++)
          {
            E22.SetDofIndex(i, j);
            EQ.EdgeMatrix22(E22, UH1, UH2, NNN2[j], NNN2[i]);
          }
      }
    }

  }
  
  template class DGIntegrator<2>;
  template class DGIntegrator<3>;
}
