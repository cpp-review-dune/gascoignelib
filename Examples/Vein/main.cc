
#include "boundaryfunctional.h"
#include "dirichletdatabycolor.h"
#include "domainfunctional.h"
#include "local.h"
#include "loop.h"
#include "residualfunctional.h"
#include "weightedpointfunctional.h"

#include "../eigen3/Eigen/Dense"
#include "boundaryfunctional.h"
#include "filescanner.h"
#include "paramfile.h"
using namespace Gascoigne;
using namespace std;

class Drag : public virtual ResidualFunctional {
  std::string GetName() const { return "drag"; }

public:
  Drag() {
    __comps.push_back(1);
    __scales.push_back(1.0);
    __cols.insert(80);
    __cols.insert(81);
    __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
  }
};
class Lift : public virtual ResidualFunctional {
  std::string GetName() const { return "lift"; }

public:
  Lift() {
    __comps.push_back(2);
    __scales.push_back(1.0);
    __cols.insert(80);
    __cols.insert(81);
    __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
  }
};

class DragInterface : public virtual ResidualFunctional {
  std::string GetName() const { return "draginterface"; }

public:
  DragInterface() {
    __comps.push_back(1);
    __scales.push_back(1.0);
    __cols.insert(80);
    __cols.insert(90);

    __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
  }
};
class LiftInterface : public virtual ResidualFunctional {
  std::string GetName() const { return "liftinterface"; }

public:
  LiftInterface() {
    __comps.push_back(2);
    __scales.push_back(1.0);
    __cols.insert(80);
    __cols.insert(90);

    __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
  }
};

extern double __DT;
class DomainFunctionalAddResidualDrag : public Gascoigne::DomainFunctional {
private:
protected:
  mutable Gascoigne::FemFunction *OLD;
  mutable double rho_s, domain;

public:
  DomainFunctionalAddResidualDrag(const Gascoigne::ParamFile *paramfile) {
    DataFormatHandler DFH;
    DFH.insert("rho_s", &rho_s, 0.0);
    FileScanner FS(DFH, paramfile, "Equation");
  }
  ~DomainFunctionalAddResidualDrag() {}

  void SetFemData(Gascoigne::FemData &q) const {
    assert(q.find("OLD") != q.end());
    OLD = &q["OLD"];
  }

  double J(const Gascoigne::FemFunction &U,
           const Gascoigne::Vertex2d &v) const {
    if (domain > 0) // solid
    {
      return rho_s * (U[1].m() - (*OLD)[1].m()) / __DT;
    } else
      return 0;
  }
  void point_cell(int material) const {
    if (material == 1)
      domain = 1;
    if (material == 2)
      domain = -1;
  }
  std::string GetName() const { return "DomainFunctionalVorticity"; }
};
class DomainFunctionalAddResidualLift : public Gascoigne::DomainFunctional {
private:
protected:
  mutable Gascoigne::FemFunction *OLD;
  mutable double rho_s, domain;

public:
  DomainFunctionalAddResidualLift(const Gascoigne::ParamFile *paramfile) {
    DataFormatHandler DFH;
    DFH.insert("rho_s", &rho_s, 0.0);
    FileScanner FS(DFH, paramfile, "Equation");
  }
  ~DomainFunctionalAddResidualLift() {}

  void SetFemData(Gascoigne::FemData &q) const {
    assert(q.find("OLD") != q.end());
    OLD = &q["OLD"];
  }

  double J(const Gascoigne::FemFunction &U,
           const Gascoigne::Vertex2d &v) const {
    if (domain > 0) // solid
    {
      return rho_s * (U[2].m() - (*OLD)[2].m()) / __DT;
    } else
      return 0;
  }
  void point_cell(int material) const {
    if (material == 1)
      domain = 1;
    if (material == 2)
      domain = -1;
  }

  std::string GetName() const { return "DomainFunctionalVorticity"; }
};

/**********************************************************/

class BoundaryFunctionalDrag : public Gascoigne::BoundaryFunctional {
private:
protected:
  typedef Eigen::Matrix<double, 2, 2> Matrix2d;
  typedef Eigen::Matrix<double, 2, 1> Vector2d;

  mutable Vector2d V, __U, _n, F_fluid, F_solid;
  mutable Matrix2d NV, NU, SIGMA_F, F, E, SIGMA_S;
  mutable double __J, P;
  double __mu_s, __nu_f, __rho_f, __rho_s, __lambda_s;
  mutable FemFunction *DEF;

public:
  BoundaryFunctionalDrag(const Gascoigne::ParamFile *paramfile) {
    DataFormatHandler DFH;
    DFH.insert("nu_f", &__nu_f);
    DFH.insert("rho_f", &__rho_f);
    FileScanner FS(DFH, paramfile, "Equation");
  }
  ~BoundaryFunctionalDrag() {}

  void SetFemData(FemData &q) const {

    assert(q.find("DEF") != q.end());
    DEF = &q["DEF"];
  }
  double J(const Gascoigne::FemFunction &U, const Gascoigne::Vertex2d &v,
           const Gascoigne::Vertex2d &n, int color) const {
    P = U[0].m();
    _n << n.x(), n.y();
    // Initialisieren der Vektoren und Matrizen
    NU << (*DEF)[1].x(), (*DEF)[1].y(), (*DEF)[2].x(), (*DEF)[2].y();
    NV << U[1].x(), U[1].y(), U[2].x(), U[2].y();

    // Deformationsgradient
    F = Matrix2d::Identity() + NU;
    __J = F.determinant();

    // Spannungstensors des Fluids J*sigma_f F^{-T}
    SIGMA_F =
        -__J * P * F.inverse().transpose() +
        __J * __rho_f * __nu_f *
            (NV * F.inverse() + F.inverse().transpose() * NV.transpose()) *
            F.inverse().transpose();

    F_fluid = -SIGMA_F * _n;

    if (color == 80 || color == 90) {
      // F_D Komponente 0
      // F_L Komponenete 1
      return F_fluid(0);
    } else {
      return 0.;
    }
  }

  std::string GetName() const { return "BoundaryFunctionalDrag"; }
};

/**********************************************************/

class BoundaryFunctionalLift : public Gascoigne::BoundaryFunctional {
private:
protected:
  typedef Eigen::Matrix<double, 2, 2> Matrix2d;
  typedef Eigen::Matrix<double, 2, 1> Vector2d;

  mutable Vector2d V, __U, _n, F_fluid, F_solid;
  mutable Matrix2d NV, NU, SIGMA_F, F, E, SIGMA_S;
  mutable double __J, P;
  double __mu_s, __nu_f, __rho_f, __rho_s, __lambda_s;
  mutable FemFunction *DEF;

public:
  BoundaryFunctionalLift(const Gascoigne::ParamFile *paramfile) {
    DataFormatHandler DFH;
    DFH.insert("nu_f", &__nu_f);
    DFH.insert("rho_f", &__rho_f);
    FileScanner FS(DFH, paramfile, "Equation");
  }
  ~BoundaryFunctionalLift() {}

  void SetFemData(FemData &q) const {

    assert(q.find("DEF") != q.end());
    DEF = &q["DEF"];
  }
  double J(const Gascoigne::FemFunction &U, const Gascoigne::Vertex2d &v,
           const Gascoigne::Vertex2d &n, int color) const {

    P = U[0].m();
    // Initialisieren der Vektoren und Matrizen
    NU << (*DEF)[1].x(), (*DEF)[1].y(), (*DEF)[2].x(), (*DEF)[2].y();
    NV << U[1].x(), U[1].y(), U[2].x(), U[2].y();

    // Deformationsgradient
    F = Matrix2d::Identity() + NU;
    __J = F.determinant();

    // Spannungstensors des Fluids J*sigma_f F^{-T}
    SIGMA_F =
        -__J * P * F.inverse().transpose() +
        __J * __rho_f * __nu_f *
            (NV * F.inverse() + F.inverse().transpose() * NV.transpose()) *
            F.inverse().transpose();

    F_fluid = -SIGMA_F * _n;

    if (color == 80 || color == 90) {
      // F_D Komponente 0
      // F_L Komponenete 1
      return F_fluid(1);
    } else {
      return 0.;
    }
  }

  std::string GetName() const { return "BoundaryFunctionalLift"; }
};

#include "FSI/FSI_reduced/boundaryfsi.h"
#include "FSI/FSI_reduced/fsi.h"
#include "FSI/FSI_reduced/fsi_dirichlet.h"

/*---------------------------------------------------*/
class ProblemDescriptorFSI2d : public ProblemDescriptorBase {
public:
  std::string GetName() const { return "fsi"; }
  void BasicInit(const ParamFile *pf) {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new FSI<2>(GetParamFile());
    //    GetBoundaryEquationPointer() = new FSI<2>(GetParamFile());
    GetDirichletDataPointer() = new MyDDFSI(GetParamFile());

    ProblemDescriptorBase::BasicInit(pf);

    //    GetComponentInformationPointer() = new FSI_CI<2>;
  }
};
/*---------------------------------------------------*/

int main(int argc, char **argv) {
  ParamFile pf("fsi-3.param");
  if (argc == 2)
    pf.SetName(argv[1]);

  ProblemDescriptorFSI2d Problem2d;
  Problem2d.BasicInit(&pf);

  ProblemContainer PC2d;
  PC2d.AddProblem("fsi", &Problem2d);

  FunctionalContainer FC2d;
  WeightedPointFunctional Pux, Puy;
  vector<Vertex2d> v2d;
  v2d.push_back(Vertex2d(0.6, 0.2));
  vector<int> cx, cy;
  cx.push_back(0);
  cy.push_back(1);
  vector<double> weigh;
  weigh.push_back(1.0);
  Pux.BasicInit(v2d, cx, weigh);
  Puy.BasicInit(v2d, cy, weigh);

  WeightedPointFunctional Pvx, Pvy;
  vector<int> cx2, cy2;
  cx2.push_back(1);
  cy2.push_back(2);
  Pvx.BasicInit(v2d, cx2, weigh);
  Pvy.BasicInit(v2d, cy2, weigh);

  Drag drag;
  Lift lift;
  DragInterface draginterface;
  LiftInterface liftinterface;
  DomainFunctionalAddResidualDrag dfAddResidualdrag(&pf);
  DomainFunctionalAddResidualLift dfAddResiduallift(&pf);

  BoundaryFunctionalDrag boundaryfunctionaldrag(&pf);
  BoundaryFunctionalLift boundaryfunctionallift(&pf);

  FC2d.AddFunctional("ux", &Pux);
  FC2d.AddFunctional("uy", &Puy);
  FC2d.AddFunctional("vx", &Pvx);
  FC2d.AddFunctional("vy", &Pvy);

  FC2d.AddFunctional("drag", &drag);
  FC2d.AddFunctional("lift", &lift);
  FC2d.AddFunctional("draginterface", &draginterface);
  FC2d.AddFunctional("liftinterface", &liftinterface);
  FC2d.AddFunctional("dfAddResidualdrag", &dfAddResidualdrag);
  FC2d.AddFunctional("dfAddResiduallift", &dfAddResiduallift);
  FC2d.AddFunctional("BoundaryFunctionalDrag", &boundaryfunctionaldrag);
  FC2d.AddFunctional("BoundaryFunctionalLift", &boundaryfunctionallift);

  Loop<2> loop;

  loop.BasicInit(&pf, &PC2d, &FC2d);
  loop.run("fsi");

  return 0;
}
