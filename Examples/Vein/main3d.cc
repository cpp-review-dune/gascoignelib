
#include "boundaryfunctional.h"
#include "dirichletdatabycolor.h"
#include "local.h"
#include "loop.h"
#include "residualfunctional.h"
#include "weightedpointfunctional.h"

using namespace Gascoigne;
using namespace std;

class Drag : public virtual ResidualFunctional {
  std::string GetName() const { return "drag"; }

public:
  Drag() {
    __comps.push_back(1);
    __scales.push_back(1.0);
    __cols.insert(84);
    __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
  }
};
class Lift : public virtual ResidualFunctional {
  std::string GetName() const { return "drag"; }

public:
  Lift() {
    __comps.push_back(2);
    __scales.push_back(1.0);
    __cols.insert(84);
    __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
  }
};

// -----------------------------------------

#include "FSI/FSI_reduced/boundaryfsi.h"
#include "FSI/FSI_reduced/fsi.h"

class ProblemDescriptorFSI3d : public ProblemDescriptorBase {
public:
  std::string GetName() const { return "fsi_reduced"; }
  void BasicInit(const ParamFile *pf) {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new FSI<3>(GetParamFile());
    GetBoundaryEquationPointer() = new BoundaryFSI<3>(GetParamFile());
    GetDirichletDataPointer() = NULL;

    ProblemDescriptorBase::BasicInit(pf);
    GetComponentInformationPointer() = NULL;
  }
};
/*---------------------------------------------------*/
#include "FSI/MeshMotion/meshmotion.h"

class ProblemDescriptormMeshMotion : public ProblemDescriptorBase {
public:
  std::string GetName() const { return "meshmotion"; }
  void BasicInit(const ParamFile *pf) {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new MeshMotion<3>(GetParamFile());
    GetBoundaryEquationPointer() = NULL;
    GetDirichletDataPointer() = NULL;

    ProblemDescriptorBase::BasicInit(pf);
  }
};
/*---------------------------------------------------*/
#include "FSI/Deformation_Solid/def_solid.h"

class ProblemDescriptormDef_Solid : public ProblemDescriptorBase {
public:
  std::string GetName() const { return "def_solid"; }
  void BasicInit(const ParamFile *pf) {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new Def_Solid<3>(GetParamFile());
    GetBoundaryEquationPointer() = NULL;
    GetDirichletDataPointer() = NULL;

    ProblemDescriptorBase::BasicInit(pf);
  }
};
/*---------------------------------------------------*/
#include "FSI/FSI_main/FSI_main.h"
#include "FSI/FSI_main/FSI_main_dirichlet.h"

class ProblemDescriptormFSI_main : public ProblemDescriptorBase {
public:
  std::string GetName() const { return "fsi_main"; }
  void BasicInit(const ParamFile *pf) {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new FSI_main<3>(GetParamFile());
    GetBoundaryEquationPointer() = NULL;
    GetDirichletDataPointer() = new FSI_main_MyDD3d(GetParamFile());

    ProblemDescriptorBase::BasicInit(pf);
    GetComponentInformationPointer() = new FSI_CI<3>;
  }
};
/*---------------------------------------------------*/

#include "Fluid_Stat/boundary_fluid_stat.h"
#include "Fluid_Stat/fluid_stat.h"
#include "Fluid_Stat/fluid_stat_dirichlet.h"

class ProblemDescriptorFluid_Stat3d : public ProblemDescriptorBase {
public:
  std::string GetName() const { return "fluid_stat"; }
  void BasicInit(const ParamFile *pf) {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new Fluid_Stat<3>(GetParamFile());
    GetBoundaryEquationPointer() = new Boundary_Fluid_Stat<3>(GetParamFile());
    GetDirichletDataPointer() = new DD_Fluid_Stat3d(GetParamFile());

    ProblemDescriptorBase::BasicInit(pf);

    GetComponentInformationPointer() = new Fluid_CI<3>;
  }
};

/*---------------------------------------------------*/

#include "Solid_Euler/boundarysolideuler.h"
#include "Solid_Euler/solid_euler_dirichlet.h"
#include "Solid_Euler/solideuler.h"

class ProblemDescriptorSolid_Euler3d : public ProblemDescriptorBase {
public:
  std::string GetName() const { return "solid_euler"; }
  void BasicInit(const ParamFile *pf) {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new SolidEuler<3>(GetParamFile());
    GetBoundaryEquationPointer() = new BoundarySolidEuler<3>(GetParamFile());
    GetDirichletDataPointer() = new DD_Solid_Euler3d(GetParamFile());

    ProblemDescriptorBase::BasicInit(pf);

    GetComponentInformationPointer() = new Solid_Euler_CI<3>;
  }
};

/*---------------------------------------------------*/

int main(int argc, char **argv) {
  ParamFile pf_fsi_red("vein_fsi_reduced.param");
  ParamFile pf_meshmotion("vein_meshmotion.param");
  ParamFile pf_def_solid("vein_def_solid.param");

  ParamFile pf_fsi_main("vein_fsi_reduced.param");
  ParamFile pf("vein_fsi_reduced.param");

  // if (argc==2)
  //  pf.SetName(argv[1]);

  ProblemDescriptorFSI3d ProblemFSI3d;
  ProblemFSI3d.BasicInit(&pf_fsi_red);

  ProblemDescriptormFSI_main ProblemFSI_main;
  ProblemFSI_main.BasicInit(&pf_fsi_main);

  ProblemDescriptormMeshMotion ProblemMeshMotion;
  ProblemMeshMotion.BasicInit(&pf_meshmotion);

  ProblemDescriptormDef_Solid ProblemDef_Solid;
  ProblemDef_Solid.BasicInit(&pf_def_solid);

  ProblemDescriptorFluid_Stat3d ProblemFluid_Stat3d;
  ProblemFluid_Stat3d.BasicInit(&pf);

  ProblemDescriptorSolid_Euler3d ProblemSolid_Euler3d;
  ProblemSolid_Euler3d.BasicInit(&pf);

  ProblemContainer PC3d;
  PC3d.AddProblem("fsi_reduced", &ProblemFSI3d);
  PC3d.AddProblem("fsi_main", &ProblemFSI_main);
  PC3d.AddProblem("meshmotion", &ProblemMeshMotion);
  PC3d.AddProblem("def_solid", &ProblemDef_Solid);

  PC3d.AddProblem("fluid_stat", &ProblemFluid_Stat3d);
  PC3d.AddProblem("solid_euler", &ProblemSolid_Euler3d);

  FunctionalContainer FC3d;
  WeightedPointFunctional Ux;
  WeightedPointFunctional Uy;
  WeightedPointFunctional Uz;
  vector<Vertex3d> v1;
  v1.push_back(Vertex3d(0.45, 0.15, 0.15));

  vector<int> cx;
  cx.push_back(4);
  vector<int> cy;
  cy.push_back(5);
  vector<int> cz;
  cz.push_back(6);

  vector<double> weigh;
  weigh.push_back(1.0);
  /*Ux.BasicInit(v1,cx,weigh);
  Uy.BasicInit(v1,cy,weigh);
  Uz.BasicInit(v1,cz,weigh);

  Drag drag;
  Lift lift;

  FC3d.AddFunctional("ux",&Ux);
  FC3d.AddFunctional("uy",&Uy);
  FC3d.AddFunctional("uz",&Uz);
  FC3d.AddFunctional("drag",&drag);
  FC3d.AddFunctional("lift",&lift);
  */

  Loop<3> loop;

  loop.BasicInit(&pf, &PC3d, &FC3d);
  loop.run("fsi");

  return 0;
}
