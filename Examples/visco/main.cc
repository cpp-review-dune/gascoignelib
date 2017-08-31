
#include  "local.h"
#include  "loop.h"
#include "weightedpointfunctional.h"
#include "boundaryfunctional.h"
#include  "residualfunctional.h"
#include  "dirichletdatabycolor.h"

using namespace Gascoigne;
using namespace std;

double BOUNDARY;


class Drag : public virtual ResidualFunctional
{ 
  std::string GetName() const { return "drag"; } 
public:
  Drag()
  {
    __comps.push_back(1);
    __scales.push_back(1.0);
    __cols.insert(80);   
    __cols.insert(81);	  
    __DD  = new DirichletDataByColor(GetComps(),GetColors(),GetScales());

  }
};
class Lift : public virtual ResidualFunctional
{ 
  std::string GetName() const { return "lift"; } 
public:
  Lift()
  {
    __comps.push_back(2);
    __scales.push_back(1.0);
    __cols.insert(80);   
    __cols.insert(81);	 
    __DD  = new DirichletDataByColor(GetComps(),GetColors(),GetScales());

  }
};

class KinEnergy : public virtual DomainFunctional
{
  std::string GetName() const { return "kinetic"; } 
public:
  double J(const FemFunction& U, const Vertex2d& v) const 
  {
    return 0.5*(U[1].m()*U[1].m()+U[2].m()*U[2].m());
  }

};
class ElEnergy : public virtual DomainFunctional
{
  std::string GetName() const { return "elastic"; } 
public:
  double J(const FemFunction& U, const Vertex2d& v) const 
  {
    return U[0].m() + U[2].m();
  }

};



/*---------------------------------------------------*/



int main(int argc, char** argv)
{
  ParamFile pf("visco.param");
  if (argc==2)
    pf.SetName(argv[1]);

  
  VelProblem vel;
  vel.BasicInit(&pf);
  StressProblem stress;
  stress.BasicInit(&pf);

  ProblemContainer PC2d;
  PC2d.AddProblem("vel", &vel);
  PC2d.AddProblem("stress", &stress);

  FunctionalContainer FC2d;
  
  // Drag drag;
  // Lift lift;
  KinEnergy kin;
  ElEnergy el;
  
  FC2d.AddFunctional("0 kin",&kin);
  FC2d.AddFunctional("1 el",&el);

  Loop<2> loop;
  loop.BasicInit(&pf,&PC2d,&FC2d);
  loop.run("vel");

  return 0;
}
