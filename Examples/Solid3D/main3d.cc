
#include "domainfunctional.h"
#include "local.h"
#include "loop.h"
using namespace Gascoigne;
using namespace std;

/*---------------------------------------------------*/

class DomainFunctionalDisp : public Gascoigne::DomainFunctional
{
private:
protected:
public:
  DomainFunctionalDisp(const Gascoigne::ParamFile* paramfile) {}
  ~DomainFunctionalDisp() {}

  double J(const Gascoigne::FemFunction& U, const Gascoigne::Vertex3d& v) const
  {
    return (U[0].m() * U[0].m()) + (U[1].m() * U[1].m()) +
           (U[2].m() * U[2].m());
  }

  std::string GetName() const { return "DomainFunctionalDisp"; }
};

int
main(int argc, char** argv)
{
  ParamFile pf("box3d.param");
  if (argc == 2)
    pf.SetName(argv[1]);

  ProblemDescriptor3d Problem3d;
  Problem3d.BasicInit(&pf);

  ProblemContainer PC3d;
  PC3d.AddProblem("solid", &Problem3d);
  FunctionalContainer FC3d;
  DomainFunctionalDisp dfDisp(&pf);
  FC3d.AddFunctional("DomainFunctionalDisp", &dfDisp);

  Loop<3> loop;

  loop.BasicInit(&pf, &PC3d, &FC3d);

  loop.run("solid");

  return 0;
}
