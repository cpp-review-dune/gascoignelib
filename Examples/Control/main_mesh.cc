#include  "meshagent.h"

using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh.param");

  MeshAgent MA;
  MA.BasicInit(&paramfile);
  system("rm -rf Results");
  system("mkdir Results");
  MA.write_gup("Results/mesh");

  return 0;
}
