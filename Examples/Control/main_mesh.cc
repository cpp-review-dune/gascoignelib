#include  "meshagent.h"

using namespace std;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("mesh.param");

  MeshAgent MA;
  MA.BasicInit(&paramfile);
  MA.write_gup("Results/mesh");

  return 0;
}
