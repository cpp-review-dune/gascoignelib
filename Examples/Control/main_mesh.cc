#include  "meshagent.h"

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  std::string paramfile("mesh.param");
  if(argc>=2) {
    paramfile = argv[1];
  }

  MeshAgent MA;
  MA.BasicInit(paramfile);
  MA.write_gup("Results/mesh");

  return 0;
}
