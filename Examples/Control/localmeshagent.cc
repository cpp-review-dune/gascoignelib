#include  "localmeshagent.h"
#include  "filescanner.h"
#include  "gascoignemeshconstructor.h"

using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

void LocalMeshAgent::BasicInit(const ParamFile* paramfile)
{
  cerr << "\t§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§\n";
  cerr << "\t i don't read any paramfile!\n";
  cerr << "\t§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§\n";

  DataFormatHandler DFH;
  string gridname("Results/forward.00000.gup");
  DFH.insert("gridname" ,&gridname,"none");
  DFH.insert("dimension",&dimension,0);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(paramfile,"Mesh");

  cerr << "\t§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§\n";
  cerr << "\t i read Mesh from ";
  cerr << gridname << endl;
  cerr << "\t§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§\n";

  if (dimension==2)
    {
      HMP = new HierarchicalMesh2d;
    }
  else if (dimension==3)
    {
      HMP = new HierarchicalMesh3d;
    }
  else
    {
      cout << "dimension of Mesh ? " << dimension << endl;
    }

  int patchdepth=1;
  int epatcher=1;
  HMP->SetParameters(gridname,patchdepth,epatcher);
  assert(HMP);

  GMG = NewMultiGridMesh();

  ReInit();
}
