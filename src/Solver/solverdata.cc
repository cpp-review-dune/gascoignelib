#include  "solverdata.h"
#include  "filescanner.h"


using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

SolverData::SolverData() :
  pfilter(0)
{
}

/*-----------------------------------------*/

SolverData::~SolverData()
{
}

/*-----------------------------------------*/

void SolverData::Init(const ParamFile* pf, int ncomp)
{
  ilum.reservesize(ncomp);
  ilum.zero();

  DataFormatHandler DFH;

  DFH.insert("pfilter" , &pfilter);
  DFH.insert("exact_lu", &exact_lu, 0);
  DFH.insert("enlarge" , &enlarge , 0);
  DFH.insert("iterpre" , &iter_pre , 4);
  DFH.insert("iterpost", &iter_post, 4);
  DFH.insert("iterexact",&iter_exact , 10);
  DFH.insert("omega"   , &omega, 1.);
  DFH.insert("ilum"    , &ilum);
  
//   DFH.insert("ilusort", &ilusort,"vectordirection");
  DFH.insert("ilusort", &ilusort,"cuthillmckee");
  DFH.insert("stream_direction",&stream_direction);
  DFH.insert("vector_direction",&vector_direction);
  
  DFH.insert("linear_smooth",     &linear_smooth,      "ilu");
  DFH.insert("bicgstab_residual" ,&bicgstab_residual , "approx");
  DFH.insert("bicgstab_pstep" ,   &bicgstab_pstep ,     0);
  DFH.insert("bicgstab_miniter",  &bicgstab_miniter,    1.);
  
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(pf,"Solver");

  if ((ilusort=="streamdirection")&&(stream_direction.size()==0))
    {
      cerr << "Bei \n\tilusort\tstreamdiretion\nmuss" << endl
	   << "\tstream_direction\n"
	   << "mit Komponenten, nach denen sortiert wird, "
	   << "angegeben werden.\n";
      abort();
    }

  assert(ilum.size()==ncomp);
  assert(stream_direction.size()<=ncomp);
}
