#include  "solverdata.h"
#include  "filescanner.h"


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  
  /*-----------------------------------------*/
  
  void SolverData::BasicInit(const ParamFile* pf, const std::string label)
  {
    DataFormatHandler DFH;
    
    _pfilter.resize(0);

    DFH.insert("pfilter" , &_pfilter);
    DFH.insert("exact_lu", &exact_lu, 0);
    DFH.insert("enlarge" , &enlarge , 0);
    DFH.insert("iterpre" , &iter_pre , 4);
    DFH.insert("iterpost", &iter_post, 4);
    DFH.insert("iterexact",&iter_exact , 10);
    DFH.insert("omega"   , &omega, 1.);
    DFH.insert("ilum"    , &_ilum);
    
    //   DFH.insert("ilusort", &ilusort,"vectordirection");
    DFH.insert("ilusort", &ilusort,"cuthillmckee");
    DFH.insert("stream_direction",&stream_direction);
    DFH.insert("vector_direction",&vector_direction);
    
    DFH.insert("linear_smooth",     &linear_smooth,      "ilu");
    DFH.insert("bicgstab_residual" ,&bicgstab_residual , "approx");
    DFH.insert("bicgstab_pstep" ,   &bicgstab_pstep ,     0);
    DFH.insert("bicgstab_miniter",  &bicgstab_miniter,    1.);

    DFH.insert("cgmass_maxiter",    &_cgMassMaxIter,      100);
    DFH.insert("cgmass_tol",        &_cgMassTol,          1.e-8);
    DFH.insert("cgmass_globaltol",  &_cgMassGlobalTol,    1.e-14);
    
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf,label);
    
    if ((ilusort=="streamdirection")&&(stream_direction.size()==0))
      {
        cerr << "Bei \n\tilusort\tstreamdiretion\nmuss" << endl
             << "\tstream_direction\n"
             << "mit Komponenten, nach denen sortiert wird, "
             << "angegeben werden.\n";
        abort();
      }
  }

  
}
