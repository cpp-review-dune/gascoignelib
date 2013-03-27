#include "meshhierarchy.h"
#include <fstream> 

using namespace std;


// ==================================================

namespace Tsuchimikado
{
  // ==================================================
  
  template<int DIM>
  void MeshHierarchy<DIM>::ReInit()
  {
    this->clear();
    this->push_back(MeshLevel<DIM>(__TC));
    GetMeshLevel(0).init_active();

    do
      {
	MeshLevel<DIM> ML(__TC);
	ML.init_from_meshlevel(GetMeshLevel(this->size()-1));
	push_back(ML);
      } 
    while(GetMeshLevel(this->size()-1).size()<GetMeshLevel(this->size()-2).size());
    this->pop_back();    
  }

  // ==================================================

  template<int DIM>
  void MeshHierarchy<DIM>::print_gnuplot(const std::string& fname) const
  {
    for (int i=0;i<nlevels();++i)
      {
	char s[20];
	sprintf(s,"%s_%i",fname.c_str(),i);
	GetMeshLevel(i).print_gnuplot2(s);
      }
  }
 
  template class MeshHierarchy<2>;
  template class MeshHierarchy<3>;

}
