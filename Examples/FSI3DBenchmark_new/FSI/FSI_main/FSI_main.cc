#include  "FSI_main.h"
#include  "filescanner.h"

extern double __DT;
extern double __THETA;
extern double __TIME;
extern bool InterfaceResidual;


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  template class FSI_main<2>;
  template class FSI_main<3>;
  
}

/*-----------------------------------------*/
