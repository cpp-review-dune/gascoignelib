#include "regular_update.h"
#include <algorithm>


using namespace std;
using namespace Gascoigne;

/*----------------------------------------------*/

void regular_update(IntSet& hr, IntSet& hc, IntVector& vr, IntVector& vc)
{
  for (int i=0; i<vr.size(); i++)
    {
      hr.insert(vr[i]);
    }
  for (int i=0; i<vc.size(); i++)
    {
      hc.erase(vc[i]);
    }
  vr.insert(vr.end(),vc.begin(),vc.end());
  sort(vr.begin(),vr.end());
  int n = unique(vr.begin(),vr.end()) - vr.begin();
  vr.reserve(n); vr.resize(n);
}
