#include "regular_update.h"
#include <algorithm>

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
  std::sort(vr.begin(),vr.end());
  int n = std::unique(vr.begin(),vr.end()) - vr.begin();
  vr.reserve(n); vr.resize(n);
}
