#include  "set2vec.h"

using namespace std;

namespace Gascoigne
{

/*---------------------------------------------------*/

void Set2Vec(vector<int>& v, const set<int>& h)
{
  v.resize(h.size());
  int j = 0;
  for (set<int>::const_iterator p=h.begin();
       p!=h.end(); p++)
    {
      v[j++] = *p;
    }
}

/*---------------------------------------------------*/

void Vec2Set(set<int>& h, const vector<int>& v)
{
  h.clear();
  for (int i=0; i<v.size(); i++)
    {
      h.insert(v[i]);
    }
}
}
