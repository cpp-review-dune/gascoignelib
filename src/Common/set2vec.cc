#include  "set2vec.h"

/*---------------------------------------------------*/

void Set2Vec(std::vector<int>& v, const std::set<int>& h)
{
  v.resize(h.size());
  int j = 0;
  for (std::set<int>::const_iterator p=h.begin();
       p!=h.end(); p++)
    {
      v[j++] = *p;
    }
}

/*---------------------------------------------------*/

void Vec2Set(std::set<int>& h, const std::vector<int>& v)
{
  h.clear();
  for (int i=0; i<v.size(); i++)
    {
      h.insert(v[i]);
    }
}
