#include  "vecalgo.h"
#include  "gascoignemath.h"
#include  <algorithm>
#include  "giota.h"

using namespace std;

namespace Gascoigne
{

/*************************************************************/

void transfer(int n, vector<int>& tr, const set<int>& del)
{
  tr.resize(n,-1);

  int count = 0;
  for(int i=0;i<n;++i)
    {
      if(del.find(i)==del.end()) tr[i] = count++;
    }
}

/*************************************************************/

void transfer(int n, vector<int>& tr, vector<int>& del)
{
  tr.resize(n,-1);
  if (del.size()==0)
    {
      iota(tr.begin(),tr.end(),0);
      return;
    }

  sort(del.begin(),del.end());
  //unique(del.begin(),del.end());

  int count = 0;
  int pos   = 0;

  for(int i=0;i<n;++i)
    {
      while ((del[pos]<i) && (pos<del.size()-1))
	{
	  pos++;
	}
      if (pos==del.size())
	{
	  tr[i] = count++;
	}
      else if (del[pos]!=i)
	{
	  tr[i] = count++;
	}
    }
}
}
