#include  <iomanip>
#include  <vector>
#include  <algorithm>

#include  "giota.h"

#include  "timer.h"
#include  "compareclass.h"

#define TABWIDTH 23

using namespace std;

/*-----------------------------------------*/

Timer::Timer()
{
}

/*-----------------------------------------*/

ostream& operator<<(ostream& os, const Timer& T)
{
  StopWatch st = T.total();
  double tt = st.read();
  Timer::const_iterator p = T.begin();
  os.precision(2);
  os.setf(ios::fixed, ios::floatfield);
  os << " ------------------------------------\n";
  os << " --  Timing  ------------------------\n";
  os << " ------------------------------------\n";
  vector<double>  x; 
  vector<string>  s; 
  while(p!=T.end())
    {
      s.push_back(p->first);
      x.push_back(p->second.read());
      p++;
    }
  nvector<int> C(x.size()); 
  iota(C.begin(),C.end(),0);
  sort(C.begin(),C.end(),CompareObjectBigToSmall<vector<double> > (x));
  
  for (int i=0; i<x.size(); i++)
    {
      os.setf(ios::left);
      int l = s[C[i]].size();
      os << setw(TABWIDTH-l) << s[C[i]] <<"  ";
      os << T.Get(s[C[i]]).GetTime() << "  " << static_cast<int>(100.*x[C[i]]/tt) <<" \%"<<endl;
    }
  os << " ------------------------------------\n";
  os << "   Total       : " << st.GetTime() << endl << endl;
  return os;
} 
