#include  <iomanip>
#include  <vector>
#include  <algorithm>

#include  "giota.h"

#include  "timer.h"
#include  "compareclass.h"

#define TABWIDTH 23

/*-----------------------------------------*/

Timer::Timer()
{
}

/*-----------------------------------------*/

std::ostream& operator<<(std::ostream& os, const Timer& T)
{
  StopWatch st = T.total();
  double tt = st.read();
  Timer::const_iterator p = T.begin();
  os.precision(2);
  os.setf(std::ios::fixed, std::ios::floatfield);
  os << " ------------------------------------\n";
  os << " --  Timing  ------------------------\n";
  os << " ------------------------------------\n";
  std::vector<double>  x; 
  std::vector<std::string>  s; 
  while(p!=T.end())
    {
      s.push_back(p->first);
      x.push_back(p->second.read());
      p++;
    }
  nvector<int> C(x.size()); 
  iota(C.begin(),C.end(),0);
  std::sort(C.begin(),C.end(),CompareObjectBigToSmall<std::vector<double> > (x));
  
  for (int i=0; i<x.size(); i++)
    {
      os.setf(std::ios::left);
      int l = s[C[i]].size();
      os << std::setw(TABWIDTH-l) << s[C[i]] <<"  ";
      os << T.Get(s[C[i]]).GetTime() << "  " << int(100.*x[C[i]]/tt) <<" \%"<<std::endl;
    }
  os << " ------------------------------------\n";
  os << "   Total       : " << st.GetTime() << std::endl << std::endl;
  return os;
} 
