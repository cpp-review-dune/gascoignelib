#ifndef  __hang_h
#define  __hang_h

#include  <vector>
#include  "fixarray.h"
#include  <string>

/*------------------------------------------------------*/

namespace Gascoigne
{
class Hang : public fixarray<3,int>
{
 public:

  Hang(); 
  Hang(const Hang& h);
  Hang(int nh, int nr, int nc) ;
  
  int  hanging   () const { return (*this)[0]; }
  int& hanging   ()       { return (*this)[0]; }
  int  rneighbour() const { return (*this)[1]; }
  int& rneighbour()       { return (*this)[1]; }
  int  cneighbour() const { return (*this)[2]; }
  int& cneighbour()       { return (*this)[2]; }

  const fixarray<3,int>& operator()() const { return (*this);}

  friend std::ostream& operator<<(std::ostream &s, const Hang& A);
  friend std::istream& operator>>(std::istream &s, Hang& A);
};
}

#endif
