#include "hang.h"

using namespace std;

/*********************************************************************/

Hang::Hang() : fixarray<3,int>(-1)
{}

/*********************************************************************/

Hang::Hang(const Hang& h) 
  : fixarray<3,int>(h) {}

/*********************************************************************/

Hang::Hang(int nh, int nr, int nc)
{
  hanging() = nh;
  rneighbour() = nr;
  cneighbour() = nc;
}
  
/*********************************************************************/

ostream& operator<<(ostream &s, const Hang& A)
{
  s << A.hanging()    << " : ";
  s << A.rneighbour() << " ";
  s << A.cneighbour() << " ";
  return s;
}

/*********************************************************************/

istream& operator>>(istream &s, Hang& A)
{
  char symbol;
  s >> A.hanging() >> symbol;
  s >> A.rneighbour() >> A.cneighbour();
  return s;
}
