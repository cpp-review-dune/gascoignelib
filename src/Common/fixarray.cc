#include  "fixarray.h"

/*------------------------------------------------*/
/*------------------------------------------------*/
/*------------------------------------------------*/

template<int N,class T>
std::ostream& operator<<(std::ostream &s, const fixarray<N,T>& A)
{
  copy(A.begin(),A.end(),std::ostream_iterator<T>(s," "));
  return s;
}


template<int N,class T> std::istream& operator>>(std::istream &s, fixarray<N,T>& A)
{
  typename fixarray<N,T>::iterator p = A.begin();
  while(p!=A.end())
    s >> *p++;
  return s;
}


template std::ostream& operator<<(std::ostream &s, const fixarray<1,float>& A);
template std::ostream& operator<<(std::ostream &s, const fixarray<1,double>& A);
template std::ostream& operator<<(std::ostream &s, const fixarray<2,double>& A);
template std::ostream& operator<<(std::ostream &s, const fixarray<3,double>& A);
template std::ostream& operator<<(std::ostream &s, const fixarray<4,double>& A);
template std::ostream& operator<<(std::ostream &s, const fixarray<6,double>& A);
template std::ostream& operator<<(std::ostream &s, const fixarray<8,double>& A);
template std::ostream& operator<<(std::ostream &s, const fixarray<9,float>& A);
template std::ostream& operator<<(std::ostream &s, const fixarray<16,float>& A);
template std::ostream& operator<<(std::ostream &s, const fixarray<25,float>& A);
template std::ostream& operator<<(std::ostream &s, const fixarray<400,float>& A);

template std::istream& operator>>(std::istream &s,  fixarray<1,float>& A);
template std::istream& operator>>(std::istream &s,  fixarray<1,double>& A);
template std::istream& operator>>(std::istream &s,  fixarray<2,double>& A);
template std::istream& operator>>(std::istream &s,  fixarray<3,double>& A);
template std::istream& operator>>(std::istream &s,  fixarray<4,double>& A);
template std::istream& operator>>(std::istream &s,  fixarray<6,double>& A);
template std::istream& operator>>(std::istream &s,  fixarray<8,double>& A);
template std::istream& operator>>(std::istream &s,  fixarray<9,float>& A);
template std::istream& operator>>(std::istream &s,  fixarray<16,float>& A);
template std::istream& operator>>(std::istream &s,  fixarray<25,float>& A);
template std::istream& operator>>(std::istream &s,  fixarray<400,float>& A);

template std::ostream& operator<<(std::ostream &, const fixarray<1,int>&);
template std::ostream& operator<<(std::ostream &, const fixarray<2,int>&);
template std::ostream& operator<<(std::ostream &, const fixarray<3,int>&);
template std::ostream& operator<<(std::ostream &, const fixarray<4,int>&);
template std::ostream& operator<<(std::ostream &, const fixarray<6,int>&);
template std::ostream& operator<<(std::ostream &, const fixarray<8,int>&);
template std::ostream& operator<<(std::ostream &, const fixarray<9,int>& );

template std::istream& operator>>(std::istream &, fixarray<1,int>& );
template std::istream& operator>>(std::istream &, fixarray<2,int>& );
template std::istream& operator>>(std::istream &, fixarray<3,int>& );
template std::istream& operator>>(std::istream &, fixarray<4,int>& );
template std::istream& operator>>(std::istream &, fixarray<6,int>& );
template std::istream& operator>>(std::istream &, fixarray<8,int>& );
template std::istream& operator>>(std::istream &, fixarray<9,int>& );

