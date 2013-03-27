#include "mesharray.h"

namespace Tsuchimikado
{
  
  mesharray<9,int> getmesharray(int i1,int i2,int i3,int i4,int i5,int i6,int i7,int i8,int i9)
  {
    mesharray<9,int> t;
    t[0]=i1; t[1]=i2; t[2]=i3; t[3]=i4; t[4]=i5; t[5]=i6; t[6]=i7; t[7]=i8; t[8]=i9;
    return t;
  }
  
  mesharray<8,int> getmesharray(int i1,int i2,int i3,int i4,int i5,int i6,int i7,int i8)
  {
    mesharray<8,int> t;
    t[0]=i1; t[1]=i2; t[2]=i3; t[3]=i4; t[4]=i5; t[5]=i6; t[6]=i7; t[7]=i8; 
    return t;
  }
  
  mesharray<6,int> getmesharray(int i1,int i2,int i3,int i4,int i5,int i6)
  {
    mesharray<6,int> t;
    t[0]=i1; t[1]=i2; t[2]=i3; t[3]=i4; t[4]=i5; t[5]=i6; 
    return t;
  }
  
  mesharray<4,int> getmesharray(int i1,int i2,int i3,int i4)
  {
    mesharray<4,int> t;
    t[0]=i1; t[1]=i2; t[2]=i3; t[3]=i4; 
    return t;
  }
  
  mesharray<3,int> getmesharray(int i1,int i2,int i3)
  {
    mesharray<3,int> t;
    t[0]=i1; t[1]=i2; t[2]=i3; 
    return t;
  }
  
  mesharray<2,int> getmesharray(int i1,int i2)
  {
    mesharray<2,int> t;
    t[0]=i1; t[1]=i2; 
    return t;
  }

  template<int N,class T>
  std::ostream& operator<<(std::ostream &s, const mesharray<N,T>& A)
  {
    copy(A.begin(),A.end(),std::ostream_iterator<T>(s," "));
    return s;
  }

  
  template<int N,class T>
  std::istream& operator>>(std::istream &s, mesharray<N,T>& A)
  {
    typename mesharray<N,T>::iterator p = A.begin();
    while(p!=A.end())
      s >> *p++;
    return s;
  }

  template std::ostream& operator<<(std::ostream &s,  const mesharray<1,double>& A);
  template std::ostream& operator<<(std::ostream &s,  const mesharray<2,double>& A);
  template std::ostream& operator<<(std::ostream &s,  const mesharray<3,double>& A);
  template std::ostream& operator<<(std::ostream &s,  const mesharray<4,double>& A);
  template std::ostream& operator<<(std::ostream &s,  const mesharray<6,double>& A);
  template std::ostream& operator<<(std::ostream &s,  const mesharray<8,double>& A);

  template std::istream& operator>>(std::istream &s,  mesharray<2,double>& A);
  template std::istream& operator>>(std::istream &s,  mesharray<3,double>& A);
  template std::istream& operator>>(std::istream &s,  mesharray<4,double>& A);
  template std::istream& operator>>(std::istream &s,  mesharray<6,double>& A);
  template std::istream& operator>>(std::istream &s,  mesharray<8,double>& A);

  template std::ostream& operator<<(std::ostream &s,  const mesharray<2,int>& A);
  template std::ostream& operator<<(std::ostream &s,  const mesharray<3,int>& A);
  template std::ostream& operator<<(std::ostream &s,  const mesharray<4,int>& A);
  template std::ostream& operator<<(std::ostream &s,  const mesharray<6,int>& A);
  template std::ostream& operator<<(std::ostream &s,  const mesharray<8,int>& A);

  template std::istream& operator>>(std::istream &s,  mesharray<2,int>& A);
  template std::istream& operator>>(std::istream &s,  mesharray<3,int>& A);
  template std::istream& operator>>(std::istream &s,  mesharray<4,int>& A);
  template std::istream& operator>>(std::istream &s,  mesharray<6,int>& A);
  template std::istream& operator>>(std::istream &s,  mesharray<8,int>& A);
  
}
