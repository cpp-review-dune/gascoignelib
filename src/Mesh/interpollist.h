#ifndef __interpollist_h
#define __interpollist_h

#include  <vector.h>
#include  "fixarray.h" 

/*---------------------------------------------------*/

template <int N>
class InterpolElement : public fixarray<N,int>
{
 public:

  int nv;

  InterpolElement() : fixarray<N,int>(), nv(-1) {}
  InterpolElement(const InterpolElement& i) : fixarray<N,int>(i), 
    nv(i.nv) {}

  InterpolElement(const fixarray<N,int>& f, int n) : 
    fixarray<N,int>(f), nv(n) {}
};

/*---------------------------------------------------*/

template <int N>
class InterpolationList : public std::vector<InterpolElement<N> >
{
  public :
    
  int  newvertex(int i)        const { return (*this)[i].nv; }
  int  oldvertex(int i, int j) const { return (*this)[i][j]; }

  void newentry(int nv, const fixarray<N,int>& w)
    {
      InterpolElement<N> I(w,nv);
      push_back(I);
    }
};

#endif
