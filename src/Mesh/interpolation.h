#ifndef  __interpolation_h
#define  __interpolation_h

#include "interpollist.h"
#include "solution.h"

class Interpolation
{
  InterpolationList<2> interpol2;
  InterpolationList<4> interpol4;
  InterpolationList<8> interpol8;

 public:

  virtual ~Interpolation(){}

  InterpolationList<2>& list2() { return interpol2; }
  InterpolationList<4>& list4() { return interpol4; }
  InterpolationList<8>& list8() { return interpol8; }
  const InterpolationList<2>& list2() const { return interpol2; }
  const InterpolationList<4>& list4() const { return interpol4; }
  const InterpolationList<8>& list8() const { return interpol8; }

  void clear() {interpol2.clear(); interpol4.clear(); interpol8.clear();}
  virtual void interpolate  (Solution& S) const
  {
    for(int i=0; i<list2().size(); i++)
      {
        int i0 = list2().newvertex(i);
        int i1 = list2().oldvertex(i,0); int i2 = list2().oldvertex(i,1);
        
        S.equ(i0, 0.5,i1, 0.5,i2);
      }  
    for(int i=0; i<list4().size(); i++)
      {
        int i0 = list4().newvertex(i);
        int i1 = list4().oldvertex(i,0);
        int i2 = list4().oldvertex(i,1);
        int i3 = list4().oldvertex(i,2);
        int i4 = list4().oldvertex(i,3);
        
        S.equ(i0, 0.25,i1, 0.25,i2, 0.25,i3, 0.25,i4);
      }  
    for(int i=0; i<list8().size(); i++)
      {
        int i0 = list8().newvertex(i);
        int i1 = list8().oldvertex(i,0);
        int i2 = list8().oldvertex(i,1);
        int i3 = list8().oldvertex(i,2);
        int i4 = list8().oldvertex(i,3);
        int i5 = list8().oldvertex(i,4);
        int i6 = list8().oldvertex(i,5);
        int i7 = list8().oldvertex(i,6);
        int i8 = list8().oldvertex(i,7);
        
        S.equ(i0, 0.125,i1, i2, i3, i4, i5, i6, i7, i8);
      }
  }
};

#endif
