#ifndef __base2d_h
#define __base2d_h

#include  "base.h"

/**************************************************/

class Base2d : public Base
{
 protected:

  mutable nvector<double>             N;
  mutable std::vector<Vertex2d>       DN;

  mutable Vertex2d  bn, bt;

 public:
  
  const Vertex2d&  normal () const { return bn;}
  const Vertex2d&  tangent() const { return bt;}

  void point_boundary(int ie, const Vertex1d& s1) const
    {
      Vertex2d s;
      if     (ie==0)      
	{
	  s.x() = s1.x(); s.y() = 0.;     
	  bn.x() =  0.; bn.y() = -1.;
	  bt.x() =  1.; bt.y() =  0.;
	}
      else if(ie==1)      
	{
	  s.x() = 1.    ; s.y() = s1.x(); 
	  bn.x() =  1.; bn.y() =  0.;
	  bt.x() =  0.; bt.y() =  1.;
	}
      else if(ie==2)      
	{
	  s.x() = s1.x(); s.y() = 1.;     
	  bn.x() =  0.; bn.y() =  1.;
	  bt.x() = -1.; bt.y() =  0.;
	}
      else                
	{
	  s.x() = 0.    ; s.y() = s1.x(); 
	  bn.x() = -1.; bn.y() =  0.;
	  bt.x() =  0.; bt.y() = -1.;
	}
      point(s);
    }
  
};

#endif
