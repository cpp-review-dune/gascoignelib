#ifndef __base3d_h
#define __base3d_h

#include  "base.h"

/**************************************************/

class Base3d : public Base
{
 protected:

  mutable nvector<double>             N;
  mutable std::vector<Vertex3d>       DN;

  mutable Vertex3d  bn, bt;

 public:
  
  Base3d()  {}
  const Vertex3d&  normal () const {return bn;}
  const Vertex3d&  tangent() const {return bt;}

  void point_boundary(int ie, const Vertex2d& s1) const
    {
      Vertex3d s;
      if     (ie==0)      
	{
	   s.x() = s1.x(); s.y() = s1.y(); s.z() = 0.;
	  bn.x() = 0.;    bn.y() = 0.;    bn.z() = -1.;
	}
      else if(ie==1)      
	{
	   s.x() = 1.;  s.y() = s1.y(); s.z() = s1.x();
	  bn.x() = 1.; bn.y() = 0.;    bn.z() = 0.;
	}
      else if(ie==2)      
	{
	   s.x() = s1.x(); s.y() = 1.;  s.z() = s1.y();
	  bn.x() = 0.;    bn.y() = 1.; bn.z() = 0.;
	}
      else if(ie==3)
	{
	   s.x() = 0.; s.y() = s1.y(); s.z() = 1.-s1.x();
	  bn.x() = -1.; bn.y() =  0.; bn.z() = 0.;
	}
      else if(ie==4)      
	{
	   s.x() = s1.x(); s.y() = 0.;   s.z() = 1.-s1.y();
	  bn.x() = 0.;    bn.y() = -1.; bn.z() = 0.;
	}
      else
	{
	   s.x() = 1.-s1.x(); s.y() = s1.y(); s.z() = 1.;
	  bn.x() = 0.;       bn.y() = 0.;    bn.z() = 1.;
	}
      point(s);
    }
};

#endif
