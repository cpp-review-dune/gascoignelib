#ifndef __runderkreis_h
#define __runderkreis_h

/* ----------------------------------------- */

class RunderKreis : public BoundaryFunction<2>
{
  double squareradius;
  Vertex2d center;
public :

  string GetName() const { return "RunderKreis";}

  void BasicInit(Vertex2d c, double r) 
    {
      center = c; 
      squareradius = r;
    }
  double operator()(const Vertex2d& c) const 
    {
      double r = - squareradius;
      for (int i=0; i<2; i++)
	{
	  double dx = c[i]-center[i];
	  r += dx * dx;
	}
      return r;
    }
};

#endif
