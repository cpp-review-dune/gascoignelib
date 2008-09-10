#ifndef  __BenchmarkMeshAgent_h
#define  __BenchmarkMeshAgent_h

#include  "meshagent.h"

/* ----------------------------------------- */

class RunderKreis : public Gascoigne::BoundaryFunction<2>
{
  double              _r;
  Gascoigne::Vertex2d _c;

public :

  std::string GetName() const { return "RunderKreis";}
  void BasicInit(Gascoigne::Vertex2d c, double r) 
  {
      _c = c;
      _r = r;
    }
  double operator()(const Gascoigne::Vertex2d& c) const 
  {
      double r = - _r;
      for (int i=0; i<2; i++)
        {
          double dx = c[i]-_c[i];
          r += dx * dx;
        }
      return r;
    }
};

/*----------------------------------------------------------------------------*/

class CurvedMeshAgent : public Gascoigne::MeshAgent
{
 protected:

  RunderKreis RK;

 public:

  CurvedMeshAgent() : MeshAgent()
    {
      double r = 0.25;
      Gascoigne::Vertex2d v(2.,2.);
      RK.BasicInit(v,r);

      AddShape(80,&RK);
    }
};

/*-------------------------------------------------------------*/

#endif
