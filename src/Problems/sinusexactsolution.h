#ifndef  __SinusExactSolution_h
#define  __SinusExactSolution_h

#include  "exactsolution.h"
#include  "stringutil.h"

/*-----------------------------------------*/

class SinusExactSolution : public ExactSolution
{
protected:

  double nx,ny,nz;

public:

  std::string GetName() const {return "Sinus"+Double2String(nx)+Double2String(ny)+Double2String(nz);}

  int GetNcomp() const { return 1;}

  SinusExactSolution(const std::vector<std::string>& args)
    {
      if(args.size()==0)
	{
	  nx = ny = nz = 1.;
	}
      else if(args.size()==2)
	{
	  nx = atof(args[0].c_str());
	  ny = atof(args[1].c_str());
	  nz = 1.;
	}
      else if(args.size()==3)
	{
	  nx = atof(args[0].c_str());
	  ny = atof(args[1].c_str());
	  nz = atof(args[2].c_str());
	}
      else assert(0);
    }

  double operator()(int c, const Vertex2d& v)const {
    return sin(GascoigneMath::pi()*nx*v.x())*sin(GascoigneMath::pi()*ny*v.y());
  }
  double operator()(int c, const Vertex3d& v)const {
    return sin(GascoigneMath::pi()*nx*v.x()) * 
           sin(GascoigneMath::pi()*ny*v.y()) * 
           sin(GascoigneMath::pi()*nz*v.z());
  }

};


#endif
