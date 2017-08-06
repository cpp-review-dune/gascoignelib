#include "chi.h"

namespace Gascoigne
{
  
  
  int Chi::operator()(const Vertex3d& v) const
  {
    std::string __type = "bench3d";
    
    double eps = 1.e-7;
    

  

    if (__type == "bench3d")
      {
	if ((v.x()>0.5) && (v.x()<0.9-eps) &&
	    (v.y()>0.19+eps) && (v.y()<0.21-eps) &&
	    (v.z()>0.1+eps)  && (v.z()<0.3 - eps) ) return 1;
	if ((v.x()>0.5) && (v.x()<0.9+eps) &&
	    (v.y()>0.19-eps) && (v.y()<0.21 + eps) &&
	    (v.z()>0.1 - eps)  && (v.z()<0.3 + eps) ) return 0;
	return -1;
      }
    else 
      {
	std::cerr << "solid type " << __type << " not known!" << std::endl;
    	abort();
      }
  }
  
  int Chi::operator()(const Vertex2d& v) const
  {
    
    double eps = 1.e-12;
    std::string __type = "benchmark";
    //    __type=="wall";
    
    if (__type=="benchmark")                     // BENCHMARK
      {
	if ((v.x()>0.2) && (v.x()<0.6-eps) && (fabs(v.y()-0.2)<0.01-eps)) return 1;
	if ((v.x()>0.2) && (v.x()<0.6+eps) && (fabs(v.y()-0.2)<0.01+eps)) return 0;
	return -1;
      }
    else if (__type=="wall")
      {
	if ((v.x()>5.) && (v.x()<9.5-eps) && (fabs(v.y()-6)<0.03-eps)) return 1;
	if ((v.x()>5.) && (v.x()<9.5+eps) && (fabs(v.y()-6)<0.01+eps)) return 0;
	return -1;
      }
    else if (__type=="art")
      {
	if (v.y()<0.5-eps) return 1;
	if (v.y()>1.5+eps) return 1;

	if (v.y()<0.5+eps) return 0;
	if (v.y()>1.5-eps) return 0;

	return -1;
      }
    else if (__type=="ch")
      {
	if ((v.x()>-0.5+eps)&&(v.x()<0.5-eps)) return 1;
	if ((v.x()>-0.5-eps)&&(v.x()<0.5+eps)) return 0;
	return -1;
      }
    else if (__type=="squares")
      {
	if ((v.x()>0.82+eps)&&(v.x()<0.956667-eps)&&(v.y()<0.136667-eps)) return 1;
	if ((v.x()>0.82-eps)&&(v.x()<0.956667+eps)&&(v.y()<0.136667+eps)) return 0;
	return -1;
      }
    else if (__type=="moving") return -1;
    else 
      {
	std::cerr << "Solid type: " << __type << " not known! " << std::endl;
	abort();
      }
    
  }
  
  
  
}


