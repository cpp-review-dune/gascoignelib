#include "chi.h"

namespace Gascoigne
{
  
  
  int Chi::operator()(const Vertex3d& v) const
  {
    double eps = 1.e-12;
    
    if (__type=="bench3d")
      {
	
	// Benchmark
	if ( (v.x()>0.45+eps)&&(v.x()<0.55-eps)&&(v.y()<0.25-eps)&&(v.y()>0.15+eps)) return 1; 
	if ( (v.x()>0.45-eps)&&(v.x()<0.55+eps)&&(v.y()<0.25+eps)&&(v.y()>0.15-eps)) return 0; 
	return -1; // fluid
      }
    else if (__type=="new")
      {
	if ( (v.x()>0.4+eps)&&(v.x()<0.5-eps)&&
	     (v.y()<0.2-eps)&&(v.z()<0.2-eps)) return 1; 
	if ( (v.x()>0.4-eps)&&(v.x()<0.5+eps)&&
	     (v.y()<0.2+eps)&&(v.z()<0.2+eps)) return 0; 
	return -1; // fluid
      }

    else abort();
  }
  
  int Chi::operator()(const Vertex2d& v) const
  {
    return -1;
    

    double eps = 1.e-12;

    if (__type=="benchmark")                     // BENCHMARK
      {
	if ((v.x()>0.2) && (v.x()<0.6-eps) && (fabs(v.y()-0.2)<0.01-eps)) return 1;
	if ((v.x()>0.2) && (v.x()<0.6+eps) && (fabs(v.y()-0.2)<0.01+eps)) return 0;
	return -1;
      }
    else if (__type=="simple")
      {
	if (v.y()<1.0-eps) return 1;
	if (v.y()<1.0+eps) return 0;
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
    else 
      {
	std::cerr << "Solid type: " << __type << " not known! " << std::endl;
	abort();
      }
    
  }
  
  
  
}


