#include "ghostvectoragent.h"
#include "stlio.h"

namespace Gascoigne
{

/*-------------------------------------------------*/

GhostVectorAgent::GhostVectorAgent() {}

/*-------------------------------------------------*/

GhostVectorAgent::~GhostVectorAgent()
{
  for (iterator p=begin(); p!=end(); p++)
    { 
      if(p->second) 
	{
	  //  Loesche p->first
	  delete p->second; 
	  p->second = NULL;
	} 
    }
}
  
/*-------------------------------------------------*/

void GhostVectorAgent::Register(const VectorInterface& mg) 
{
  iterator p = find(mg);
  if(p==end())
    {
      insert(std::make_pair(mg,static_cast<GlobalVector*>(NULL)));
    }
}
  
/*-------------------------------------------------*/

void GhostVectorAgent::Delete(VectorInterface& mg) 
{
  iterator p=find(mg);
  if (p!=end())
    {
      delete p->second; 
      erase(p);
    }
}
    
/*-------------------------------------------------*/
  
GlobalVector& GhostVectorAgent::operator()(const VectorInterface& g) 
{
  iterator p = find(g);
  if (p==end())
    {
      std::cerr << __FILE__ << ":" << __LINE__;
      std::cerr << ": GhostVectorAgent::operator(): ERROR"<<std::endl;
      std::cerr << __FILE__ << ":" << __LINE__;
      std::cerr << ": Ghostvector '"<< g <<"' not found in list of: "<<std::endl;
      std::cerr << " "<< *this << std::endl;
      abort();
    }
  GlobalVector* vp = p->second;
  if (vp==NULL) 
    {
      std::cerr <<  "GhostVectorAgent  GlobalVector* NULL\t" << p->first;
      std::cerr << "\n" << *this << std::endl;
      abort();
    }
  return *vp;
}

}
