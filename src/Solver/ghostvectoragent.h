#ifndef  __GhostVectorAgent_h
#define  __GhostVectorAgent_h

#include  <string>
#include  "gascoigne.h"
#include  "ghostvector.h"
#include  "compvector.h"
#include  "stlio.h"

using namespace Gascoigne;

/*-----------------------------------------*/

class GhostVectorAgent : public std::map<GhostVector,GlobalVector*>
{
protected:

  GlobalVector& Get(const GhostVector& g) 
    {
      return operator()(g);
    }

public:

  typedef std::map<GhostVector,GlobalVector*>::const_iterator const_iterator;
  typedef std::map<GhostVector,GlobalVector*>::iterator       iterator;

  GhostVectorAgent() : std::map<GhostVector,GlobalVector*>() {}

  ~GhostVectorAgent() 
    {
      iterator p=begin();
      while(p!=end())
	{ 
	  if(p->second) 
	    {
	      delete p->second; 
	      p->second=NULL;
	    } 
	  p++;
	}
    }

  void Register(std::string g) 
    {
      if(find(g)!=end()) return;
      insert(std::make_pair(g,(GlobalVector*) NULL));
    }

  GlobalVector& operator()(const GhostVector& g) 
    {
      iterator p = find(g);
      if(p==end())
	{
	  std::cerr << "GhostVectorAgent::operator():\tnotfound\n";
	  std::cerr << g << "\t out of \n";
	  std::cerr << *this << "\n";
	  abort();
	}
      GlobalVector* vp = p->second;
      if(vp==NULL) 
	{
	  std::cerr <<  "GhostVectorAgent scheisse " << p->first << std::endl;
	  std::cerr << *this << std::endl;
	}
      assert(vp);
      return *vp;
    }
  
  void GVequ(GhostVector Gx, double d, GhostVector Gy)
  {
    Get(Gx).equ(d,Get(Gy));
  }
  
  void GVzero(GhostVector Gx)
  {
    Get(Gx).zero();
  }

  void GVadd(GhostVector Gx, double d, GhostVector Gy)
  {
    Get(Gx).add(d,Get(Gy));
  }

  void GVsadd(double d0, GhostVector Gx, double d, GhostVector Gy)
  {
    Get(Gx).sadd(d0,d,Get(Gy));
  }
  
  double GVscp(GhostVector Gx, GhostVector Gy)
  {
    return Get(Gx)*Get(Gy);
  }

};


#endif
