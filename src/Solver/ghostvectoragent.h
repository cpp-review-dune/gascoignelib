#ifndef  __GhostVectorAgent_h
#define  __GhostVectorAgent_h

#include  <string>
#include  "gascoigne.h"
#include  "ghostvector.h"
#include  "stlio.h"


namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments GhostVectorAgent

////
////
/////////////////////////////////////////////

class GhostVectorAgent : public std::map<GhostVector,GlobalVector*>
{
public:

//
////  Con(De)structor 
//

  typedef std::map<GhostVector,GlobalVector*>::const_iterator const_iterator;
  typedef std::map<GhostVector,GlobalVector*>::iterator       iterator;

  GhostVectorAgent() {}
  ~GhostVectorAgent()
    {
      iterator p=begin();
      while(p!=end())
	{ 
// 	  std::cerr << "GhostVectorAgent loesche:\t"<<p->first << endl;
	  if(p->second) 
	    {
	      delete p->second; 
	      p->second=NULL;
	    } 
	  p++;
	}
    }

  void Register(const BasicGhostVector& mg, const SolverInterface* S) 
    {
      const std::string& name = mg.GetName();
      GhostVector g(S,name,mg.GetType());
      iterator p = find(mg);
      if(p!=end())
	{
// 	  std::cerr << "GhostVectorAgent::Register():\talrteady registered\n";
	  assert(p->first.GetSolver()==S);
// 	  std::cerr << g << "\t out of \n";
// 	  std::cerr << *this << "\n";
// 	  abort();
	}
      else
	{
	  insert(std::make_pair(g,static_cast<GlobalVector*>(NULL)));
	}
    }
  void Delete(BasicGhostVector& mg) 
    {
      iterator p=find(mg);
      if(p==end()) return;
      delete p->second; 
      erase(p);
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
	  std::cerr <<  "GhostVectorAgent  GlobalVector* NULL\t" << p->first << std::endl;
	  std::cerr << *this << std::endl;
	  abort();
	}
      assert(vp);
      return *vp;
    }

};
}

#endif
