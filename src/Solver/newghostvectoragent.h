#ifndef  __NewGhostVectorAgent_h
#define  __NewGhostVectorAgent_h

/////////////////////////////////////////////
////
////@brief
////  ... comments NewGhostVectorAgent

////
////
/////////////////////////////////////////////

#include  <string>
#include  "gascoigne.h"
#include  "newghostvector.h"
#include  "newmultilevelghostvector.h"
#include  "stlio.h"

using namespace std;
using namespace Gascoigne;

class NewGhostVectorAgent : public map<NewGhostVector,GlobalVector*>
{
public:

//
////  Con(De)structor 
//

  typedef std::map<NewGhostVector,GlobalVector*>::const_iterator const_iterator;
  typedef std::map<NewGhostVector,GlobalVector*>::iterator       iterator;

  NewGhostVectorAgent() {}
  ~NewGhostVectorAgent()
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
      const string& name = mg.GetName();
      NewGhostVector g(S,name,mg.GetType());
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
	  insert(std::make_pair(g,(GlobalVector*) NULL));
	}
    }
  void Delete(BasicGhostVector& mg) 
    {
      iterator p=find(mg);
      if(p==end()) return;
      delete p->second; 
      erase(p);
    }


  GlobalVector& operator()(const NewGhostVector& g) 
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


#endif
