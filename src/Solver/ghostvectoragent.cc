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

void GhostVectorAgent::Register(const BasicGhostVector& mg, const SolverInterface* S) 
{
  const std::string& name = mg.GetName();
  GhostVector g(S,name,mg.GetType());
  g.SetNcomp(mg.GetNcomp());
  iterator p = find(mg);
  if(p!=end())
    {
      assert(p->first.GetSolver()==S);
    }
  else
    {
      insert(std::make_pair(g,static_cast<GlobalVector*>(NULL)));
    }
}
  
/*-------------------------------------------------*/

void GhostVectorAgent::Delete(BasicGhostVector& mg) 
{
  iterator p=find(mg);
  if (p!=end())
    {
      delete p->second; 
      erase(p);
    }
}
    
/*-------------------------------------------------*/
  
GlobalVector& GhostVectorAgent::operator()(const GhostVector& g) 
{
  iterator p = find(g);
  if (p==end())
    {
      std::cerr << __FILE__ << ":" << __LINE__;
      std::cerr << ": GhostVectorAgent::operator(): ERROR:\n";
      std::cerr << __FILE__ << ":" << __LINE__ << ": Ghostvector '"<< g;
      std::cerr <<"' not found in list of: '"<< *this <<"'\n";
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
