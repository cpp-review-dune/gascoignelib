#ifndef  __GhostVectorAgent_h
#define  __GhostVectorAgent_h

#include  <string>
#include  "gascoigne.h"
#include  "ghostvector.h"


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

  typedef std::map<GhostVector,GlobalVector*>::const_iterator const_iterator;
  typedef std::map<GhostVector,GlobalVector*>::iterator       iterator;

//
////  Con(De)structor 
//

  GhostVectorAgent();
  ~GhostVectorAgent();

  void Register(const BasicGhostVector& mg, const SolverInterface* S);
  void Delete(BasicGhostVector& mg);

  GlobalVector& operator()(const GhostVector& g);
};
}

#endif
