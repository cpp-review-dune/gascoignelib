#ifndef  __GhostVectorAgent_h
#define  __GhostVectorAgent_h

#include  <string>
#include  "gascoigne.h"
#include  "vectorinterface.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments GhostVectorAgent
////
////
/////////////////////////////////////////////

class GhostVectorAgent : public std::map<VectorInterface,GlobalVector*>
{
public:

  typedef std::map<VectorInterface,GlobalVector*>::const_iterator const_iterator;
  typedef std::map<VectorInterface,GlobalVector*>::iterator       iterator;

//
////  Con(De)structor 
//

  GhostVectorAgent();
  ~GhostVectorAgent();

  void Register(const VectorInterface& mg);
  void Delete(VectorInterface& mg);

  GlobalVector& operator()(const VectorInterface& g);
};
}

#endif
