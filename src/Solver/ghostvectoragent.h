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

  friend std::ostream& operator<<(std::ostream& os, const GhostVectorAgent& gva) {
    int i=0,n=gva.size();
    os << "GhostVectorAgent: size=" << n << ", ";
    for (const_iterator p=gva.begin(); p!=gva.end(); p++,i++){
      os << "VectorInterface("<<i<<")=('"<< p->first.GetName() << "',"<< p->second <<")";
      if( i <n-1 ) os << ", "; else os << ". ";
    }
    return os;
  }

};
}

#endif
