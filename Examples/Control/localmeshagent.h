#ifndef  __LocalMeshAgent_h
#define  __LocalMeshAgent_h


/////////////////////////////////////////////
////
////@brief
////  ... comments LocalMeshAgent

////
////
/////////////////////////////////////////////


#include  "meshagent.h"

class LocalMeshAgent : public MeshAgent
{
private:


protected:


public:


//
////  Con(De)structor 
//
  LocalMeshAgent() : MeshAgent() {}
  void BasicInit(const std::string& paramfile);

};


#endif
