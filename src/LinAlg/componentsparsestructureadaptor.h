#ifndef  __ComponentSparseStructureAdaptor_h
#define  __ComponentSparseStructureAdaptor_h


/////////////////////////////////////////////
///
///@brief
///  ... comments ComponentSparseStructureAdaptor

///
///
/////////////////////////////////////////////


#include  "nodesparsestructureadaptor.h"


class ComponentSparseStructureAdaptor : public NodeSparseStructureAdaptor
{
public:


private:


protected:


public:


  ComponentSparseStructureAdaptor(int ncomp) : NodeSparseStructureAdaptor(ncomp) {}

  int index(int i, int c) const {return i+c*n_base();}

  std::string GetName() const {return "Component";}

};


#endif
