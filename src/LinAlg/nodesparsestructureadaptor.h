#ifndef  __NodeSparseStructureAdaptor_h
#define  __NodeSparseStructureAdaptor_h


/////////////////////////////////////////////
///
///@brief
///  ... comments NodeSparseStructureAdaptor

///
///
/////////////////////////////////////////////


#include  "sparsestructureadaptor.h"


class NodeSparseStructureAdaptor : public SparseStructureAdaptor
{
public:


private:


protected:

  int _ncomp;

public:


  NodeSparseStructureAdaptor(int ncomp) : _ncomp(ncomp), SparseStructureAdaptor() {}

  std::string GetName() const {return "Node";}

  int n() const { return _ncomp*n_base();}
  int nentries() const { return _ncomp*_ncomp*nentries_base();}

  int index(int i, int c) const {return i*_ncomp+c;}

  void FillStencil(ColumnDiagStencil& S) const;
  nvector<int> GetIndicesDirichlet(int inode, const std::vector<int>& cv) const;
};


#endif
