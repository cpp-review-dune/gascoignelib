#ifndef  __SimpleSparseStructureAdaptor_h
#define  __SimpleSparseStructureAdaptor_h


/////////////////////////////////////////////
///
///@brief
///  ... comments SimpleSparseStructureAdaptor

///
///
/////////////////////////////////////////////

#include  "sparsestructureadaptor.h"

namespace Gascoigne
{
class SimpleSparseStructureAdaptor : public SparseStructureAdaptor
{
public:


private:


protected:


public:

  SimpleSparseStructureAdaptor() : SparseStructureAdaptor() {}

  std::string GetName() const {return "Simple";}

  int n() const {return n_base();} 
  int nentries() const {return nentries_base();} 
  void FillStencil(ColumnDiagStencil& S) const;

  int index(int i, int c) const {return i;}

  IntVector GetIndicesDirichlet(int inode, const std::vector<int>& cv) const{return IntVector(1,inode);}
};
}

#endif
