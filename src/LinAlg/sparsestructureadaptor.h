#ifndef  __SparseStructureAdaptor_h
#define  __SparseStructureAdaptor_h



/////////////////////////////////////////////
///
///@brief
///  ... comments SparseStructureAdaptor

///
///
/////////////////////////////////////////////


#include  "columndiagstencil.h"
#include  "sparsestructure.h"
#include  <string>

namespace Gascoigne
{
class SparseStructureAdaptor
{
public:


private:

protected:

  int _nnodes;
  const SparseStructure* SSP;

  int n_base() const {assert(SSP);   return SSP->n();} 
  int nentries_base() const {assert(SSP); return SSP->ntotal();} 

public:


//
///  Constructor 
//
    SparseStructureAdaptor() : SSP(NULL) {}
    virtual ~SparseStructureAdaptor() {}

    virtual std::string GetName() const=0;

    void InitStructure(const SparseStructureInterface* SS) {
      SSP = dynamic_cast<const SparseStructure*>(SS);
      assert(SSP);
      _nnodes = SSP->n();
    }

    int nnodes() const {return _nnodes;}
    virtual int index(int i, int c) const=0;

    virtual int n() const=0; 
    virtual int nentries() const=0; 
    virtual void FillStencil(ColumnDiagStencil& S) const=0;
    virtual IntVector GetIndicesDirichlet(int inode, const std::vector<int>& cv)const=0;
};
}

#endif
