#ifndef  __ColumnDiagStencil_h
#define  __ColumnDiagStencil_h

#include  "columnstencil.h"


namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments ColumnDiagStencil

////
////
/////////////////////////////////////////////

class ColumnDiagStencil : public ColumnStencil
{
protected:

  IntVector   sdiag;

public:

//
////  Con(De)structor 
//
  ColumnDiagStencil() : ColumnStencil() {}
  ~ColumnDiagStencil() {}

  const IntVector&  diag() const { return sdiag; }
        IntVector&  diag()       { return sdiag; }
        int&           diag(int i)       { return sdiag[i]; } 
  const int&           diag(int i) const { return sdiag[i]; } 

  void memory(int n, int nt);
  void memory(const SparseStructureInterface*);
};
}

#endif
