#ifndef  __TwinStencil_h
#define  __TwinStencil_h

#include  "columnstencil.h"


namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments TwinStencil

////
////
/////////////////////////////////////////////

class TwinStencil : public ColumnStencil
{
 protected:

 public:

  TwinStencil() : ColumnStencil() {}
  ~TwinStencil() {}

  int  diag(int i)      const { return sstart[i]; } 
  int  half(int i)      const; 
  void memory(int n, int nt);
  void memory(const SparseStructureInterface*);

  // sollte protected sein, wird aber im moment in "constructstructure" von ilu benutzt
  void diagonalfirst();
};
}

#endif
