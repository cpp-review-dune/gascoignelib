#ifndef  __TwinStencil_h
#define  __TwinStencil_h


/////////////////////////////////////////////
////
////@brief
////  ... comments TwinStencil

////
////
/////////////////////////////////////////////

#include  "columnstencil.h"

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


#endif
