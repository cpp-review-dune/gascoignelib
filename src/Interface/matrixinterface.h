#ifndef __matrixinterface_h
#define __matrixinterface_h

#include  "entrymatrix.h"
#include  "compvector.h"
#include  "sparsestructureinterface.h"
#include  "stencilinterface.h"
#include  "timepattern.h"
#include  <string>

/*-------------------------------------------------------------*/

class MatrixInterface
{
 public:

  MatrixInterface() { }
  virtual ~MatrixInterface() { }

  virtual std::string GetName() const=0;
  
  virtual const StencilInterface* GetStencil() const { assert(0); return NULL;}

  virtual void ReInit(const SparseStructureInterface* S)=0;

  virtual void AddMassWithDifferentStencil(const MatrixInterface* M, const TimePattern& TP, double s=1.) { assert(0);}

  virtual void zero()=0;
  virtual void transpose() { assert(0); }

  virtual std::ostream& Write(std::ostream& os) const=0;

  //
  /// for matrix assembling
  //
  typedef nvector<int>::const_iterator niiterator;

  virtual void entry(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.) { assert(0);}

  virtual void entrydual(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.) { assert(0);}

  //
  /// for hanging nodes
  //
  virtual void entry_diag(int i, const nmatrix<double>& M)    { assert(0);}
  //
  /// for boundary conditions
  //
  virtual void dirichlet (int i, const std::vector<int>& cv)  { assert(0);}

  //
  ///
  //
  virtual void vmult(CompVector<double>& y, const CompVector<double>& x, double s=1.) const  { assert(0);}
  virtual void vmult_transpose(CompVector<double>& y, const CompVector<double>& x, double s=1.) const  { assert(0);}
  virtual void vmult_time(CompVector<double>& y, const CompVector<double>& x, const TimePattern& TP, double s=1.) const { assert(0);}
};

/*-------------------------------------------------------------*/

#endif
