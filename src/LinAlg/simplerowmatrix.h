#ifndef  __SimpleRowMatrix_h
#define  __SimpleRowMatrix_h

#include  "matrixinterface.h"
#include  "rowcolumnstencil.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments SimpleRowMatrix

////
////
/////////////////////////////////////////////

class SimpleRowMatrix : virtual public MatrixInterface
{
private:
  
  
protected:
  
  RowColumnStencil  ST;
  nvector<double> value; 
  
public:
  
  
  //
  ////  Con(De)structor 
  //
  
  SimpleRowMatrix() : MatrixInterface() {}
  ~SimpleRowMatrix() {}
  
  std::string GetName() const {return "SimpleRowMatrix";}
  
  std::ostream& Write(std::ostream& os) const;
  
  const StencilInterface* GetStencil() const { return &ST;}
  double& GetValue(int pos) {return value[pos];}
  const double& GetValue(int pos) const {return value[pos];}
  const double& GetValue(int i, int j) const {return value[ST.Find(i,j)];}
  
  void zero() {value.zero();}
  void ReInit(int n, int nentries);
  void ReInit(const SparseStructureInterface* S);
  void entry(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.);
  void vmult(nvector<double>& y, const nvector<double>& x, double d=1.) const;
  
};
}

#endif
