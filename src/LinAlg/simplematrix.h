#ifndef  __SimpleMatrix_h
#define  __SimpleMatrix_h



/////////////////////////////////////////////
///
///@brief
///  ... comments SimpleMatrix

///
///
/////////////////////////////////////////////


#include  "matrixinterface.h"
#include  "columndiagstencil.h"
#include  "sparsestructureadaptor.h"
#include  "compvector.h"


namespace Gascoigne
{
class SimpleMatrix : virtual public MatrixInterface
{
protected:

  ColumnDiagStencil  ST;
  DoubleVector value; 

public:

//
///  Constructor 
//
    SimpleMatrix() : MatrixInterface() {}
    ~SimpleMatrix() {}

    std::string GetName() const {return "SimpleMatrix";}

    std::ostream& Write(std::ostream& os) const;

    const StencilInterface* GetStencil() const { return &ST;}
    double& GetValue(int pos) {return value[pos];}
    const double& GetValue(int pos) const {return value[pos];}
    const double& GetValue(int i, int j) const {return value[ST.Find(i,j)];}

    void zero() {value.zero();}
    void ReInit(const SparseStructureInterface* S);
    void ReInit(int n, int nentries);
    void entry(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.);
    void vmult(DoubleVector& y, const DoubleVector& x, double d=1.) const;
    void vmult_transpose(DoubleVector& y, const DoubleVector& x, double d=1.) const;
    void vmult_comp(int c, int d, GlobalVector& y, const GlobalVector& x, double s=1.) const;
    void vmult_comp_trans(int c, int d, GlobalVector& y, const GlobalVector& x, double s=1.) const;
    void vmult_time(GlobalVector& y, const GlobalVector& x, const TimePattern& TP, double s=1.) const;
    void dirichlet(const IntVector& indices);

    void transpose();
    void entry_diag(int i, const nmatrix<double>& M);
};
}

#endif
