#ifndef  __PointMatrix_h
#define  __PointMatrix_h



/////////////////////////////////////////////
///
///@brief
///  ... comments PointMatrix

///
///
/////////////////////////////////////////////


#include  "matrixinterface.h"
#include  "simplematrix.h"
#include  "sparsestructureadaptor.h"
//#include  "mginterpolatormatrix.h"


class PointMatrix : public SimpleMatrix, virtual public MatrixInterface
{
public:


private:


protected:

  int _ncomp;
  SparseStructureAdaptor* SSAP;

public:


//
///  Constructor 
//
    PointMatrix(int ncomp, std::string type);
    virtual ~PointMatrix();

    std::string GetName() const {return "PointMatrix";}

    void zero() {
      SimpleMatrix::zero();
    }
    void vmult(CompVector<double>& y, const CompVector<double>& x, double d=1.) const;
    void vmult_transpose(CompVector<double>& y, const CompVector<double>& x, double d=1.) const;

    const StencilInterface* GetStencil() const { return SimpleMatrix::GetStencil();}
    void ReInit(const SparseStructureInterface* S);

    void entry(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.);
    void entry_diag(int i, const nmatrix<double>& M);
    void dirichlet (int i, const std::vector<int>& cv);

    void transpose() {
      SimpleMatrix::transpose();
    }

    void AddMassWithDifferentStencil(const MatrixInterface* M, const TimePattern& TP, double s=1.);

    //    void RestrictMatrix(const MgInterpolatorMatrix& I, const PointMatrix& Ah);
};


#endif
