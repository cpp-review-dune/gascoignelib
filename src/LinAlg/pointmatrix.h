#ifndef  __PointMatrix_h
#define  __PointMatrix_h

#include  "matrixinterface.h"
#include  "simplematrix.h"
#include  "sparsestructureadaptor.h"
#include  "mginterpolatormatrix.h"

namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments PointMatrix

///
///
/////////////////////////////////////////////

class PointMatrix : public SimpleMatrix, virtual public MatrixInterface
{
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
    void vmult(GlobalVector& y, const GlobalVector& x, double d=1.) const;
    void vmult_transpose(GlobalVector& y, const GlobalVector& x, double d=1.) const;

    const StencilInterface* GetStencil() const { return SimpleMatrix::GetStencil();}
    void ReInit(const SparseStructureInterface* S);

    void entry(niiterator start, niiterator stop, const EntryMatrix& M, double s=1.);
    void entry_diag(int i, const nmatrix<double>& M);
    void dirichlet (int i, const std::vector<int>& cv);
    void dirichlet_only_row (int i, const std::vector<int>& cv);
    void periodic (const std::map<int,int> &m_PeriodicPairs, const IntVector &iv_Components);

    void transpose() {
      SimpleMatrix::transpose();
    }

    void AddMassWithDifferentStencil(const MatrixInterface* M, const TimePattern& TP, double s=1.);
    void AddMassWithDifferentStencilJacobi(const MatrixInterface* M, const TimePattern& TP, double s=1.);

    void RestrictMatrix(const MgInterpolatorMatrix& I, const PointMatrix& Ah);
};
}

#endif
