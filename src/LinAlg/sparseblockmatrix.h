#ifndef __sparseblockmatrix_h
#define __sparseblockmatrix_h

#include  "sparsestructure.h"
#include  "columndiagstencil.h"
#include  "matrixinterface.h"
#include  "gascoigne.h"

using namespace std;

/*-------------------------------------------------------------*/

namespace Gascoigne
{
template<class B>
class SparseBlockMatrix : public MatrixInterface
{
 protected:

  typedef typename vector<B>::const_iterator      const_iterator;
  typedef typename vector<B>::iterator            iterator;
  typedef typename std::pair<int,int>                  IntPair;
  
  ColumnDiagStencil  US;
  vector<B>      smat;
  int            nc;
  
  void matrix_vector_trans(int p, double* yp, const double* xp, double s=1.) const;
  
  int size() const { return smat.size();}

 public:

  SparseBlockMatrix<B>();
  SparseBlockMatrix<B>(const SparseBlockMatrix<B>& A);
  virtual ~SparseBlockMatrix<B>() {}
  
  void transpose();

  string GetName() const {return "SparseBlockMatrix";}

  /////// Zugriff //////////////////////

  const_iterator  mat(int pos)            const { return smat.begin()+pos; }
        iterator  mat(int pos)                  { return smat.begin()+pos; }

  const StencilInterface* GetStencil() const { return &US;}

  int   n()          const { return US.n();};
  int   nentries()   const { return US.nentries();};
  int   ntotal()     const { return smat.size();};

  int  rowsize(int i)     const { return US.start(i+1)-US.start(i);}
  const vector<B>& mat()  const { return smat; }

  ///// Methods //////////////////////

  void AddMassWithDifferentStencil(const MatrixInterface* M, 
				   const TimePattern& TP, double s=1.);
  void copy_entries(const MatrixInterface& S);

  SparseBlockMatrix& operator=(const SparseBlockMatrix<B>& S); 

  void ReInit   (const SparseStructureInterface*);
  void dirichlet(int i, const vector<int>& cv);

/*-------------------------------------------------------------*/
  
/*   void FillInterfaceList(const nvector<int>&, nvector<int>&, nvector<float>& ) const; */
  
/*-------------------------------------------------------------*/

  void zero();
  void entry(nvector<int>::const_iterator start, nvector<int>::const_iterator stop, const EntryMatrix& M, double s=1.);
  void entrydual(nvector<int>::const_iterator start, nvector<int>::const_iterator stop, const EntryMatrix& M, double s=1.);

  void GaussSeidel      (GlobalVector& y, const GlobalVector& x) const;
  void Jacobi           (GlobalVector& x) const;

  void vmult(GlobalVector& y, const GlobalVector& x, double s=1.) const;
  void vmult(GlobalVector& y, const GlobalVector& x, const TimePattern& TP, double s=1.)const;
  void entry_diag(int i, const nmatrix<double>& M);
 
/*-----------------------------------------------*/

  ostream& Write(ostream &s) const;
  friend   ostream& operator<<(ostream &s, const SparseBlockMatrix<B>& A) {}
};
}

#endif