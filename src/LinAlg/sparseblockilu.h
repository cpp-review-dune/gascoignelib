#ifndef __sparseblockilu_h
#define __sparseblockilu_h

#include  "sparseblockmatrix.h"
#include  "iluinterface.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{
template<class B>
class SparseBlockIlu: public virtual IluInterface, public SparseBlockMatrix<B>
{
protected:

  nvector<int>          p,q;
  GlobalVector*   yp;

  void backward() const;
  void forward () const;
  virtual void hin(const GlobalVector& x) const;
  virtual void her(GlobalVector& x) const;

        int   n()          const { return US.n();};
  const int&  start(int i) const { return US.start(i); }; 
  const int&  stop (int i) const { return US.stop (i); }; 
  const int&  col(int pos) const { return US.col(pos); };
  const int&  diag(int i)  const { return US.diag(i); }; 

  public:

  SparseBlockIlu<B>();
  SparseBlockIlu<B>(const SparseBlockIlu<B>& I);
  ~SparseBlockIlu();

  string GetName() const {return "SparseBlockIlu";}
  
  nvector<int>&       GetP() {return p;}
  nvector<int>&       GetQ() {return q;}
  const nvector<int>& GetP() const {return p;}
  const nvector<int>& GetQ() const {return q;}

  void modify(int c, double s);
  void zero() { SparseBlockMatrix<B>::zero(); }

  void compute_ilu ();
  void ReInit      (const SparseStructureInterface* SI);
  void ConstructStructure(const nvector<int>& perm, const MatrixInterface& A);
  void copy_entries(const MatrixInterface* A);
  void solve       (GlobalVector& x) const;
  void solvetrans  (GlobalVector& x) const { assert(0);};
  ostream& Write(ostream &s) const;
};
}

#endif
