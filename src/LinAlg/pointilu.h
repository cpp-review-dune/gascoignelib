#ifndef  __PointIlu_h
#define  __PointIlu_h

/////////////////////////////////////////////
///
///@brief
///  ... comments PointIlu

///
///
/////////////////////////////////////////////

#include  "iluinterface.h"
#include  "simpleilu.h"
#include  "sparsestructureadaptor.h"

class PointIlu : public SimpleIlu, virtual public IluInterface
{
protected:

  int _ncomp;
  SparseStructureAdaptor* SSAP;

public:

//
///  Constructor 
//
    
    PointIlu(int ncomp, std::string type);
    ~PointIlu();

    std::string GetName() const {return "PointIlu";}
    
    void ReInit(const SparseStructureInterface* S);
    
    void zero() {
      SimpleIlu::zero();
    }
    
    void ConstructStructure(const nvector<int>& perm, const MatrixInterface& A);
    void modify(int c, double s);
    void copy_entries(const MatrixInterface*  A) {
      SimpleIlu::copy_entries(A);
    }
    void compute_ilu() {
      SimpleIlu::compute_ilu();
    }
    void solve(CompVector<double>& x) const {
      SimpleIlu::solve(x);
    }
    void solve_transpose(CompVector<double>& x) const {
      SimpleIlu::solve_transpose(x);
    }
};


#endif
