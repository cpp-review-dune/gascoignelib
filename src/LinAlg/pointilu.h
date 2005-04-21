#ifndef  __PointIlu_h
#define  __PointIlu_h

#include  "iluinterface.h"
#include  "simpleilu.h"
#include  "sparsestructureadaptor.h"


namespace Gascoigne
{

/////////////////////////////////////////////
///
///@brief
///  ... comments PointIlu

///
///
/////////////////////////////////////////////

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
    
    int   n()          const { return GetStencil()->n();};
    void zero() {
      SimpleIlu::zero();
    }
    
    void ConstructStructure(const IntVector& perm, const MatrixInterface& A);
    void modify(int c, double s);
    void copy_entries(const MatrixInterface*  A) {
      SimpleIlu::copy_entries(A);
    }
    void compute_ilu() {
      SimpleIlu::compute_ilu();
    }
    void solve(GlobalVector& x) const {
      SimpleIlu::solve(x);
    }
    void solve_transpose(GlobalVector& x) const {
      SimpleIlu::solve_transpose(x);
    }
};
}

#endif
