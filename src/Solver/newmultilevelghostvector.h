#ifndef  __NewMultiLevelGhostVector_h
#define  __NewMultiLevelGhostVector_h


/////////////////////////////////////////////
////
////@brief
////  ... comments NewMultiLevelGhostVector

////
////
/////////////////////////////////////////////

#include  "newghostvector.h"
#include  "basicghostvector.h"
#include  "multilevelsolverinterface.h"

class NewMultiLevelGhostVector : public BasicGhostVector
{
private:

  const MultiLevelSolverInterface* __S;

protected:


public:


//
////  Con(De)structor 
//

  NewMultiLevelGhostVector() : __S(NULL), BasicGhostVector() {}
  NewMultiLevelGhostVector(const std::string& name) : __S(NULL), BasicGhostVector(name) {}
  NewMultiLevelGhostVector(const std::string& name, const std::string& type) : __S(NULL), BasicGhostVector(name,type) {}
  NewMultiLevelGhostVector(const std::string& name, const std::string& type, const MultiLevelSolverInterface* S) : __S(S), BasicGhostVector(name,type) {}
  NewMultiLevelGhostVector(const NewMultiLevelGhostVector& v) : BasicGhostVector(v) {
    SetMultiLevelSolver(v.GetMultiLevelSolver());
  }
  ~NewMultiLevelGhostVector() {}

  void SetMultiLevelSolver(const MultiLevelSolverInterface* S) {__S=S;}
  const MultiLevelSolverInterface* GetMultiLevelSolver() const {return __S;}

  int n() const {assert(__S); return __S->nlevels();}

  const Gascoigne::GlobalVector& Vector(int l) const {
    assert(__S);
    return __S->GetSolver(l)->GetGV(*this);
  }

  Gascoigne::GlobalVector& Vector(int l) {
    assert(__S);
    return __S->GetSolver(l)->GetGV(*this);
  }

  const BasicGhostVector& operator()(int l) const {
    return *this;
  }
  BasicGhostVector& operator()(int l) {
    return *this;
  }
  const BasicGhostVector& finest() const {
    return *this;
  }
  BasicGhostVector& finest() {
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const NewMultiLevelGhostVector& g) {
    os << "size:\t" << g.n() << std::endl;
    os << dynamic_cast<const BasicGhostVector&>(g);
    return os;
  }

  void zero() {
    for(int l=0;l<n();l++) {
      Vector(l).zero();
    }
  }
  void equ(double d, const NewMultiLevelGhostVector& v) {
    assert(GetMultiLevelSolver()==v.GetMultiLevelSolver());
    for(int l=0;l<n();l++) {
      Vector(l).equ(d,v.Vector(l));
    }
  }
};


#endif
