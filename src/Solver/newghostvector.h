#ifndef  __NewGhostVector_h
#define  __NewGhostVector_h


/////////////////////////////////////////////
////
////@brief
////  ... comments NewGhostVector

////  identifier for a GlobalVector storde in Solver
////  identifier for a GlobalVector storde in Solver
////
/////////////////////////////////////////////

#include  "basicghostvector.h"
#include  "gascoigne.h"
#include  "solverinterface.h"

using namespace Gascoigne;

class NewGhostVector : public BasicGhostVector
{
private:

  const SolverInterface* __S;

protected:


public:


//
////  Con(De)structor 
//

  NewGhostVector() : __S(NULL), BasicGhostVector() {}
  NewGhostVector(const std::string& name) : __S(NULL), BasicGhostVector(name) {}
  NewGhostVector(const std::string& name, const std::string& type) : __S(NULL), BasicGhostVector(name,type) {}
  NewGhostVector(const SolverInterface* S, const std::string& name) : __S(S), BasicGhostVector(name) {}
  NewGhostVector(const SolverInterface* S, const std::string& name, const std::string& type) : __S(S), BasicGhostVector(name,type) {}
  NewGhostVector(const NewGhostVector& v) : BasicGhostVector(v) {
    SetSolver(v.GetSolver());
  }
  ~NewGhostVector() {}

  void SetSolver(const SolverInterface* S) {__S=S;}
  const SolverInterface* GetSolver() const {return __S;}

  GlobalVector& GetGlobalVector() const {
    assert(__S);
    __S->GetGV(*this);
  }

  friend std::ostream& operator<<(std::ostream& os, const NewGhostVector& g) {
    os << "Solver:\t" << g.GetSolver() << std::endl;
    os << dynamic_cast<const BasicGhostVector&>(g);
    return os;
  }

  double norm() {
    GlobalVector& v = GetSolver()->GetGV(*this);
    return v.norm();
  }
  void zero() {
    GlobalVector& v = GetSolver()->GetGV(*this);
    v.zero();
  }
  void equ(double d, const NewGhostVector& gw) {
    GlobalVector& v = GetSolver()->GetGV(*this);
    GlobalVector& w = gw.GetGlobalVector();
    v.equ(d,w);
  }
  void add(double d1, const NewGhostVector& gw) {
    GlobalVector& v = GetSolver()->GetGV(*this);
    GlobalVector& w = gw.GetGlobalVector();
    v.add(d1,w);
  }
  void sadd(double d1, double d2, const NewGhostVector& gw) {
    GlobalVector& v = GetSolver()->GetGV(*this);
    GlobalVector& w = gw.GetGlobalVector();
    v.sadd(d1,d2,w);
  }
};


#endif
