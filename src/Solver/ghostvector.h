#ifndef  __GhostVector_h
#define  __GhostVector_h

#include  "basicghostvector.h"
#include  "gascoigne.h"
#include  "solverinterface.h"

namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments GhostVector

////  identifier for a GlobalVector storde in Solver
////  identifier for a GlobalVector storde in Solver
////
/////////////////////////////////////////////

class GhostVector : public BasicGhostVector
{
private:

  const SolverInterface* __S;

protected:


public:


//
////  Con(De)structor 
//

  GhostVector() : BasicGhostVector(), __S(NULL) {}
  GhostVector(const std::string& name) : BasicGhostVector(name), __S(NULL) {}
  GhostVector(const std::string& name, const std::string& type) : BasicGhostVector(name,type), __S(NULL) {}
  GhostVector(const SolverInterface* S, const std::string& name) : BasicGhostVector(name), __S(S) {}
  GhostVector(const SolverInterface* S, const std::string& name, const std::string& type) : BasicGhostVector(name,type), __S(S) {}
  GhostVector(const GhostVector& v) : BasicGhostVector(v) {
    SetSolver(v.GetSolver());
  }
  ~GhostVector() {}

  void SetSolver(const SolverInterface* S) {__S=S;}
  const SolverInterface* GetSolver() const {return __S;}

  GlobalVector& GetGlobalVector() const {
    assert(__S);
    __S->GetGV(*this);
  }

  friend std::ostream& operator<<(std::ostream& os, const GhostVector& g) {
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
  void equ(double d, const GhostVector& gw) {
    GlobalVector& v = GetSolver()->GetGV(*this);
    GlobalVector& w = gw.GetGlobalVector();
    v.equ(d,w);
  }
  void add(double d1, const GhostVector& gw) {
    GlobalVector& v = GetSolver()->GetGV(*this);
    GlobalVector& w = gw.GetGlobalVector();
    v.add(d1,w);
  }
  void sadd(double d1, double d2, const GhostVector& gw) {
    GlobalVector& v = GetSolver()->GetGV(*this);
    GlobalVector& w = gw.GetGlobalVector();
    v.sadd(d1,d2,w);
  }
};
}

#endif
