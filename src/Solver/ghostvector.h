#ifndef  __GhostVector_h
#define  __GhostVector_h


/////////////////////////////////////////////
////
////@brief
////  ... comments GhostVector

////  identifier for a GlobalVector storde in Solver
////  identifier for a GlobalVector storde in Solver
////
/////////////////////////////////////////////

#include  "basicghostvector.h"
#include  "gascoigne.h"
#include  "solverinterface.h"

class GhostVector : public BasicGhostVector
{
private:

  const SolverInterface* __S;

protected:


public:


//
////  Con(De)structor 
//

  GhostVector() : __S(NULL), BasicGhostVector() {}
  GhostVector(const std::string& name) : __S(NULL), BasicGhostVector(name) {}
  GhostVector(const std::string& name, const std::string& type) : __S(NULL), BasicGhostVector(name,type) {}
  GhostVector(const SolverInterface* S, const std::string& name) : __S(S), BasicGhostVector(name) {}
  GhostVector(const SolverInterface* S, const std::string& name, const std::string& type) : __S(S), BasicGhostVector(name,type) {}
  GhostVector(const GhostVector& v) : BasicGhostVector(v) {
    SetSolver(v.GetSolver());
  }
  ~GhostVector() {}

  void SetSolver(const SolverInterface* S) {__S=S;}
  const SolverInterface* GetSolver() const {return __S;}

  Gascoigne::GlobalVector& GetGlobalVector() const {
    assert(__S);
    __S->GetGV(*this);
  }

  friend std::ostream& operator<<(std::ostream& os, const GhostVector& g) {
    os << "Solver:\t" << g.GetSolver() << std::endl;
    os << dynamic_cast<const BasicGhostVector&>(g);
    return os;
  }

  double norm() {
    Gascoigne::GlobalVector& v = GetSolver()->GetGV(*this);
    return v.norm();
  }
  void zero() {
    Gascoigne::GlobalVector& v = GetSolver()->GetGV(*this);
    v.zero();
  }
  void equ(double d, const GhostVector& gw) {
    Gascoigne::GlobalVector& v = GetSolver()->GetGV(*this);
    Gascoigne::GlobalVector& w = gw.GetGlobalVector();
    v.equ(d,w);
  }
  void add(double d1, const GhostVector& gw) {
    Gascoigne::GlobalVector& v = GetSolver()->GetGV(*this);
    Gascoigne::GlobalVector& w = gw.GetGlobalVector();
    v.add(d1,w);
  }
  void sadd(double d1, double d2, const GhostVector& gw) {
    Gascoigne::GlobalVector& v = GetSolver()->GetGV(*this);
    Gascoigne::GlobalVector& w = gw.GetGlobalVector();
    v.sadd(d1,d2,w);
  }
};


#endif
