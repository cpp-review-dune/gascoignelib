#ifndef  __IntegrationFormulaBase_h
#define  __IntegrationFormulaBase_h

#include  "integrationformulainterface.h"
#include  "gascoigne.h"

/*-----------------------------------------*/


namespace Gascoigne
{
template<int DIM>
class IntegrationFormulaBase : public IntegrationFormulaInterface
{
private:

  typedef Vertex<DIM>   VERTEX;

  int                 _in;
  DoubleVector     _iw;
  std::vector<VERTEX> _ic;

protected:

  void ReInit(int n) {
    _in = n;
    _iw.reserve(n);      
    _iw.resize (n);
    _ic.reserve(n);      
    _ic.resize (n);
  }

public:

  IntegrationFormulaBase<DIM>() : IntegrationFormulaInterface() {}
  IntegrationFormulaBase<DIM>(int n) {
    ReInit(n);
  } 
  IntegrationFormulaBase<DIM>(const IntegrationFormulaBase<DIM>& IF) 
    : _iw(IF.w()), _in(IF.n()), _ic(IF.c()) {}
 
  int    n()                 const { return _in;}
  double w(int k)            const { return _iw[k];}
  const DoubleVector& w() const { return _iw;}

  double& w(int k) { return _iw[k];}
  VERTEX& c(int k) { return _ic[k];}

  void xi(VERTEX& v, int k)  const { v = _ic[k];}
  const std::vector<VERTEX>& c()  const { return _ic;}

};
}

#endif
