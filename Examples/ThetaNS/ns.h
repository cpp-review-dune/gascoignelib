/*----------------------------   ns.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __ns_H
#define __ns_H
/*----------------------------   ns.h     ---------------------------*/

#include  "equation.h"
#include  "lpsequation.h"
#include  "boundaryequation.h"
#include  "paramfile.h"
#include <cstdlib>

/*-----------------------------------------*/

namespace Gascoigne
{
  class NS : public Equation
  {
    
  protected:
    
    // problem parameter
    double __nu, __lpsp0,__lpsv0;
    
    mutable double __lpsp,__lpsv, __h;

    mutable FemFunction* OLD;

          
    void SetFemData(FemData& q) const
    {
      assert(q.find("old")!=q.end());
      OLD = &q["old"];
    }

    
  public:
    ~NS() { }
    NS() { abort(); }
    NS(const ParamFile* pf);
    
    
    std::string GetName() const { return "NS"; }
    
    int    GetNcomp  () const { return 3; }
    
    void point(double h, const FemFunction& U, const Vertex2d& v) const;
    void point_M(int j, const FemFunction& U, const TestFunction& M) const;
    
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
    void MatrixBlock(EntryMatrix& A, const FemFunction& U,  const FemFunction& NNN) const;


    //////////
    void lpspoint(double h, const FemFunction& U, const Vertex2d& v) const {
      
    }
    void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
    void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const;    


    
  };
  
  
  class AdjointNS : public Equation
  {
    
  protected:
    
    // problem parameter
    double __nu, __lpsp0,__lpsv0;
    
    mutable double __lpsp,__lpsv,__h;

    mutable FemFunction *ZOLD, *_U;

          
    void SetFemData(FemData& q) const
    {
      assert(q.find("zold")!=q.end());
      ZOLD = &q["zold"];
      assert(q.find("uu")!=q.end());
      _U = &q["uu"];
    }

    
  public:
    ~AdjointNS() { }
    AdjointNS() { abort(); }
    AdjointNS(const ParamFile* pf);
    
    
    std::string GetName() const { return "AdjointNS"; }
    
    int    GetNcomp  () const { return 3; }
    
    void point(double h, const FemFunction& U, const Vertex2d& v) const;
    void point_M(int j, const FemFunction& U, const TestFunction& M) const;
    
    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;
    
    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
    void MatrixBlock(EntryMatrix& A, const FemFunction& U,  const FemFunction& NNN) const;

        //////////
    void lpspoint(double h, const FemFunction& U, const Vertex2d& v) const {
    }
    void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
    
    void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const;    
    

  };
  

}



/*----------------------------   ns.h     ---------------------------*/
/* end of #ifndef __ns_H */
#endif
/*----------------------------   ns.h     ---------------------------*/
