/*----------------------------   fsi.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __fsi_H
#define __fsi_H
/*----------------------------   fsi.h     ---------------------------*/

#include  "equation.h"
#include  "boundaryequation.h"
#include  "paramfile.h"
#include  "lpsequation.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  class VPSEQ : public LpsEquation
  {

  protected:

    mutable double    __h;
    mutable Vertex2d  __v;
    mutable FemFunction  *OLD;
    double mu_e,mu_v,lambda,lps0;
    mutable double lps;

    //not sure what to do here
    void SetFemData(FemData& q) const
    {
      assert(q.find("old")!=q.end());
      OLD = &q["old"];
    }

  public:
    ~VPSEQ() { }
    VPSEQ() { abort(); }
    VPSEQ(const ParamFile* pf);

    std::string GetName() const { return "VPSEQ"; }
    int    GetNcomp  () const { return 6; }

    void point(double h, const FemFunction& U, const Vertex<2>& v) const;
    void point_M(int j, const FemFunction& U, const TestFunction& M) const;

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const;

    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const;
    void MatrixBlock(EntryMatrix& A, const FemFunction& U,  const FemFunction& NNN) const;

  // LPS
    void lpspoint(double h, const FemFunction& U, const Vertex<2>& v) const;
    void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const;
    void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const;

  };

}

/*----------------------------   fsi.h     ---------------------------*/
/* end of #ifndef __fsi_H */
#endif
/*----------------------------   fsi.h     ---------------------------*/
