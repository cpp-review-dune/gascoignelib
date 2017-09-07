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
  class VelEQ : public LpsEquation
  {

  protected:

    mutable double    __h;
    mutable Vertex2d  __v;
    mutable FemFunction  *OLD,*SIGMA;
    double mu_e,mu_v,lambda,lps0;
    mutable double lps;


    void SetFemData(FemData& q) const
    {
      assert(q.find("old")!=q.end());
      OLD = &q["old"];
      assert(q.find("sigma")!=q.end());
      SIGMA = &q["sigma"];
    }

  public:
    ~VelEQ() { }
    VelEQ() { abort(); }
    VelEQ(const ParamFile* pf);


    std::string GetName() const { return "VelEQ"; }
    int    GetNcomp  () const { return 3; }  // TODO to 6

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



  class StressEQ : public LpsEquation
  {

  protected:

    mutable double    __h;
    mutable Vertex2d  __v;
    mutable FemFunction  *V,*SIGMAOLD;
    double lambda,lpsstress0;
    mutable double lpsstress;


    void SetFemData(FemData& q) const
    {
      assert(q.find("V")!=q.end());
      V = &q["V"];
      assert(q.find("sigmaold")!=q.end());
      SIGMAOLD = &q["sigmaold"];
    }

  public:
    ~StressEQ() { }
    StressEQ() { abort(); }
    StressEQ(const ParamFile* pf);


    std::string GetName() const { return "StressEQ"; }
    int    GetNcomp  () const { return 3; }

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
