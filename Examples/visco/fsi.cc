#include  "fsi.h"
#include  "filescanner.h"

extern double __DT;
extern double __THETA;



using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{

  VelEQ::VelEQ(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("mu_v" ,    &mu_v , 0.0);
    DFH.insert("mu_e" ,    &mu_e , 0.0);
    DFH.insert("lambda" ,    &lambda , 0.0);
    DFH.insert("lps" ,    &lps0 , 0.0);
    FileScanner FS(DFH, pf, "Equation");
  }

  void VelEQ::point(double h, const FemFunction& U, const Vertex<2>& v) const
  {
    __h = h;
    __v = v;
  }



  void VelEQ::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    for (int i=0;i<2;++i){
      //p?
      b[0]   += U[i+1][i+1] * N.m();

      //b[i+1] += (U[i+1].m()-(*OLD)[i+1].m())/__DT * N.m();

      for (int j=0;j<2;++j){
        //?
        b[i+1] += mu_v * 1.0 * (U[i+1][j+1] + U[j+1][i+1]) * N[j+1];
        //      b[i+1] += mu_v * 0.5 * ((*OLD)[i+1][j+1] + (*OLD)[j+1][i+1]) * N[j+1];
      }

      b[i+1] -= U[0].m() * N[i+1];
      }

      b[1] += mu_e/lambda * (*SIGMA)[0].m() * N.x();
      b[1] += mu_e/lambda * (*SIGMA)[1].m() * N.y();
      b[2] += mu_e/lambda * (*SIGMA)[1].m() * N.x();
      b[2] += mu_e/lambda * (*SIGMA)[2].m() * N.y();
      b[1] += mu_e/lambda * (-1.0) * N.x();
      b[2] += mu_e/lambda * (-1.0) * N.y();
  }


  void VelEQ::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    for (int i=0;i<2;++i)
      {
  A(0,i+1)  += M[i+1] * N.m();

  //A(i+1,i+1) += M.m()/__DT * N.m();

  for (int j=0;j<2;++j)
    {
      A(i+1,i+1) += mu_v * 1.0 * M[j+1] * N[j+1];
      A(i+1,j+1) += mu_v * 1.0 * M[i+1] * N[j+1];
    }

  A(i+1,0) -= M.m() * N[i+1];
      }
  }




  void VelEQ::point_M(int j, const FemFunction& U, const TestFunction& M) const
  {

  }


  void VelEQ::MatrixBlock(EntryMatrix& A, const FemFunction& U, const FemFunction& N) const
  {
    for (int j=0; j<N.size();++j)  // trial
      {
#define M N[j]
  point_M(j,U,M);


  for (int i=0; i<N.size();++i) // test
    {
      A.SetDofIndex(i,j);
      Matrix(A,U,M , N[i]);
    }
#undef M
      }
  }



  void VelEQ::lpspoint(double h, const FemFunction& U, const Vertex<2>& v) const
  {
    double vel = 1.0;
    lps = lps0 / ( mu_v/h/h + vel/h);
  }

  void VelEQ::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
  {
    for (int i=0;i<2;++i)
      b[0] += lps * UP[0][i+1] * N[i+1];
  }


  void VelEQ::StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
  {
    for (int i=0;i<2;++i)
      A(0,0) += lps *Mp[i+1] * Np[i+1];
  }

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////


  StressEQ::StressEQ(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("lambda" ,    &lambda , 0.0);
    DFH.insert("lpsstress" ,    &lpsstress0 , 0.0);
    FileScanner FS(DFH, pf, "Equation");
  }

  void StressEQ::point(double h, const FemFunction& U, const Vertex<2>& v) const
  {
    __h = h;
    __v = v;
  }



  void StressEQ::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    for (int i=0;i<3;++i)
      {
      //lambda (B - Bp) / dt Bt
      b[i] += lambda * (U[i].m()-(*SIGMAOLD)[i].m())/__DT * N.m();

      for (int j=0;j<2;++j) {
        //lambda v grad(B)
        b[i] += lambda * (*V)[j+1].m()*(*SIGMAOLD)[i][j+1] * N.m();
      //  b[i] += lambda * (*V)[j+1].m()*U[i][j+1] * N.m();
      //  b[i] += (*SIGMAOLD)[i].m() * N.m();
      }
      //B
      b[i] += U[i].m() * N.m();
      }
    // - I
    b[0] += -1.0 * N.m();
    b[2] += -1.0 * N.m();

    // b[0] -= lambda * (2.0 * ((*SIGMAOLD)[0].m()*(*V)[1].x() + 2.0 * (*SIGMAOLD)[1].m()*(*V)[1].y())) * N.m();
    // b[1] -= lambda * ((*SIGMAOLD)[0].m()*(*V)[2].x() + (*SIGMAOLD)[1].m()*(*V)[1].x()
    //            + (*SIGMAOLD)[1].m()*(*V)[2].y() + (*SIGMAOLD)[2].m()*(*V)[1].y()) * N.m();
    // b[2] -= lambda * (2.0 * ((*SIGMAOLD)[1].m()*(*V)[2].x() + 2.0 * (*SIGMAOLD)[2].m()*(*V)[2].y())) * N.m();

    // - lambda (grad(v)B+Bgrad(v)^T) componentwise
    b[0] -= lambda * (2.0 * ((U)[0].m()*(*V)[1].x() + 2.0 * (U)[1].m()*(*V)[1].y())) * N.m();
    b[1] -= lambda * ((U)[0].m()*(*V)[2].x() + (U)[1].m()*(*V)[1].x()
                + (U)[1].m()*(*V)[2].y() + (U)[2].m()*(*V)[1].y()) * N.m();
    b[2] -= lambda * (2.0 * ((U)[1].m()*(*V)[2].x() + 2.0 * (U)[2].m()*(*V)[2].y())) * N.m();


    // ADD STABILIZATION
    double test = (*V)[1].m()*N.x() + (*V)[2].m() * N.y();
    for (int i=0;i<3;++i)
      for (int j=0;j<2;++j)
  {
    b[i] += lambda * __DT / 3.0 * (*V)[j+1].m()*(*SIGMAOLD)[i][j+1] * test;
    b[i] += lambda * __DT / 6.0 * (*V)[j+1].m()*U[i][j+1] * test;
  }
  }


  void StressEQ::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    for (int i=0;i<3;++i)
      {
  A(i,i) += lambda * M.m()/__DT * N.m();

  // for (int j=0;j<2;++j)
  //   //    b[i] += lambda * (*V)[j+1].m()*(*SIGMAOLD)[i][j+1] * N.m();
  //   A(i,i) += lambda * (*V)[j+1].m()*M[j+1] * N.m();

  A(i,i) += M.m() * N.m();
      }

    A(0,0) -= lambda * (2.0 * M.m()*(*V)[1].x() ) * N.m();
    A(0,1) -= lambda * (2.0 * M.m()*(*V)[1].y() ) * N.m();
    A(1,0) -= lambda * (M.m()*(*V)[2].x()) * N.m();

    A(1,1) -= lambda * (M.m()*(*V)[1].x()  + M.m()*(*V)[2].y()) * N.m();
    A(1,2) -= lambda * (M.m()*(*V)[1].y()) * N.m();

    A(2,1) -= lambda * (2.0 * (M.m()*(*V)[2].x())) * N.m();
    A(2,2) -= lambda * (2.0 * M.m()*(*V)[2].y()) * N.m();


    // ADD STABILIZATION
    double test = (*V)[1].m()*N.x() + (*V)[2].m() * N.y();
    for (int i=0;i<3;++i)
      for (int j=0;j<2;++j)
  {
    //b[i] += lambda * __DT / 3.0 * (*V)[j+1].m()*(*SIGMAOLD)[i][j+1] * test;
    A(i,i) += lambda * __DT / 6.0 * (*V)[j+1].m()*M[j+1] * test;
  }
  }




  void StressEQ::point_M(int j, const FemFunction& U, const TestFunction& M) const
  {

  }






  void StressEQ::MatrixBlock(EntryMatrix& A, const FemFunction& U, const FemFunction& N) const
  {
    for (int j=0; j<N.size();++j)  // trial
      {
#define M N[j]
  point_M(j,U,M);


  for (int i=0; i<N.size();++i) // test
    {
      A.SetDofIndex(i,j);
      Matrix(A,U,M , N[i]);
    }
#undef M
      }
  }



  void StressEQ::lpspoint(double h, const FemFunction& U, const Vertex<2>& v) const
  {
    lpsstress = lpsstress0 * h*h;
    lpsstress = 0.0;

  }

  void StressEQ::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
  {
    for (int i=0;i<3;++i)
      //      for (int j=0;j<2;++j)
      b[i] += lpsstress * UP[i].m() * N.m();
    //    als nächstes stabilisieren von nabla v in der anderen gleichung
  }


  void StressEQ::StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
  {
    for (int i=0;i<3;++i)
      //      for (int j=0;j<2;++j)
      A(i,i) += lpsstress *Mp.m() * Np.m();
  }



}

/*-----------------------------------------*/
