#include "myequation.h"
#include "filescanner.h"

extern double DT;
extern double DTM;
extern double CHIGAUSS;

namespace Gascoigne
{
/*----------------------------------------------------------------------------*/

MyEquation::MyEquation(const ParamFile* pf) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("epsilon",   &epsilon, 0.0);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(pf, "Equation");
  assert(epsilon>0.0);
}

/*----------------------------------------------------------------------------*/

void MyEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
{
  beta_x = (v.y()-0.5)*2.0*M_PI;
  beta_y =-(v.x()-0.5)*2.0*M_PI;
}

/*----------------------------------------------------------------------------*/

void MyEquation::Form(VectorIterator b,
                      const FemFunction &U,
                      const TestFunction &N) const
{
  //Zeit // unabhaengig vom Gausspunkt, da (d_t Uk) konstant.
  b[0] += (U[0].m() - (*oldu)[0].m()) * N.m();
  
  // Lineare und nichtlineare Anteile...
  // U(chigauss) ausrechnen, lineare Interpolation
  FemFunction Ug(U.size());
  for (int c=0;c<U.size();++c)
    Ug[c].equ(1.0-CHIGAUSS, (*oldu)[c], CHIGAUSS , U[c]);
  
  //// Laplace
  b[0] +=DT * epsilon*( Ug[0].x() * N.x() + Ug[0].y() * N.y());
  b[0] +=DT * ( Ug[0].x() * Ug[0].x() + Ug[0].y() * Ug[0].y())*N.m();

  //Transport
  b[0]+= DT * (beta_x * Ug[0].x() + beta_y * Ug[0].y()) * N.m();
}

/*----------------------------------------------------------------------------*/

void MyEquation::Matrix(EntryMatrix &A,
                        const FemFunction &U,
                        const TestFunction &M,
                        const TestFunction &N) const
{
  //Zeit, unabhaengig vom Gausspunkt
  A(0,0) += M.m() * N.m() ;

  // Lineare und nichtlineare Anteile...
  // U(chigauss) ausrechnen, lineare Interpolation
  FemFunction Ug(U.size());
  for (int c=0;c<U.size();++c)
    Ug[c].equ(1.0-CHIGAUSS, (*oldu)[c], CHIGAUSS , U[c]);

  //Laplace
  A(0, 0) +=DT * CHIGAUSS *  epsilon* (M.x() * N.x() + M.y() * N.y());
  
  A(0,0) += DT* CHIGAUSS * (2.*Ug[0].x()*M.x()+2.*Ug[0].y()*M.y())*N.m();
 
  //Transport
   A(0,0)+= DT * CHIGAUSS * (beta_x * M.x() + beta_y * M.y()) * N.m();
}


//dual problem
/*----------------------------------------------------------------------------*/
MyDualEquation::MyDualEquation(const ParamFile* pf) : Equation()
{
    DataFormatHandler DFH;
  DFH.insert("epsilon",   &epsilon, 0.);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(pf, "Equation");
  assert(epsilon>0.0);
}

/*----------------------------------------------------------------------------*/

void MyDualEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
{
   beta_x = (v.y()-0.5)*2.0*M_PI;
  beta_y =-(v.x()-0.5)*2.0*M_PI;

}

/*----------------------------------------------------------------------------*/

void MyDualEquation::Form(VectorIterator b,
                      const FemFunction &U,
			  const TestFunction &N) const
{

  
  
  
  // Laplace
  b[0] +=epsilon *  DTM*( U[0].x() * N.x() + U[0].y() * N.y());
  b[0] += (U[0].m()-(*oldz)[0].m()) * N.m();
  
  //Transport
  b[0]+=DTM*(beta_x * U[0].m()*N.x() + beta_y * U[0].m()* N.y());
  b[0] += DTM * 2.0* ( (*Pu)[0].x() * N.x() + (*Pu)[0].y() * N.y())*U[0].m();
  // b[0] += DTM * N.m() * N.m() *U[0].m();
  
}

/*----------------------------------------------------------------------------*/

void MyDualEquation::Matrix(EntryMatrix &A,
                        const FemFunction &U,
                        const TestFunction &M,
                        const TestFunction &N) const
{


  //Laplace
  A(0,0) +=epsilon *  DTM*(M.x() * N.x() + M.y() * N.y());
  A(0,0) += M.m()*N.m();

  //Transport
   A(0,0)+= DTM*(beta_x * M.m()*N.x() +beta_y * M.m() * N.y());
   A(0,0) += DTM * 2.0* ( (*Pu)[0].x() * N.x() + (*Pu)[0].y() * N.y())*M.m();
   // A(0,0) += DTM* (N.x()*N.x()+N.y()*N.y())*M.m();
   // A(0,0) += DTM * N.m() * N.m() *M.m();
}

/*----------------------------------------------------------------------------*/
  MyDualHelpEquation::MyDualHelpEquation(const ParamFile* pf) : Equation()
{
    DataFormatHandler DFH;
  DFH.insert("epsilon",   &epsilon, 0.);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(pf, "Equation");
  assert(epsilon>0.0);

}

/*----------------------------------------------------------------------------*/

void MyDualHelpEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
{
   beta_x = (v.y()-0.5)*2.0*M_PI;
  beta_y =-(v.x()-0.5)*2.0*M_PI;

}




void MyDualHelpEquation::Form(VectorIterator b,
                      const FemFunction &U,
                      const TestFunction &N) const
{
  b[0] += epsilon * DTM * (U[0].x() * N.x() + U[0].y() * N.y());
  b[0]+= DTM * (beta_x * U[0].x() + beta_y * U[0].y()) * N.m();
  
  b[0] += ((*U_h)[0].m()-(*Uold)[0].m())*N.m();
  
   b[0] += DTM * ( U[0].x() * U[0].x() + U[0].y() * U[0].y())*N.m();
  // b[0] += DTM * U[0].m() * U[0].m() *N.m();

}

/*----------------------------------------------------------------------------*/

void MyDualHelpEquation::Matrix(EntryMatrix &A,
                        const FemFunction &U,
                        const TestFunction &M,
                        const TestFunction &N) const
{

  abort();
}


 
}
