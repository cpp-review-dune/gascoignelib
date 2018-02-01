#include "myequation.h"
#include "filescanner.h"

extern double DT;
extern double DTM;
extern double CHIGAUSS;
extern bool FIRSTDUAL;
extern bool LASTDUAL;

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
 
}

/*----------------------------------------------------------------------------*/

void MyEquation::Form(VectorIterator b,
                      const FemFunction &U,
                      const TestFunction &N) const
{
  //Zeit // unabhaengig vom Gausspunkt, da (d_t Uk) konstant.
  b[0] += (U[0].m() - (*oldu)[0].m()) * N.m();
  b[1] += (U[1].m() - (*oldu)[1].m()) * N.m();
  
  // Lineare und nichtlineare Anteile...
  // U(chigauss) ausrechnen, lineare Interpolation
  FemFunction Ug(U.size());
  for (int c=0;c<U.size();++c)
    Ug[c].equ(1.0-CHIGAUSS, (*oldu)[c], CHIGAUSS , U[c]);
  
  //// Laplace
  b[0] +=DT * epsilon*( Ug[0].x() * N.x() + Ug[0].y() * N.y());
  b[1] +=DT * epsilon*( Ug[1].x() * N.x() + Ug[1].y() * N.y());

  //Nichtlineraritaet

  b[0]+=DT*(Ug[0].m()*Ug[0].x()+Ug[1].m()*Ug[0].y())*N.m();
  b[1]+=DT*(Ug[0].m()*Ug[1].x()+Ug[1].m()*Ug[1].y())*N.m();

    
}

/*----------------------------------------------------------------------------*/

void MyEquation::Matrix(EntryMatrix &A,
                        const FemFunction &U,
                        const TestFunction &M,
                        const TestFunction &N) const
{
  //Zeit, unabhaengig vom Gausspunkt
  A(0,0) += M.m() * N.m() ;
  A(1,1) += M.m() * N.m() ;

  // Lineare und nichtlineare Anteile...
  // U(chigauss) ausrechnen, lineare Interpolation
  FemFunction Ug(U.size());
  for (int c=0;c<U.size();++c)
    Ug[c].equ(1.0-CHIGAUSS, (*oldu)[c], CHIGAUSS , U[c]);

  //Laplace
  A(0, 0) +=DT *  CHIGAUSS *epsilon* (M.x() * N.x() + M.y() * N.y());
  A(1, 1) +=DT * CHIGAUSS * epsilon* (M.x() * N.x() + M.y() * N.y());

  //  b[0]+=DT*(U[0].m()*U[0].x()+U[1].m()*U[0].y())*N.m();
  // b[1]+=DT*(U[0].m()*U[1].x()+U[1].m()*U[1].y())*N.m();
  
  //Nichtlinear
  
  A(0,0)+=DT* CHIGAUSS *(M.m()*Ug[0].x()+Ug[0].m()*M.x()+Ug[1].m()*M.y())*N.m();
  A(0,1)+=DT* CHIGAUSS *(M.m()*Ug[0].y())*N.m();
  
 
 
  A(1,1)+=DT*CHIGAUSS *(Ug[0].m()*M.x()+M.m()*Ug[1].y()+Ug[1].m()*M.y())*N.m();
  A(1,0)+=DT*CHIGAUSS*(M.m()*Ug[1].x())*N.m();
  
  

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
  

}

/*----------------------------------------------------------------------------*/

  void MyDualEquation::Nonlinear(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w) const
  {
    b[0] += w*DTM * 0.5 * ( (s*U1[0].m()+(1-s)*U2[0].m()) * N.x() + (s*U1[1].m()+(1-s)*U2[1].m()) * N.y()) * Z[0].m();
    b[1] += w*DTM * 0.5 * ( (s*U1[0].m()+(1-s)*U2[0].m()) * N.x() + (s*U1[1].m()+(1-s)*U2[1].m()) * N.y()) * Z[1].m();

    b[0] += w*DTM * 0.5 * N.m() * (  (s*U1[0].x()+(1-s)*U2[0].x()) * Z[0].m() + (s*U1[1].x()+(1-s)*U2[1].x()) * Z[1].m());
    b[1] += w*DTM * 0.5 * N.m() * (  (s*U1[0].y()+(1-s)*U2[0].y()) * Z[0].m() + (s*U1[1].y()+(1-s)*U2[1].y()) * Z[1].m());
  }
  
void MyDualEquation::Form(VectorIterator b,
			  const FemFunction &Z,
			  const TestFunction &N) const
{  
  //Zeit
  b[0] += (Z[0].m()-(*oldz)[0].m()) * N.m();
  b[1] += (Z[1].m()-(*oldz)[1].m()) * N.m();
   
  // Laplace.
  if (!LASTDUAL)
    {
      b[0] +=epsilon * DTM/2*( Z[0].x() * N.x() + Z[0].y() * N.y());
      b[1] +=epsilon * DTM/2*( Z[1].x() * N.x() + Z[1].y() * N.y());
    }
  if (!FIRSTDUAL)
    {
          b[0] +=epsilon * DTM/2*( (*oldz)[0].x() * N.x() + (*oldz)[0].y() * N.y());
         b[1] +=epsilon * DTM/2*( (*oldz)[1].x() * N.x() + (*oldz)[1].y() * N.y());
    }

  // Nichtlinearitaet. u1 zu t_m-1, u2 zu t_m und u3 zu t_m+1
  if (!LASTDUAL)
    {
        Nonlinear(b, 0.5+0.5/sqrt(3.0), (*u1), (*u2), Z, N,        0.5-0.5/sqrt(3.0));
         Nonlinear(b, 0.5-0.5/sqrt(3.0), (*u1), (*u2), Z, N,        0.5+0.5/sqrt(3.0));
    }
  if (!FIRSTDUAL)
    {
         Nonlinear(b, 0.5+0.5/sqrt(3.0), (*u2), (*u3), (*oldz), N,  0.5+0.5/sqrt(3.0));
         Nonlinear(b, 0.5-0.5/sqrt(3.0), (*u2), (*u3), (*oldz), N,  0.5-0.5/sqrt(3.0));
    }
}

/*----------------------------------------------------------------------------*/

void MyDualEquation::Nonlinear_Matrix(EntryMatrix&A, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z,const TestFunction &M, const TestFunction& N,double w) const
  {
    //  b[0] += w*DTM * 0.5 * ( (s*U1[0].m()+(1-s)*U2[0].m()) * N.x() + (s*U1[1].m()+(1-s)*U2[1].m()) * N.y()) * Z[0].m();

    A(0,0) += w*DTM * 0.5 * ( (s*U1[0].m()+(1-s)*U2[0].m()) * N.x() + (s*U1[1].m()+(1-s)*U2[1].m()) * N.y()) * M.m();
    
    //b[1] += w*DTM * 0.5 * ( (s*U1[0].m()+(1-s)*U2[0].m()) * N.x() + (s*U1[1].m()+(1-s)*U2[1].m()) * N.y()) * Z[1].m();
    A(1,1)+= w*DTM * 0.5 * ( (s*U1[0].m()+(1-s)*U2[0].m()) * N.x() + (s*U1[1].m()+(1-s)*U2[1].m()) * N.y()) * M.m();

    // b[0] += w*DTM * 0.5 * N.m() * (  (s*U1[0].x()+(1-s)*U2[0].x()) * Z[0].m() + (s*U1[1].x()+(1-s)*U2[1].x()) * Z[1].m());

    A(0,0)+= w*DTM * 0.5 * N.m() *  (s*U1[0].x()+(1-s)*U2[0].x()) * M.m();
    A(0,1) +=w*DTM * 0.5 * N.m() *  (s*U1[1].x()+(1-s)*U2[1].x()) * M.m();
    
				     
    // b[1] += w*DTM * 0.5 * N.m() * (  (s*U1[0].y()+(1-s)*U2[0].y()) * Z[0].m() + (s*U1[1].y()+(1-s)*U2[1].y()) * Z[1].m());

    A(1,1) += w*DTM * 0.5 * N.m() *(s*U1[1].y()+(1-s)*U2[1].y()) * M.m();
    A(1,0)+= w*DTM * 0.5 * N.m() * (s*U1[0].y()+(1-s)*U2[0].y()) * M.m(); 


    
  }

void MyDualEquation::Matrix(EntryMatrix &A,
                        const FemFunction &Z,
                        const TestFunction &M,
                        const TestFunction &N) const
{

  //Zeit
  A(0,0) += M.m()*N.m();
  A(1,1) += M.m()*N.m();
  
  // Laplace.
  if (!LASTDUAL)
    {
      A(0,0) +=epsilon * DTM/2*( M.x() * N.x() + M.y() * N.y());
      A(1,1) +=epsilon * DTM/2*( M.x() * N.x() + M.y() * N.y());
    }
 
 // Nichtlinearitaet. u1 zu t_m-1, u2 zu t_m und u3 zu t_m+1

 
  if (!LASTDUAL)
    {
        Nonlinear_Matrix(A, 0.5+0.5/sqrt(3.0), (*u1), (*u2), Z, M,N,        0.5-0.5/sqrt(3.0));
        Nonlinear_Matrix(A, 0.5-0.5/sqrt(3.0), (*u1), (*u2), Z, M,N,        0.5+0.5/sqrt(3.0));
    }
 
   
}

/*----------------------------------------------------------------------------*/
  

 
}
