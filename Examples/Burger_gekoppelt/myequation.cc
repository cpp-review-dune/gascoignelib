#include "myequation.h"
#include "filescanner.h"

extern double DT;
extern double DTM1;
extern double DTM2;
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
  b[0] += (U[0].m() - (*oldu)[0].m())*(1.0+(*H)[0].m()) * N.m();
  b[1] += (U[1].m() - (*oldu)[1].m())*(1.0+(*H)[0].m()) * N.m();
  
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
  A(0,0) += M.m() *(1.0+(*H)[0].m())* N.m() ;
  A(1,1) += M.m() *(1.0+(*H)[0].m())* N.m() ;

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

  void MyDualEquation::Nonlinear(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w,int DTM) const
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
  b[0] += (Z[0].m()-(*oldz)[0].m()) *(1.0)*N.m();
  b[1] += (Z[1].m()-(*oldz)[1].m()) *(1.0)* N.m();
  
if(!LASTDUAL){
 b[0] += (Z[0].m()-(*oldz)[0].m())*(*h)[0].m() *N.m();
 b[1] += (Z[1].m()-(*oldz)[1].m())*(*h)[1].m() *N.m();
}

if(LASTDUAL){
 b[0] += ((*oldz)[0].m())*(*h)[0].m() *N.m();
 b[1] += ((*oldz)[1].m())*(*h)[1].m() *N.m();
}

  // Laplace.
  if (!LASTDUAL)
    {
      b[0] +=epsilon * DTM1/2*( Z[0].x() * N.x() + Z[0].y() * N.y());
      b[1] +=epsilon * DTM1/2*( Z[1].x() * N.x() + Z[1].y() * N.y());
    }
  if (!FIRSTDUAL)
    {
     b[0] +=epsilon * DTM2/2*( (*oldz)[0].x() * N.x() + (*oldz)[0].y() * N.y());
     b[1] +=epsilon * DTM2/2*( (*oldz)[1].x() * N.x() + (*oldz)[1].y() * N.y());
    }

  // Nichtlinearitaet. u1 zu t_m-1, u2 zu t_m und u3 zu t_m+1
  if (!LASTDUAL)
    {
     Nonlinear(b, 0.5+0.5/sqrt(3.0), (*u1), (*u2), Z, N,        0.5-0.5/sqrt(3.0),DTM1);
     Nonlinear(b, 0.5-0.5/sqrt(3.0), (*u1), (*u2), Z, N,        0.5+0.5/sqrt(3.0),DTM1);
    }
  if (!FIRSTDUAL)
    {
      Nonlinear(b, 0.5+0.5/sqrt(3.0), (*u2), (*u3), (*oldz), N,  0.5+0.5/sqrt(3.0),DTM2);
     Nonlinear(b, 0.5-0.5/sqrt(3.0), (*u2), (*u3), (*oldz), N,  0.5-0.5/sqrt(3.0),DTM2);
    }
   //div(phi) hw
   if (!LASTDUAL)
   {
    b[0]+=DTM1/2.*((*h)[0].m()*(*w)[0].m())*(N.x());   
    b[1]+=DTM1/2.*((*h)[0].m()*(*w)[0].m())*(N.y());  
   }
   //v grad hw
   if (!LASTDUAL)
   {
    b[0]+=DTM1/2.*((*h)[0].x()*(*w)[0].m())*(N.m());   
    b[1]+=DTM1/2.*((*h)[0].y()*(*w)[0].m())*(N.m());  
   }
    //div(phi) hw
   if(!FIRSTDUAL)
   {
    b[0]+=DTM2/2.*((*newH)[0].m()*(*oldW)[0].m())*(N.x());   
    b[1]+=DTM2/2.*((*newH)[0].m()*(*oldW)[0].m())*(N.y());
   }
   
    //v grad hw
   if (!FIRSTDUAL)
   {
    b[0]+=DTM2/2.*((*newH)[0].x()*(*oldW)[0].m())*(N.m());   
    b[1]+=DTM2/2.*((*newH)[0].y()*(*oldW)[0].m())*(N.m());  
   }
   
}

/*----------------------------------------------------------------------------*/

void MyDualEquation::Nonlinear_Matrix(EntryMatrix&A, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z,const TestFunction &M, const TestFunction& N,double w,int DTM) const
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
  A(0,0) += M.m()*(1.0)*N.m();
  A(1,1) += M.m()*(1.0)*N.m();
  
  
  if(!LASTDUAL){
 A(0,0)+= (M.m())*(*h)[0].m() *N.m();
 A(1,1)+= (M.m())*(*h)[1].m() *N.m();
}
// Laplace.
  if (!LASTDUAL)
    {
      A(0,0) +=epsilon * DTM1/2*( M.x() * N.x() + M.y() * N.y());
      A(1,1) +=epsilon * DTM1/2*( M.x() * N.x() + M.y() * N.y());
    }
 
 // Nichtlinearitaet. u1 zu t_m-1, u2 zu t_m und u3 zu t_m+1

 
  if (!LASTDUAL)
    {
     Nonlinear_Matrix(A, 0.5+0.5/sqrt(3.0), (*u1), (*u2), Z, M,N, 0.5-0.5/sqrt(3.0),DTM1);
     Nonlinear_Matrix(A, 0.5-0.5/sqrt(3.0), (*u1), (*u2), Z, M,N, 0.5+0.5/sqrt(3.0),DTM1);
    }
 
   
}

/*----------------------------------------------------------------------------*/

// Transport

  MyTransportEquation::MyTransportEquation(const ParamFile* pf) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("epsilon",   &epsilon, 0.0);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(pf, "Equation");
  assert(epsilon>0.0);
}

/*----------------------------------------------------------------------------*/

void MyTransportEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
{
 
}

/*----------------------------------------------------------------------------*/



void MyTransportEquation::Form(VectorIterator b,
                      const FemFunction &U,
                      const TestFunction &N) const
{
  //Zeit
  b[0] += (U[0].m() -(*oldh)[0].m()) * N.m();
  b[1] += (U[1].m() -(*oldh)[1].m()) * N.m();  
  // ganz einfache stabilisierung...
  b[0] += 0.1*DT * (U[0].x()*N.x() + U[0].y()*N.y());
  b[1] += 0.1*DT * (U[1].x()*N.x() + U[1].y()*N.y());


  // //// div v*h
  b[0]+=DT/2.* ((*V)[0].x()+(*V)[1].y())*U[0].m()*N.m();
  b[0]+=DT/2.* ((*V)[0].x()+(*V)[1].y())*(*oldh)[0].m()*N.m();
  
  //   //// div v*h
  b[1]+=DT/2.* ((*V)[0].x()+(*V)[1].y())*U[1].m()*N.m();
  b[1]+=DT/2.* ((*V)[0].x()+(*V)[1].y())*(*oldh)[1].m()*N.m();
  // v nabla h
  
  b[0]+=DT/2.*((*V)[0].m()*U[0].x()      + (*V)[1].m()*U[0].y()      )*N.m();
  b[0]+=DT/2.*((*V)[0].m()*(*oldh)[0].x()+ (*V)[1].m()*(*oldh)[0].y())*N.m();

  //2te Komponente


  // v nabla h
  b[1]+=DT/2.*((*V)[0].m()*U[1].x()+(*V)[1].m()*U[1].y())*N.m();
  b[1]+=DT/2.*((*V)[0].m()*(*oldh)[1].x()+(*V)[1].m()*(*oldh)[1].y())*N.m(); 

}

/*----------------------------------------------------------------------------*/

void MyTransportEquation::Matrix(EntryMatrix &A,
                        const FemFunction &U,
                        const TestFunction &M,
                        const TestFunction &N) const
{
  //Zeit
  A(0,0) += M.m() * N.m() ;
  A(1,1) += M.m() * N.m() ;


  A(0,0) += 0.1*DT * (M.x()*N.x() + M.y()*N.y());
  A(1,1) += 0.1*DT * (M.x()*N.x() + M.y()*N.y());

  

  A(0,0)+=DT/2.*((*V)[0].x()+(*V)[1].y())*M.m()*N.m();
  A(1,1)+=DT/2.*((*V)[0].x()+(*V)[1].y())*M.m()*N.m();


  A(0,0)+=DT/2.*((*V)[0].m()*M.x()+(*V)[1].m()*M.y())*N.m();
  A(1,1)+=DT/2.*((*V)[0].m()*M.x()+(*V)[1].m()*M.y())*N.m();
  
  
}


  
// Transport_DUaL

  MyDualTransportEquation::MyDualTransportEquation(const ParamFile* pf) : Equation()
{
  DataFormatHandler DFH;
  DFH.insert("epsilon",   &epsilon, 0.0);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(pf, "Equation");
  assert(epsilon>0.0);
}

/*----------------------------------------------------------------------------*/

void MyDualTransportEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
{
 
}

/*----------------------------------------------------------------------------*/

 void MyDualTransportEquation::Nonlinear(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w,int DTM) const
  {
    //div v h 
    b[0] += w*DTM * 0.5 * ( (s*(U1[0].x()+U1[1].y())+(1-s)*(U2[0].x()+U2[1].y())) * N.m()*Z[0].m();
    b[1] += w*DTM * 0.5 * ( (s*(U1[0].x()+U1[1].y())+(1-s)*(U2[0].x()+U2[1].y())) * N.m()*Z[1].m();
			    //v nabla h
    b[0] += w*DTM * 0.5*((s*U1[0].m()+(1-s)*U2[0].m())*N.x()+ (s*U1[1].m()+(1-s)*U2[1].m())*N.y())*Z[0].m());
    b[1] += w*DTM * 0.5*((s*U1[0].m()+(1-s)*U2[0].m())*N.x()+ (s*U1[1].m()+(1-s)*U2[1].m())*N.y())*Z[1].m());
  }

void MyDualTransportEquation::Kopplung(VectorIterator b, double s, const FemFunction &U1, const FemFunction& U2, const FemFunction& Z, const TestFunction& N,double w,int DTM) const
  {
    //div v h 
    b[0] += w* 0.5 * ( (s*(U1[0].m()-U2[0].m())+(1-s)*(U1[0].m()-U0[0].y())) * N.m()*OLDZ[0].m();
    b[1] += w*DTM * 0.5 * ( (s*(U1[0].x()+U1[1].y())+(1-s)*(U2[0].x()+U2[1].y())) * N.m()*Z[1].m();
			    //v nabla h
    b[0] += w*DTM * 0.5*((s*U1[0].m()+(1-s)*U2[0].m())*N.x()+ (s*U1[1].m()+(1-s)*U2[1].m())*N.y())*Z[0].m());
    b[1] += w*DTM * 0.5*((s*U1[0].m()+(1-s)*U2[0].m())*N.x()+ (s*U1[1].m()+(1-s)*U2[1].m())*N.y())*Z[1].m());
  }
  

void MyDualTransportEquation::Form(VectorIterator b,
                      const FemFunction &U,
                      const TestFunction &N) const
{
  //Zeit
  b[0] += (U[0].m() -(*oldw)[0].m()) * N.m();
  b[1] += (U[1].m() -(*oldw)[1].m()) * N.m();


   // ganz einfache stabilisierung...
  b[0] += 0.01*DT * (U[0].x()*N.x() + U[0].y()*N.y());
  b[1] += 0.01*DT * (U[1].x()*N.x() + U[1].y()*N.y());

  // //// div v*h  +v nabla h
  
  // Nichtlinearitaet. u1 zu t_m-1, u2 zu t_m und u3 zu t_m+1
  if (!LASTDUAL){

   Nonlinear(b, 0.5+0.5/sqrt(3.0), (*u1), (*u2), Z, N,        0.5-0.5/sqrt(3.0),DTM1);
   Nonlinear(b, 0.5-0.5/sqrt(3.0), (*u1), (*u2), Z, N,        0.5+0.5/sqrt(3.0),DTM1);

  }
  

 if (!LASTDUAL){
   // koppolungsterm
  b[0]+=DTM1/2*(((*V)[0].m()-(*newV)[0].m())/DTM1*(*z)[0].m()+((*V)[1].m()-(*newV)[1].m())/DTM1*(*z)[1].m())*N.m();
  b[1]+=DTM1/2*(((*V)[0].m()-(*newV)[0].m())/DTM1*(*z)[0].m()+((*V)[1].m()-(*newV)[1].m())/DTM1*(*z)[1].m())*N.m();
 }

  if (!FIRSTDUAL){
  //   //// div v*h
  b[0]+=DTM2/2.* ((*newV)[0].x()+(*newV)[1].y())*(*oldw)[0].m()*N.m();
  b[1]+=DTM2/2.* ((*newV)[0].x()+(*newV)[1].y())*(*oldw)[1].m()*N.m();
  
  // v nabla h
  
  b[0]+=DTM2/2.*((*newV)[0].m()*N.x()+ (*newV)[1].m()*N.y())*(*oldw)[0].m();
  b[1]+=DTM2/2.*((*newV)[0].m()*N.x()+(*newV)[1].m()*N.y())*(*oldw)[1].m();


  //koppolung
 b[0]+=DTM2/2*(((*newV)[0].m()-(*newnewV)[0].m())/DTM2*(*oldz)[0].m()+((*newV)[1].m()-(*newnewV)[1].m())/DTM2*(*oldz)[1].m())*N.m();
 b[1]+=DTM2/2*(((*newV)[0].m()-(*newnewV)[0].m())/DTM2*(*oldz)[0].m()+((*newV)[1].m()-(*newnewV)[1].m())/DTM2*(*oldz)[1].m())*N.m();
  }
}

/*----------------------------------------------------------------------------*/

void MyDualTransportEquation::Matrix(EntryMatrix &A,
                        const FemFunction &U,
                        const TestFunction &M,
                        const TestFunction &N) const
{
  //Zeit
  A(0,0) += M.m() * N.m() ;
  A(1,1) += M.m() * N.m() ;

  A(0,0) += 0.01*DT * (M.x()*N.x() + M.y()*N.y());
  A(1,1) += 0.01*DT * (M.x()*N.x() + M.y()*N.y());
  
 if (!LASTDUAL){
  //b[0]+=DTM1/2.* ((*V)[0].x()+(*V)[1].y())*U[0].m()*N.m();
  
  A(0,0)+=DTM1/2.*((*V)[0].x()+(*V)[1].y())*M.m()*N.m();
  A(1,1)+=DTM1/2.*((*V)[0].x()+(*V)[1].y())*M.m()*N.m();

  //b[0]+=DTM1/2.*((*V)[0].m()*N.x()      + (*V)[1].m()*N.y()      )*U[0].m();
  A(0,0)+=DTM1/2.*((*V)[0].m()*N.x()+(*V)[1].m()*N.y())*M.m();
  A(1,1)+=DTM1/2.*((*V)[0].m()*N.x()+(*V)[1].m()*N.y())*M.m();

 }
  
}



  

 
}
