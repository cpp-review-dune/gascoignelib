#include  "myequation.h"
#include  "filescanner.h"


/*-----------------------------------------*/
extern double __GLOBAL_TIME;



namespace Gascoigne

{

  double p =2.0;
  
  
  
  
  
  

  
  /*-----------------------------------------*/

  MyEquation::MyEquation(const ParamFile* pf) : Equation()
  {
    DataFormatHandler DFH;
    DFH.insert("eps",&__eps,1.);
    FileScanner FS(DFH,pf,"Equation");
  }

  /*-----------------------------------------*/

  void MyEquation::point(double h, const FemFunction& U, const Vertex2d& v) const
  {
    __b0 =  v.y();
    __b1 = -v.x();

    // if (__GLOBAL_TIME>M_PI)
    //   {
    // 	__b1 = 1.0 / M_PI;
    // 	__b0 = 1.0/M_PI * sin(8.0 *__GLOBAL_TIME);
    //   }

    // if (__GLOBAL_TIME<M_PI/2.0)
    //   { __b0 = 1.0; __b1 = -1.0; }
    // else if (__GLOBAL_TIME<M_PI)
    //   { __b0 = -1.0; __b1 = -1.0; }
    // else if (__GLOBAL_TIME<3.0*M_PI/2.0)
    //   { __b0 = -1.0; __b1 = 1.0; }
    // else
    //   { __b0 = 1.0; __b1 = 1.0; }
    // __b0 /= M_PI; __b1 /= M_PI;
  }
 
  /*-----------------------------------------*/

  void MyEquation::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    
    double P = 1.e-5+U[0].x()*U[0].x()+U[0].y()*U[0].y();

    double X = 1.0;
    if (p!=2) X = pow(P,(p-2.0)/2.0);
    
    b[0] += __eps * X * (U[0].x() * N.x() + U[0].y() * N.y());
    b[0] += (__b0 * U[0].x() + __b1 * U[0].y()) * N.m();

    //b[0] += pow(U[0].m(),3.0) * N.m();
    b[0] += exp(U[0].m()) * N.m();
  }

  /*-----------------------------------------*/

  void MyEquation::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    double P = 1.e-5 + U[0].x()*U[0].x()+U[0].y()*U[0].y();

    double X = 1.0;
    if (p!=2) X = pow(P,(p-2.0)/2.0);

    double PD = 2.0 * U[0].x() * M.x() +  2.0 * U[0].y() * M.y();
    
    A(0,0) += __eps *  X * (M.x() * N.x() + M.y() * N.y());
    
    if (p!=2.0)
      A(0,0) += __eps *  (p-2.0)/2.0 * pow(P,(p-2.0)/2.0-1) * PD * (U[0].x() * N.x() + U[0].y() * N.y());

    A(0,0) += (__b0 * M.x() + __b1 * M.y()) * N.m();


    //A(0,0) += 3.0 * pow(U[0].m(),2.0) * M.m() * N.m();
    A(0,0) += exp(U[0].m()) * M.m() * N.m();
  }
  
}

/*-----------------------------------------*/
