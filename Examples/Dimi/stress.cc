#include  "stress.h"
#include  "filescanner.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  Stress::Stress(const ParamFile* paramfile) : Equation()
  {
    DataFormatHandler DFH;
    // DFH.insert("tau", &tau0, 0.);
    // DFH.insert("shock", &shock0, 0.);
    DFH.insert("rho", &rho, 0.);
    // DFH.insert("rhow", &rhow, 0.);
    DFH.insert("Tref",&Tref,0.0);
    // DFH.insert("Lref",&Lref,0.0);
     DFH.insert("Pstern",&Pstern,2.75e4);
    DFH.insert("ellipse",&ellipse,2.0);
    DFH.insert("C",&C,20.0);
    // DFH.insert("Cdw",&Cdw,5.2e-3);
    //  DFH.insert("f",&f,0.0);
    DFH.insert("deltamin",&deltamin,2.e-9);
    // DFH.insert("theta_w",&theta_w,0.0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(paramfile, "Equation");
  }
  
  
  void Stress::point(double h, const FemFunction &U, const Vertex2d &v) const
  {
  }
 
  /*-----------------------------------------*/

  void Stress::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
   
    
    b[0] += U[0].m()*N.m();
    b[1] += U[1].m()*N.m();
    b[2] += U[2].m()*N.m();
    b[3] += U[3].m()*N.m();

     // Exponentialfaktor
    double ef = exp(-C*(1.0-(*V)[3].m()));
double P= ef*(*V)[2].m()*Pstern;
    double dmin = deltamin* Tref;
    double DELTAsquare = //  dmin*dmin+
      (1.0 + pow(ellipse,-2.0)) * ((*V)[0].x()*(*V)[0].x() + (*V)[1].y()*(*V)[1].y())
      + pow(ellipse,-2.0) * pow((*V)[0].y() + (*V)[1].x(),2.0)
      + 2.0 * (1.0-pow(ellipse,-2.0)) * (*V)[0].x()*(*V)[1].y();
    DELTAsquare = std::max(DELTAsquare,dmin*dmin);
    
    double DELTA = sqrt(DELTAsquare);
    
    double sigma_0= (*V)[2].m() * pow(ellipse,-2.0) * ef / DELTA * ( (*V)[0][1] + (*V)[0][1] ) 
    +  (*V)[2].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * (*V)[0][1]
    +- (*V)[2].m() * ef ;
    double sigma_1=(*V)[2].m() * pow(ellipse,-2.0) * ef / DELTA * ( (*V)[1][2] + (*V)[1][2] ) 
    +(*V)[2].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * (*V)[1][2]
    - (*V)[2].m() * ef;
	
    b[0]-= ( ((sigma_0+sigma_1+P)/P)*((sigma_0+sigma_1+P)/P)+(ellipse*(sigma_0-sigma_1)/P)*(ellipse*(sigma_0-sigma_1)/P)-1.0 ) * N.m();
    
    
    ///////  (sigma, nabla phi)	
    //hier. sigma phi
    // nabla u, nabla phi
    // nable u^T, nabla phi
    //for (int i=0;i<2;++i)
      //for (int j=0;j<2;++j)
	//{
	//  sigma[i] -=  U[2].m() * pow(ellipse,-2.0) * ef / DELTA * ( U[i][j+1] + U[j][i+1] ) * N.m();
	  //sigma[i] -=  U[2].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * U[j][j+1]*N.m();
	//}
    
    // // ef, siehe oben, exponentialfaktor
  //  sigma[0] -= - U[2].m() * ef * N.m();
    //sigma[1] -= - U[2].m() * ef * N.m();
  }
  /*-----------------------------------------*/
  
  void Stress::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    A(0,0)  += M.m()*N.m();
    A(1,1)  += M.m()*N.m();
    A(2,2)  += M.m()*N.m();
    A(3,3)  += M.m()*N.m();
  }
}

/*-----------------------------------------*/
