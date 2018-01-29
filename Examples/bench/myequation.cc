#include  "myequation.h"
#include  "filescanner.h"


extern double TIME;
extern double DT,DDD;
extern double DELTAMIN;
extern double STEUERUNG_MU,zahl;


/*-----------------------------------------*/

namespace Gascoigne
{
  MyEquation::MyEquation(const ParamFile* paramfile) : Equation()
  {
    DataFormatHandler DFH;



    DFH.insert("v_in", &vin0, 0.);
    DFH.insert("gamma", &gamma0, 0.);

    DFH.insert("rho", &rho, 0.);
    DFH.insert("rhow", &rhow, 0.);
    DFH.insert("Tref",&Tref,0.0);
    DFH.insert("Lref",&Lref,0.0);
    DFH.insert("Pstern",&Pstern,2.75e4);
    DFH.insert("ellipse",&ellipse,2.0);
    DFH.insert("C",&C,20.0);
    DFH.insert("Cdw",&Cdw,5.2e-3);
    DFH.insert("f",&f,0.0);
    DFH.insert("theta_w",&theta_w,0.0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(paramfile, "Equation");
    assert(rho>0);
    assert(Tref>0);
    assert(Lref>0);
    assert(Pstern>0);
    assert(ellipse>0);
    assert(C>0);
    //assert(f>0);

    MZ = 0.5*Tref*Tref * Pstern / rho / Lref / Lref;
    std::cout << "Mehlmann-Zahl " << MZ << std::endl;
  }
  
  
  void MyEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
  {
    // Stabilisierung muss auch auf Dimensionslose groessen transformiert werden
    h_=h;
    
    double X = v.x()*Lref;
    double Y = v.y()*Lref;
    double Lx = 0.5e6;
    double Ly = 0.5e6;

    double Uw_x=0.1*(2*Y-Ly)/Ly;
    double Uw_y=(-0.1)*(2*X-Lx)/Lx;

    uwx= Uw_x *Tref/Lref/10;
    uwy= Uw_y *Tref/Lref/10;
  }
 
  /*-----------------------------------------*/

  // u[0], u[1] v
  
  // u[2] h    u[3] a
  
  void MyEquation::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    // Wasser-Tensor tau_w
    double WZ = rhow/rho * Cdw * Lref;
    
    double vd_x = U[0].m() - uwx;
    double vd_y = U[1].m() - uwy;
    
    double vd = sqrt(vd_x*vd_x+vd_y*vd_y+1.e-8);

      
    b[0] += WZ * vd *(vd_x * cos(theta_w) - vd_y*sin(theta_w)) * N.m();
    b[1] += WZ * vd *(vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
    
     
    //zustatz term Correolis 1.46 e -4 /s
    b[0] += -Tref*(*H)[0].m()*f*(vd_y)*N.m();
    b[1] +=  Tref*(*H)[0].m()*f*(vd_x)* N.m();
    

    // Zeitableitung
    b[0] += (*H)[0].m()*(U[0].m() - (*oldu)[0].m()) / DT * N.m();
    b[1] += (*H)[0].m()*(U[1].m() - (*oldu)[1].m()) / DT * N.m();

    // Exponentialfaktor
    double ef = exp(-C*(1.0-(*H)[1].m()));

    double dmin = DELTAMIN*Tref;
  
    // DELTA implizit
    double DELTAsquare =
      (1.0 + pow(ellipse,-2.0)) * (U[0].x()*U[0].x() + U[1].y()*U[1].y())
      + pow(ellipse,-2.0) * pow(U[0].y() + U[1].x(),2.0)
      + 2.0 * (1.0-pow(ellipse,-2.0)) * U[0].x()*U[1].y();

    DELTAsquare +=   dmin*dmin;
    
    double DELTA = sqrt(DELTAsquare);
    
    ///////  (sigma, nabla phi)	
    // nabla u, nabla phi
    // nable u^T, nabla phi
    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
	{
	  b[i] += MZ * (*H)[0].m() * pow(ellipse,-2.0) * ef / DELTA * ( U[i][j+1] + U[j][i+1] ) * N[j+1];
	  b[i] += MZ * (*H)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * U[j][j+1]*N[i+1];
	}
    
    // // ef, siehe oben, exponentialfaktor
    b[0] += -MZ * (*H)[0].m() * ef * N.x();
    b[1] += -MZ * (*H)[0].m() * ef *N.y();
  }
  /*-----------------------------------------*/
  
  void MyEquation::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    // Wasser-Tensor tau_w
    double WZ = rhow/rho * Cdw * Lref;
   
    double vd_x = U[0].m() - uwx;
    double vd_y = U[1].m() - uwy;

    double vd = sqrt(vd_x*vd_x+vd_y*vd_y+1.e-8);
    double vd_0 = 1.0/(2.0*vd) * 2.0 * vd_x * M.m();
    double vd_1 = 1.0/(2.0*vd) * 2.0 * vd_y * M.m();
 
    //b[0] += WZ * vd   *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * N.m();
    A(0,0) += WZ * vd_0 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * N.m();
    A(0,1) += WZ * vd_1 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * N.m();
    A(0,0) += WZ * vd   *(M.m() * cos(theta_w)                    ) * N.m();
    A(0,1) += WZ * vd   *(                    -M.m()*sin(theta_w) ) * N.m();
  
    //b[1] += WZ * vd   * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
    A(1,0) += WZ * vd_0 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
    A(1,1) += WZ * vd_1 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
    A(1,1) += WZ * vd   * (M.m() * cos(theta_w)                   ) * N.m();
    A(1,0) += WZ * vd   * (                     M.m()*sin(theta_w)) * N.m();
			   

    //zustatz term
    //b[0] += -Tref*(*H)[0].m()*f*(vd_y)*N.m();
    A(0,1) += -Tref*(*H)[0].m()*f*M.m() * N.m();

    //b[1] +=  Tref*(*H)[0].m()*f*(vd_x)* N.m();
    A(1,0) +=  Tref*(*H)[0].m()*f*M.m() * N.m();
 
 
    
    // Zeitableitung
    A(0,0) +=  (*H)[0].m()* M.m() / DT * N.m();
    A(1,1) +=  (*H)[0].m()* M.m() / DT * N.m();

     double dmin = DELTAMIN*Tref;
    // double dmin = DELTAMIN*DDD*Tref+(1.0-DDD)*Tref*0.05e-9;
   
    // implizit
    double DELTAsquare =
      (1.0 + pow(ellipse,-2.0)) * (U[0].x()*U[0].x() + U[1].y()*U[1].y())
      + pow(ellipse,-2.0) * pow(U[0].y() + U[1].x(),2.0)
      + 2.0 * (1.0-pow(ellipse,-2.0)) * U[0].x()*U[1].y();
      
    DELTAsquare +=  dmin*dmin;
     
    //------------------------------------------------------
    double DELTA = sqrt(DELTAsquare);
    
    // implizit
    double DELTAsquare_0 = 
      (1.0 + pow(ellipse,-2.0)) * (2.0*U[0].x()*M.x())
      + pow(ellipse,-2.0) * 2.0 * (U[0].y() + U[1].x()) * M.y()
      + 2.0 * (1.0-pow(ellipse,-2.0)) * M.x()*U[1].y();

    double DELTAsquare_1 = 
      (1.0 + pow(ellipse,-2.0)) * (2.0 * U[1].y()*M.y())
      + pow(ellipse,-2.0) * 2.0*(U[0].y() + U[1].x()) * M.x()
      + 2.0 * (1.0-pow(ellipse,-2.0)) * U[0].x()*M.y();
   
    double DELTA_0 = -0.5 * pow(DELTAsquare,-1.5) * DELTAsquare_0*DDD;
    double DELTA_1 = -0.5 * pow(DELTAsquare,-1.5) * DELTAsquare_1*DDD;

 
    // Exponentialfaktor
    double ef   = exp(-C*(1.0-(*H)[1].m()));
   
    // nabla u, nabla phi
    // nable u^T, nabla phi
    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
	{
	  //	  b[i] += MZ * (*H)[0].m() * pow(ellipse,-2.0) * ef / DELTA * ( U[i][j+1] + U[j][i+1] ) * N[j+1];
	  A(i,0) += MZ * (*H)[0].m() * pow(ellipse,-2.0) * ef * DELTA_0 * ( U[i][j+1] + U[j][i+1] ) * N[j+1];
	  A(i,1) += MZ * (*H)[0].m() * pow(ellipse,-2.0) * ef * DELTA_1 * ( U[i][j+1] + U[j][i+1] ) * N[j+1];

	  A(i,i) += MZ * (*H)[0].m() * pow(ellipse,-2.0) * ef / DELTA * M[j+1] * N[j+1];
	  A(i,j) += MZ * (*H)[0].m() * pow(ellipse,-2.0) * ef / DELTA * M[i+1] * N[j+1];

	  

	  // // b[i] += MZ * (*H)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * U[j][j+1]*N[i+1];

	  A(i,j) += MZ * (*H)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * M[j+1]*N[i+1];
	  
	  A(i,0) += MZ * (*H)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_0 * U[j][j+1]*N[i+1];
	  A(i,1) += MZ * (*H)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_1 * U[j][j+1]*N[i+1];
	}
  }



  //////////////////////////////////////////////////
  ///////////////// Other

 
  OtherEquation::OtherEquation(const ParamFile* paramfile) : Equation()
  {
    DataFormatHandler DFH;
    DFH.insert("gamma", &gamma, 0.);
    DFH.insert("rho", &rho, 0.);
    DFH.insert("rhow", &rhow, 0.);
    DFH.insert("Tref",&Tref,0.0);
    DFH.insert("Lref",&Lref,0.0);
    DFH.insert("Pstern",&Pstern,2.75e4);
    DFH.insert("ellipse",&ellipse,2.0);
    DFH.insert("C",&C,20.0);
    DFH.insert("Cdw",&Cdw,5.2e-3);
    DFH.insert("f",&f,0.0);
    DFH.insert("theta_w",&theta_w,0.0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(paramfile, "Equation");
    assert(rho>0);
    assert(Tref>0);
    assert(Lref>0);
    assert(Pstern>0);
    assert(ellipse>0);
    assert(C>0);
    //assert(f>0);

    MZ = 0.5*Tref*Tref * Pstern / rho / Lref / Lref;
    std::cout << "Mehlmann-Zahl " << MZ << std::endl;
  }
  
 
  /*-----------------------------------------*/

  // u[0], u[1] v
  
  // u[2] h    u[3] a
  void OtherEquation::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
   b[0] += U[0].m()*N.m();
    b[1] += U[1].m()*N.m();
    

    // Exponentialfaktor
    double ef = exp(-C*(1.0-(*HH)[1].m()));

    double dmin = DELTAMIN * Tref;
    double DELTAsquare =   dmin*dmin;
    
    DELTAsquare += 
      (1.0 + pow(ellipse,-2.0)) * ((*UU)[0].x()*(*UU)[0].x() + (*UU)[1].y()*(*UU)[1].y())
      + pow(ellipse,-2.0) * pow((*UU)[0].y() + (*UU)[1].x(),2.0)
      + 2.0 * (1.0-pow(ellipse,-2.0)) * (*UU)[0].x()*(*UU)[1].y();
      
    // DELTAsquare = std::max(DELTAsquare,dmin*dmin);
    
    double DELTA = sqrt(DELTAsquare);
    double a=((*UU)[0].x()-(*UU)[1].y())*((*UU)[0].x()-(*UU)[1].y())+((*UU)[0].y()+(*UU)[1].x())*((*UU)[0].y()+(*UU)[1].x());
   
    double eta= Pstern*(*HH)[0].m()*ef/(2.0*DELTA)/4.0;
    //    cout << Pstern << " "<< (*HH)[0].m() << " "
    //    	 << ef << " "<< DELTA << endl;
    //    cout << eta << endl;
    b[0] -= log(sqrt(a))/log(10)* N.m();
    //  b[0]-=sqrt(a)*N.m();
    b[1] -= log(eta*Tref)/log(10)*N.m();
    //b[1] -= eta*N.m();
    
  }
  /*-----------------------------------------*/
  
  void OtherEquation::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    A(0,0)+= M.m()*N.m();
    A(1,1)+= M.m()*N.m();
  }

}
