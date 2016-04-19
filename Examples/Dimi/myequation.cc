#include  "myequation.h"
#include  "filescanner.h"


extern double TIME;
extern double DT,DDD;
extern double DELTAMIN;


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

   // Wx.init("uvel.txt", 80,42,-20.0,0.0, 20.0,20.0);
    //Wy.init("vvel.txt", 80,42,-20.0,0.0, 20.0,20.0);
    

  }
  
  
  void MyEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
  {
    // Stabilisierung muss auch auf Dimensionslose groessen transformiert werden
    h_=h;
    
    double X = v.x()*Lref;
    double Y = v.y()*Lref;
    double Lx = 1.28e6;
    double Ly = 1.28e6;
    
    double Uw_x=0.1*(2*Y-Ly)/Ly;
       double Uw_y=(-0.1)*(2*X-Lx)/Lx;
    // double Uw_y=0.0;
    uwx= Uw_x *Tref/Lref;
    uwy= Uw_y *Tref/Lref;
    
  }
 
  /*-----------------------------------------*/

  // u[0], u[1] v
  
  // u[2] h    u[3] a
  
  void MyEquation::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    // Wasser-Tensor tau_w
       double WZ = rhow/rho * Cdw * Lref;
       //  double uwx = (*W)[0].m() ;
       //double uwy = (*W)[1].m() ;
   
    
    double vd_x = U[0].m() - uwx;
    double vd_y = U[1].m() - uwy;
    
    double vd = sqrt(vd_x*vd_x+vd_y*vd_y+1.e-8);
    double time_factor = 0.5*(1.0-cos(TIME * M_PI/10.0));
    if (TIME>10) time_factor = 1.0;

      
    b[0] += WZ * time_factor * vd *(vd_x * cos(theta_w) -vd_y*sin(theta_w))* N.m();
    b[1] += WZ * time_factor * vd *(vd_y * cos(theta_w) + vd_x*sin(theta_w))*N.m();
    
     
    //zustatz term Correolis 1.46 e -4 /s
    b[0] += -Tref*(*H)[0].m()*f*(vd_y)*N.m();
    b[1] +=  Tref*(*H)[0].m()*f*(vd_x)* N.m();
    

    // Zeitableitung
    b[0] += (*H)[0].m()*(U[0].m() - (*oldu)[0].m()) / DT * N.m();
    b[1] += (*H)[0].m()*(U[1].m() - (*oldu)[1].m()) / DT * N.m();

    // Exponentialfaktor
    double ef = exp(-C*(1.0-(*H)[1].m()));

    double dmin = DELTAMIN * Tref;
    double DELTAsquare =   dmin*dmin;
    
     
          // DELTA implizit
    DELTAsquare +=
    (1.0 + pow(ellipse,-2.0)) * (U[0].x()*U[0].x() + U[1].y()*U[1].y())
     + pow(ellipse,-2.0) * pow(U[0].y() + U[1].x(),2.0)
    + 2.0 * (1.0-pow(ellipse,-2.0)) * U[0].x()*U[1].y();
    
    
    /*     //delta explizit
    DELTAsquare += 
      (1.0 + pow(ellipse,-2.0)) * ((*extu)[0].x()*(*extu)[0].x() + (*extu)[1].y()*(*extu)[1].y())
      + pow(ellipse,-2.0) * pow((*extu)[0].y() + (*extu)[1].x(),2.0)
      + 2.0 * (1.0-pow(ellipse,-2.0)) * (*extu)[0].x()*(*extu)[1].y();
    */
   // DELTAsquare = std::max(DELTAsquare,dmin*dmin);
    
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
    b[1] += -MZ * (*H)[0].m() * ef * N.y();
  }
  /*-----------------------------------------*/
  
  void MyEquation::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {

    // Wasser-Tensor tau_w
    double WZ = rhow/rho * Cdw * Lref;
   
    //  double uwx = (*W)[0].m();
    //double uwy = (*W)[1].m();
    
    double vd_x = U[0].m() - uwx;
    double vd_y = U[1].m() - uwy;
    
    

    double vd = sqrt(vd_x*vd_x+vd_y*vd_y+1.e-8);
    double vd_0 = 1.0/(2.0*vd) * 2.0 * vd_x * M.m();
    double vd_1 = 1.0/(2.0*vd) * 2.0 * vd_y * M.m();
    double time_factor = 0.5*(1.0-cos(TIME * M_PI/10.0));
    if (TIME>10) time_factor = 1.0;
 
 
    A(0,0) += WZ * time_factor *  vd_0 *( vd_x * cos(theta_w) -vd_y*sin(theta_w) )* N.m();
    A(0,1) += WZ * time_factor * vd_1 *(vd_y * cos(theta_w) + vd_x*sin(theta_w)) *N.m();

    A(1,0) += WZ * time_factor * vd_0 * (vd_x * cos(theta_w) - vd_y*sin(theta_w))* N.m();
    A(1,1) += WZ * time_factor * vd_1 * (vd_y * cos(theta_w) + vd_x*sin(theta_w))*N.m();

    A(0,0) += WZ * time_factor * vd * M.m() * cos(theta_w) * N.m();
    A(1,1) += WZ * time_factor * vd * M.m() * cos(theta_w) * N.m();
    A(0,1) += WZ * time_factor * vd * (-M.m())*sin(theta_w) *N.m();
    A(1,0) += WZ * time_factor * vd * M.m()*sin(theta_w)   *N.m();


    // b[0] += -Tref*(*H)[0].m()*f*(U[1].m()-vd_y)*N.m();
    // b[1] +=  Tref*(*H)[0].m()*f*(U[0].m()-vd_x)* N.m();
    
 //zustatz term
  A(0,0) +=-Tref*(*H)[0].m()*f*M.m() * N.m();
  A(1,1) +=Tref*(*H)[0].m()*f*M.m() * N.m();
 
 
    
    // Zeitableitung
    A(0,0) +=  (*H)[0].m()* M.m() / DT * N.m();
    A(1,1) +=  (*H)[0].m()* M.m() / DT * N.m();

    double dmin = DELTAMIN*Tref;

    double DELTAsquare =  dmin*dmin;
        // implizit
 DELTAsquare +=
      (1.0 + pow(ellipse,-2.0)) * (U[0].x()*U[0].x() + U[1].y()*U[1].y())
     + pow(ellipse,-2.0) * pow(U[0].y() + U[1].x(),2.0)
      + 2.0 * (1.0-pow(ellipse,-2.0)) * U[0].x()*U[1].y();
    
    /*   
  // explizit
  DELTAsquare += 
      (1.0 + pow(ellipse,-2.0)) * ((*extu)[0].x()*(*extu)[0].x() + (*extu)[1].y()*(*extu)[1].y())
      + pow(ellipse,-2.0) * pow((*extu)[0].y() + (*extu)[1].x(),2.0)
      + 2.0 * (1.0-pow(ellipse,-2.0)) * (*extu)[0].x()*(*extu)[1].y();
  */
   // DELTAsquare = std::max(DELTAsquare,dmin*dmin);

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
    if (DELTAsquare==dmin*dmin)
      {
	DELTA_0 = 0.0;
	DELTA_1 = 0.0;
      }
    
    // explizit
    //  double DELTA_0 = 0.0;
    // double  DELTA_1 = 0.0;
    


 
    // Exponentialfaktor
    double ef   = exp(-C*(1.0-(*H)[1].m()));
    //    double ef_3 = exp(-C*(1.0-(*H)[1].m())) * (-C) * (-M.m());

    
    
   
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
  ///////////////// Transport



  TransportEquation::TransportEquation(const ParamFile* paramfile) : Equation()
  {
    DataFormatHandler DFH;
    DFH.insert("gamma", &gamma0, 0.);
    DFH.insert("Tref",&Tref,0.0);
    DFH.insert("tau", &tau0, 0.);
    DFH.insert("shock", &shock0, 0.);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(paramfile, "Equation");
  }
  
  
  void TransportEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
  {
    // Stabilisierung muss auch auf Dimensionslose groessen transformiert werden
    h_=h;

    double R=fabs(1.0/DT + ((*extu)[0].x()+(*extu)[1].y()));
    R = 1.0;
    tau=tau0 /(2. *sqrt((*extu)[0].m()*(*extu)[0].m()+(*extu)[1].m()*(*extu)[1].m())/h_ +R);
    shock = shock0 * h*h;
  }
 
  /*-----------------------------------------*/

  // u[0], u[1] v
  
  // u[2] h    u[3] a
  
  void TransportEquation::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    // Gleichung fuer h
    double test = N.m() + tau*((*extu)[0].m()*N.x()+(*extu)[1].m()*N.y());
    b[0] += (U[0].m() - (*oldh)[0].m()) / DT * test;
    b[0] += ((*extu)[0].m()*U[0].x()+(*extu)[1].m()*U[0].y()) *  test;
    b[0] += U[0].m() * ((*extu)[0].x()+(*extu)[1].y()) * test;
    
    
    // Gleichung für A
    b[1] += (U[1].m() - (*oldh)[1].m()) / DT      * test;
    b[1] += ((*extu)[0].m()*U[1].x()+(*extu)[1].m()*U[1].y())* test;
    b[1] += U[1].m() * ((*extu)[0].x()+(*extu)[1].y())       * test;
    
   
 
    
    // shock capt.
    double gh = sqrt(U[0].x()*U[0].x()+U[0].y()*U[0].y()+1.e-8);
    b[0] += shock * gh * (U[0].x()*N.x()+U[0].y()*N.y());
    
    // shock capt.
    double gA = sqrt(U[1].x()*U[1].x()+U[1].y()*U[1].y()+1.e-8);
    b[1] += shock * gA * (U[1].x()*N.x()+U[1].y()*N.y());
  }
  /*-----------------------------------------*/
  
  void TransportEquation::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
	
	// Stabilisierungsfaktor:
      // Gleichung fuer h
    double test = N.m() + tau*((*extu)[0].m()*N.x()+(*extu)[1].m()*N.y());
    
    //b[2] += (U[2].m() - (*old)[2].m()) / DT * test;
    A(0,0) += M.m() / DT * test;
 
    // b[2] += ((*extu)[0].m()*U[2].x()+(*extu)[1].m()*U[2].y()) *  test;
    A(0,0) += ((*extu)[0].m()*M.x()+(*extu)[1].m()*M.y()) *  test;

    //    b[2] += U[2].m() * ((*extu)[0].x()+(*extu)[1].y()) * test;
    A(0,0) += M.m() * ((*extu)[0].x()+(*extu)[1].y()) * test;


    //Gleichung für A
    //    b[3] += ((*H)[1].m() - (*old)[3].m()) / DT      * test;
    A(1,1) += M.m() / DT      * test;
    
    //b[3] += ((*extu)[0].m()*(*H)[1].x()+(*extu)[1].m()*(*H)[1].y())* test;
    A(1,1) += ((*extu)[0].m()*M.x()+(*extu)[1].m()*M.y())* test;
    
    //b[3] += (*H)[1].m() * ((*extu)[0].x()+(*extu)[1].y())       * test;
    A(1,1) += M.m() * ((*extu)[0].x()+(*extu)[1].y())       * test;


   
    
    
   
    //stabilisierung    
    //shock capt.
    double gh = sqrt(U[0].x()*U[0].x()+U[0].y()*U[0].y()+1.e-8);
    double Dgh = 1.0/gh *(U[0].x()*M.x() + U[0].y()*M.y());
    A(0,0) += shock * Dgh * (U[0].x()*N.x()+U[0].y()*N.y());
    A(0,0) += shock * gh * (M.x()*N.x()+M.y()*N.y());
    

    
    // // stabilisierung
    double gA = sqrt(U[1].x()*U[1].x()+U[1].y()*U[1].y()+1.e-8);
    double DgA = 1.0/gA *(U[1].x()*M.x() + U[1].y()*M.y());
    A(1,1) += shock * DgA * (U[1].x()*N.x()+U[1].y()*N.y());
    A(1,1) += shock * gA * (M.x()*N.x()+M.y()*N.y());
  }




  






//////////////////////////////////////////////////
  ///////////////// Other



  OtherEquation::OtherEquation(const ParamFile* paramfile) : Equation()
  {
    DataFormatHandler DFH;
    //    DFH.insert("theta_w",&theta_w,0.0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(paramfile, "Equation");
  }
  
 
  /*-----------------------------------------*/


  
  void OtherEquation::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    b[0]+=U[0].m()*N.m();
 
    
    b[1] +=U[1].m()*N.m();
   

 
  }
  /*-----------------------------------------*/
  
  void OtherEquation::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    

    A(0,0)+= M.m()*N.m();
    A(1,1) += M.m()*N.m();
  
   
  }



 MyEq::MyEq(const ParamFile* paramfile) : Equation()
  {
    DataFormatHandler DFH;
    DFH.insert("gamma", &gamma0, 0.);
    DFH.insert("Tref",&Tref,0.0);
    DFH.insert("tau", &tau0, 0.);
    DFH.insert("shock", &shock0, 0.);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(paramfile, "Equation");
  }
  
  
  void MyEq::point(double h, const FemFunction &U, const Vertex2d &v) const
  {
    // Stabilisierung muss auch auf Dimensionslose groessen transformiert werden
    h_=h;

    
  }
 
  /*-----------------------------------------*/

 
  
  void MyEq::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const 
  {
    // Gleichung fuer h
     b[0] +=   (U[0].x()+U[0].y())*N.m();
    

    // Gleichung für A
     b[1] +=   (U[1].x()+U[1].y())*N.m();
    

    
    
    
  }
  /*-----------------------------------------*/
  
  void MyEq::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {
    //    Das ist die  C- Matrix phi_i * nabla phi_j 
    A(0,0) += N.m() * M.x();
    A(1,1) += N.m() * M.y();
  }














}
 /*
    double div=((*oldu)[0].x()+(*oldu)[1].y());
    double d_t1= (*oldu)[0].m()*(*oldh)[0].x()+(*oldu)[1].m()*(*oldh)[0].y() ; 
    double d_t2= (*oldh)[0].m() * div;
    double d_t= d_t1+d_t2;
    
    b[0] += U[0].m() * N.m();
    
    
    b[0] += - (*oldh)[0].m() * N.m();
     
    b[0] +=   (U[0].x()+U[0].y())*N.m();

    double dtV1 = ((*oldu)[0].m()-(*oldoldu)[0].m())/DT;
    double dtV2 = ((*oldu)[1].m()-(*oldoldu)[1].m())/DT;

    b[0] += DT*DT/2.0 *   (-dtV1 * (*oldh)[0].m()+ (*oldu)[0].m()*d_t )   * N.x();
    b[0] += DT*DT/2.0 *   (-dtV2 * (*oldh)[0].m()+ (*oldu)[1].m()*d_t )   * N.y();
     */
 /*
    double d_At1= (*oldu)[0].m()*(*oldh)[1].x()+(*oldu)[1].m()*(*oldh)[1].y() ; 
    double d_At2= (*oldh)[1].m() * div;
    double d_At= d_At1+d_At2;
    
    b[1] += U[1].m() * N.m();
    
    b[1] += - (*oldh)[1].m() * N.m();
    b[1] +=   DT*d_At*N.m();

   

    b[1] += DT*DT/2.0 *   (-dtV1 * (*oldh)[1].m()+(*oldu)[0].m()*d_At )   * N.x();
    b[1] += DT*DT/2.0 *   (-dtV2 * (*oldh)[1].m()+(*oldu)[1].m()*d_At )   * N.y();
     */
    


/*-----------------------------------------*/
