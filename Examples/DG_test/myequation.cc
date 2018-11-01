#include "myequation.h"


extern double DT,DELTAMIN,DDD;
extern double CHIGAUSS;
extern bool Jump;

namespace Gascoigne
{
  /*----------------------------------------------------------------------------*/



  
  double Ice::delta(const FemFunction& U) const 
  {    
    double dmin = DELTAMIN*Tref;
    
    // DELTA implizit
    double DELTAsquare =
      (1.0 + pow(ellipse,-2.0)) * (U[0].x()*U[0].x() + U[1].y()*U[1].y())
      + pow(ellipse,-2.0) * pow(U[0].y() + U[1].x(),2.0)
      + 2.0 * (1.0-pow(ellipse,-2.0)) * U[0].x()*U[1].y();
    
    DELTAsquare += dmin * dmin;
    return sqrt(DELTAsquare);
  }
  
  double Ice::delta_0(const FemFunction& V, const TestFunction& M) const
  {
    double DELTAsquare_0 =
      (1.0 + pow(ellipse,-2.0)) * (2.0*V[0].x()*M.x())
      + pow(ellipse,-2.0) * 2.0 * (V[0].y() + V[1].x()) * M.y()
      + 2.0 * (1.0-pow(ellipse,-2.0)) * M.x()*V[1].y();
    return  -0.5 * pow(delta(V),-3) * DELTAsquare_0*DDD;
  }
  
  double Ice::delta_1(const FemFunction& V, const TestFunction& M) const
  {
    double DELTAsquare_1 =
      (1.0 + pow(ellipse,-2.0)) * (2.0 * V[1].y()*M.y())
      + pow(ellipse,-2.0) * 2.0*(V[0].y() + V[1].x()) * M.x()
      + 2.0 * (1.0-pow(ellipse,-2.0)) * V[0].x()*M.y();
    return -0.5 * pow(delta(V),-3) * DELTAsquare_1*DDD;
  }

  double Ice::SIGMA_ij(const FemFunction& V, const TestFunction& M, const TestFunction& N, int i, int j) const
  {
    abort();
  }
  
  /*----------------------------------------------------------------------------*/
  

  // Transport
  
  
  
  

  TransportEquation::TransportEquation(const ParamFile* pf) : Equation()
  {
    DataFormatHandler DFH;
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
    //   V.resize(2);
    
  }
  
  /*----------------------------------------------------------------------------*/

  void TransportEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
  {/*
     V[0].m() = -v.y(); V[0].x() = -0.0; V[0].y() = -1.0;
     V[1].m() =  v.x(); V[1].x() =  1.0; V[1].y() =  0.0;
     V[0]*=2.0*M_PI;
     V[1]*=2.0*M_PI;
   */
  }

  void TransportEquation::Form(VectorIterator b,
			       const FemFunction &U,
			       const TestFunction &N) const
  {
    // 2te-Gleichung: Null setzen
      
    if(Jump)
      { 
	//Zeit
	b[0] += (U[0].m() -(*oldH)[0].m()) * N.m();
	b[1] += (U[1].m() -(*oldH)[1].m()) * N.m();
      }
    // ganz einfache stabilisierung...
    b[0] += 0.0001 * DT * (U[0].x()*N.x() + U[0].y()*N.y());
    b[1] += 0.0001 * DT * (U[1].x()*N.x() + U[1].y()*N.y());
    
    // div (vH)
    b[0] += DT * ( (*V)[0].m()*U[0].x() + (*V)[1].m()*U[0].y() ) * N.m();
    b[0] += DT * ( (*V)[0].x() + (*V)[1].y() ) * U[0].m() * N.m();
    
    b[1] += DT * ( (*V)[0].m()*U[1].x() + (*V)[1].m()*U[1].y() ) * N.m();
    b[1] += DT * ( (*V)[0].x() + (*V)[1].y() ) * U[1].m() * N.m();
  }
  

  /*----------------------------------------------------------------------------*/

  void TransportEquation::Matrix(EntryMatrix &A,
				 const FemFunction &U,
				 const TestFunction &M,
				 const TestFunction &N) const
  {
    

    //Zeit
    A(0,0) += M.m() * N.m() ;
    A(1,1) += M.m() * N.m() ;
    
    //stabilisierung

    A(0,0) += 0.0001 * DT * (M.x()*N.x() + M.y()*N.y());
    A(1,1) += 0.0001 * DT * (M.x()*N.x() + M.y()*N.y());
    
    // div (vH)
    A(0,0) += DT * ( (*V)[0].m()*M.x() + (*V)[1].m()*M.y() ) * N.m();
    A(0,0) += DT * ( (*V)[0].x() + (*V)[1].y() ) * M.m() * N.m();
    
    A(1,1) += DT * ( (*V)[0].m()*M.x() + (*V)[1].m()*M.y() ) * N.m();
    A(1,1) += DT * ( (*V)[0].x() + (*V)[1].y() ) * M.m() * N.m();
  }



  MyDualEquation::MyDualEquation(const ParamFile* pf) : Equation(), Ice(pf)
  {
  }

  /*----------------------------------------------------------------------------*/

  void MyDualEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
  {/*
     V[0].zero();
     V[0].m() = -v.y(); V[0].x() = -0.0; V[0].y() = -1.0;
     V[1].m() =  v.x(); V[1].x() =  1.0; V[1].y() =  0.0;
     V[0]*=2.0*M_PI;
     V[1]*=2.0*M_PI;
   */
  
    double X = v.x()*Lref;
    double Y = v.y()*Lref;
    double Lx = 1.28e6;
    double Ly = 1.28e6;

    double Uw_x=0.1*(2*Y-Ly)/Ly;
    double Uw_y=(-0.1)*(2*X-Lx)/Lx;

    uwx= Uw_x *Tref/Lref;
    uwy= Uw_y *Tref/Lref; 
  
 
  }
 

  /*----------------------------------------------------------------------------*/

  void MyDualEquation::Form(VectorIterator b,
			    const FemFunction &Z,
			    const TestFunction &N) const
  {
  
    
    if(Jump)
      {
	//Zeit
	b[0] += ( Z[0].m() -(*nextZ)[0].m() ) * N.m();
	b[1] += ( Z[1].m() -(*nextZ)[1].m() ) * N.m();
    
	// Kopplung Zeit
	b[0] += N.m() * rho*  ( (*nextV)[0].m()-(*V)[0].m() ) * (*nextQ)[0].m();
	b[0] += N.m() * rho * ( (*nextV)[1].m()-(*V)[1].m() ) * (*nextQ)[1].m();
      }
    
    // Stabilisierung
    b[0] += 0.0001*DT * (Z[0].x()*N.x() + Z[0].y()*N.y());
    b[1] += 0.0001*DT * (Z[1].x()*N.x() + Z[1].y()*N.y());
    
    // //// div v*h +v grad h
    b[0]+=DT * ( (*V)[0].x() + (*V)[1].y() ) * Z[0].m() * N.m();
    b[0]+=DT * (N.x()*(*V)[0].m()+N.y()*(*V)[1].m()) * Z[0].m();
    
    b[1]+=DT * ( (*V)[0].x() + (*V)[1].y() ) * Z[1].m() * N.m();
    b[1]+=DT * (N.x()*(*V)[0].m()+N.y()*(*V)[1].m()) * Z[1].m();
    
    // kopplung tensor
    
    // Exponentialfaktor
        
    double ef = exp(-C*(1.0-(*H)[1].m()));
   
    double DELTA = delta(*nextV);
   
    ///////  (sigma, nabla phi)       
    // nabla u, nabla phi
    // nable u^T, nabla phi

    // Kopplung nach H
    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
        {
          b[0] += DT*MZ * N.m() * pow(ellipse,-2.0) * ef / DELTA * ( (*nextV)[i][j+1] + (*nextV)[j][i+1] ) * (*nextQ)[i][j+1];
          b[0] += DT*MZ * N.m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * (*nextV)[j][j+1]*(*nextQ)[i][i+1];
        }
    // druck-Term
    b[0] += -MZ * DT*N.m() * ef * (*nextQ)[0].x();
    b[0] += -MZ * DT*N.m() * ef * (*nextQ)[1].y();

    // Kopplung nach A
    double ef_A = C * N.m() * ef;
    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
        {
          b[1] += DT*MZ * (*H)[0].m() * pow(ellipse,-2.0)       * ef_A / DELTA * ( (*nextV)[i][j+1] + (*nextV)[j][i+1] ) * (*nextQ)[i][j+1];
          b[1] += DT*MZ * (*H)[0].m() * (1.0-pow(ellipse,-2.0)) * ef_A / DELTA * (*nextV)[j][j+1]* (*nextQ)[i][i+1];
        }
    // // ef, siehe oben, exponentialfaktor
    b[1] += -MZ * DT*(*H)[0].m() * ef_A * (*nextQ)[0].x();
    b[1] += -MZ * DT*(*H)[0].m() * ef_A * (*nextQ)[1].y();

    
    double WZ = rhow/rho * Cdw * Lref;
   
    double vd_x = (*V)[0].m() - uwx;
    double vd_y = (*V)[1].m() - uwy;
   
    // Korreolistermkopplung in H
    b[0] += -Tref*DT*N.m()*f*(vd_y) * (*nextQ)[0].m();
    b[0] +=  Tref*DT*N.m()*f*(vd_x) * (*nextQ)[1].m(); 
    
    
  }
  
  void MyDualEquation::Matrix(EntryMatrix &A,
			      const FemFunction &U,
			      const TestFunction &M,
			      const TestFunction &N) const
  {
    //Zeit
    A(0,0) += M.m() * N.m();
    A(1,1) += M.m() * N.m();
    
    // ganz einfache stabilisierung...
    A(0,0) += 0.0001*DT * (M.x()*N.x() + M.y()*N.y());
    A(1,1) += 0.0001*DT * (M.x()*N.x() + M.y()*N.y());
    
    //div vh 
    
    A(0,0)+=DT*((*V)[0].x()+(*V)[1].y())*M.m()*N.m();
    A(0,0)+=DT*(N.x()*(*V)[0].m()+N.y()*(*V)[1].m())*M.m();
        
    A(1,1)+=DT*((*V)[0].x()+(*V)[1].y())*M.m()*N.m();
    A(1,1)+=DT*(N.x()*(*V)[0].m()+N.y()*(*V)[1].m())*M.m();
    
  }
   

  /// Burger


  MyBurgerEquation::MyBurgerEquation(const ParamFile* pf) : Equation(), Ice(pf)
  {
    assert(rho>0);
    assert(Tref>0);
    assert(Lref>0);
    assert(Pstern>0);
    assert(ellipse>0);
    assert(C>0);
    assert(f>0);
  }

  /*----------------------------------------------------------------------------*/

  void MyBurgerEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
  {
      
    // H[0].zero();
    double X = v.x()*Lref;
    double Y = v.y()*Lref;
    double Lx = 1.28e6;
    double Ly = 1.28e6;

    double Uw_x=0.1*(2*Y-Ly)/Ly;
    double Uw_y=(-0.1)*(2*X-Lx)/Lx;

    uwx= Uw_x *Tref/Lref;
    uwy= Uw_y *Tref/Lref; 
    
    
    
    // H[0].m() = 0.5; H[0].x() = 0.0; H[0].y() =  0.0;
    // H[1].m() = 0.5; H[1].x() = 0.0; H[1].y() =  0.0;
 

    
    
    // Spaeter erstezen durch oldH 
    
    
  }

  /*----------------------------------------------------------------------------*/
  
  void MyBurgerEquation::Form(VectorIterator b,
			      const FemFunction &V,
			      const TestFunction &N) const
  {
    if(Jump)
      {
	b[0] += rho*(1  * (*oldH)[0].m() ) * (V[0].m() - (*oldV)[0].m())  * N.m();
	b[1] += rho*(1  * (*oldH)[0].m() ) * (V[1].m() - (*oldV)[1].m())  * N.m();
      }

    double WZ = rhow/rho * Cdw * Lref;
   
    double vd_x = V[0].m() - uwx;
    double vd_y = V[1].m() - uwy;
   
    double vd = sqrt(vd_x*vd_x+vd_y*vd_y+1.e-8);
     
    b[0] += DT*WZ * vd *(vd_x * cos(theta_w) - vd_y*sin(theta_w)) * N.m();
    b[1] += DT*WZ * vd *(vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
   
    //zustatz term Correolis 1.46 e -4 /s
    b[0] += -Tref*DT*(*oldH)[0].m()*f*(vd_y)*N.m();
    b[1] +=  Tref*DT*(*oldH)[0].m()*f*(vd_x)* N.m(); 
      
    // Exponentialfaktor
    double ef = exp(-C*(1.0-(*oldH)[1].m()));

    double DELTA = delta(V);
   
    ///////  (sigma, nabla phi)       
    // nabla u, nabla phi
    // nable u^T, nabla phi
    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
        {
          b[i] += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * ( V[i][j+1] + V[j][i+1] ) * N[j+1];
          b[i] += DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * V[j][j+1]*N[i+1];
        }
   
    // // ef, siehe oben, exponentialfaktor
    b[0] += -MZ * DT*(*oldH)[0].m() * ef * N.x();
    b[1] += -MZ * DT*(*oldH)[0].m() * ef * N.y();
  }

  /*----------------------------------------------------------------------------*/

  void MyBurgerEquation::Matrix(EntryMatrix &A,
				const FemFunction &V,
				const TestFunction &M,
				const TestFunction &N) const
  {
    //Zeit,
    A(0,0) += rho*(1 * (*oldH)[0].m() ) * M.m() * N.m();
    A(1,1) += rho*(1 * (*oldH)[0].m() ) * M.m() * N.m();
    
    // Wasser-Tensor tau_w
    double WZ = rhow/rho * Cdw * Lref;
   
    double vd_x = V[0].m() - uwx;
    double vd_y = V[1].m() - uwy;

    double vd = sqrt(vd_x*vd_x+vd_y*vd_y+1.e-8);
    double vd_0 = 1.0/(2.0*vd) * 2.0 * vd_x * M.m();
    double vd_1 = 1.0/(2.0*vd) * 2.0 * vd_y * M.m();

    //b[0] += WZ * vd   *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * N.m();
    A(0,0) += DT*WZ * vd_0 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * N.m();
    A(0,1) += DT*WZ * vd_1 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * N.m();
    A(0,0) += DT*WZ * vd   *(M.m() * cos(theta_w)                    ) * N.m();
    A(0,1) += DT*WZ * vd   *(                    -M.m()*sin(theta_w) ) * N.m();
 
    //b[1] += WZ * vd   * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
    A(1,0) += DT*WZ * vd_0 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
    A(1,1) += DT*WZ * vd_1 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
    A(1,1) += DT*WZ * vd   * (M.m() * cos(theta_w)                   ) * N.m();
    A(1,0) += DT*WZ * vd   * (                     M.m()*sin(theta_w)) * N.m();
                          

    //zustatz term
    //b[0] += -Tref*(*H)[0].m()*f*(vd_y)*N.m();
    A(0,1) += -DT*Tref*(*oldH)[0].m()*f*M.m() * N.m();
    //b[1] +=  Tref*(*H)[0].m()*f*(vd_x)* N.m();
    A(1,0) +=  DT*Tref*(*oldH)[0].m()*f*M.m() * N.m();
    
    double DELTA   = delta(V);
    double DELTA_0 = delta_0(V,M);
    double DELTA_1 = delta_1(V,M);
 
    // Exponentialfaktor
    double ef   = exp(-C*(1.0-(*oldH)[1].m()));
   
    // nabla u, nabla phi
    // nable u^T, nabla phi
    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
	{
          //          b[i] += MZ * (*H)[0].m() * pow(ellipse,-2.0) * ef / DELTA * ( U[i][j+1] + U[j][i+1] ) * N[j+1];
          A(i,0) += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef * DELTA_0 * ( V[i][j+1] + V[j][i+1] ) * N[j+1];
          A(i,1) += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef * DELTA_1 * ( V[i][j+1] + V[j][i+1] ) * N[j+1];

          A(i,i) += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * M[j+1] * N[j+1];
          A(i,j) += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * M[i+1] * N[j+1];

	  // // b[i] += MZ * (*H)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * U[j][j+1]*N[i+1];

          A(i,j) += DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * M[j+1]*N[i+1];
          
          A(i,0) += DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_0 * V[j][j+1]*N[i+1];
          A(i,1) += DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_1 * V[j][j+1]*N[i+1];
        } 

  }



  //dual burger problem
  /*----------------------------------------------------------------------------*/
  MyBurgerDualEquation::MyBurgerDualEquation(const ParamFile* pf) : Equation(), Ice(pf)
  {
    assert(f>0);
    assert(Tref>0); 
    assert(rho>0);
    assert(Tref>0);
    assert(Lref>0);
    assert(Pstern>0);
    assert(ellipse>0);
    assert(C>0);
  }


  /*----------------------------------------------------------------------------*/  
  void MyBurgerDualEquation::point(double h, const FemFunction &U, const Vertex2d &v) const
  {
    //H[0].zero();
    double X = v.x()*Lref;
    double Y = v.y()*Lref;
    double Lx = 0.5e6;
    double Ly = 0.5e6;

    double Uw_x=0.1*(2*Y-Ly)/Ly;
    double Uw_y=(-0.1)*(2*X-Lx)/Lx;

    uwx= Uw_x *Tref/Lref/10;
    uwy= Uw_y *Tref/Lref/10;
    
    // H[0].m() = 0.5; H[0].x() =  0.0; H[0].y() =  0.0;
    // H[1].m() = 0.5; H[1].x() =  0.0; H[1].y() =  0.0;
  }

  /*----------------------------------------------------------------------------*/

  
  void MyBurgerDualEquation::Form(VectorIterator b,
				  const FemFunction &Q,
				  const TestFunction &N) const
  {  
    if(Jump)
      {
	//Zeit
	b[0] += rho*(1.0  * (*oldH)[0].m() ) * Q[0].m() * N.m();
	b[1] += rho*(1.0  * (*oldH)[0].m() ) * Q[1].m() * N.m();
	
	b[0] += -rho*(1.0  * (*H)[0].m() ) * (*nextQ)[0].m() * N.m();
	b[1] += -rho*(1.0  * (*H)[0].m() ) * (*nextQ)[1].m() * N.m();
      }
    // Correolisterm
    b[1] += -DT*Tref*(*oldH)[0].m()*f*N.m()*Q[0].m();
    b[0] +=  DT*Tref*(*oldH)[0].m()*f*N.m()*Q[1].m();
     
     
    /*
      A(0,0) += DT*WZ * vd_0 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * N.m();
      A(0,1) += DT*WZ * vd_1 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * N.m();
      A(0,0) += DT*WZ * vd   *(M.m() * cos(theta_w)                    ) * N.m();
      A(0,1) += DT*WZ * vd   *(                    -M.m()*sin(theta_w) ) * N.m();
 
      //b[1] += WZ * vd   * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
      A(1,0) += DT*WZ * vd_0 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
      A(1,1) += DT*WZ * vd_1 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * N.m();
      A(1,1) += DT*WZ * vd   * (M.m() * cos(theta_w)                   ) * N.m();
      A(1,0) += DT*WZ * vd   * (                     M.m()*sin(theta_w)) * N.m();
   */
    
    double WZ = rhow/rho * Cdw * Lref;
   
   
    double vd_x = (*V)[0].m() - uwx;
    double vd_y = (*V)[1].m() - uwy;

    double vd = sqrt(vd_x*vd_x+vd_y*vd_y+1.e-8);
    
    
    double vd_0 = 1.0/(2.0*vd) * 2.0 * vd_x * N.m();
    double vd_1 = 1.0/(2.0*vd) * 2.0 * vd_y * N.m();
  
    //Tau wasser
    b[0] += DT*WZ * vd_0 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * Q[0].m();
    b[1] += DT*WZ * vd_1 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * Q[0].m();
    b[0] += DT*WZ * vd   *(N.m() * cos(theta_w)                    ) * Q[0].m();
    b[1] += DT*WZ * vd   *(                    -N.m()*sin(theta_w) ) * Q[0].m();
 
    b[0] += DT*WZ * vd_0 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * Q[1].m();
    b[1] += DT*WZ * vd_1 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * Q[1].m();
    b[1] += DT*WZ * vd   * (N.m() * cos(theta_w)                   ) * Q[1].m();
    b[0] += DT*WZ * vd   * (                     N.m()*sin(theta_w)) * Q[1].m();
    
        
   
    double DELTA = delta(*V);
   
    double DELTA_0 = delta_0(*V,N);
    double DELTA_1 = delta_1(*V,N);

    double ef   = exp(-C*(1.0-(*oldH)[1].m()));
    
    // nabla u, nabla phi
    // nable u^T, nabla phi
    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
        {
	  // A(i,0) += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef * DELTA_0 * ( V[i][j+1] + V[j][i+1] ) * N[j+1];
	  b[0]+= DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef * DELTA_0 * ( (*V)[i][j+1]+(*V)[j][i+1] ) * Q[i][j+1];
            
	  //A(i,1) += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef * DELTA_1 * ( V[i][j+1] + V[j][i+1] ) * N[j+1];
           
	  b[1]+= DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef * DELTA_1 * ( (*V)[i][j+1] + (*V)[j][i+1]) * Q[i][j+1];
           
          
	  //b[i] += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * ( V[i][j+1] + V[j][i+1] ) * N[j+1];
	  //A(i,i) += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * M[j+1] * N[j+1];
          //A(i,j) += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * M[i+1] * N[j+1];
           
          b[i] += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * N[j+1] * Q[i][j+1];
          b[j] += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * N[i+1] * Q[i][j+1];

        
	  //b[i] += MZ * (*H)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * U[j][j+1]*N[i+1];
          //A(i,j) += DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * M[j+1]*N[i+1];
          b[j]+= DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * N[j+1]*Q[i][i+1];
          
	  // A(i,0) += DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_0 * V[j][j+1]*N[i+1];
          b[0]+=DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_0 * (*V)[j][j+1]*Q[i][i+1];
          
          //A(i,1) += DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_1 * V[j][j+1]*N[i+1];
          b[1]+=DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_1 * (*V)[j][j+1]*Q[i][i+1];
        } 
              
    //Kopplung H
    // div (phi) H z
    b[0] += DT * N.x() * (*H)[0].m() * (*Z)[0].m();
    b[1] += DT * N.y() * (*H)[0].m() * (*Z)[0].m();
    // phi * nabla H z
    b[0] += DT * N.m() * (*H)[0].x() * (*Z)[0].m();
    b[1] += DT * N.m() * (*H)[0].y() * (*Z)[0].m();  
    
    //Kopplung A
    // div (phi) H z
    b[0] += DT * N.x() * (*H)[1].m() * (*Z)[1].m();
    b[1] += DT * N.y() * (*H)[1].m() * (*Z)[1].m();
    // phi * nabla H z
    b[0] += DT * N.m() * (*H)[1].x() * (*Z)[1].m();
    b[1] += DT * N.m() * (*H)[1].y() * (*Z)[1].m(); 
  }
  
  void MyBurgerDualEquation::Matrix(EntryMatrix &A,
				    const FemFunction &Z,
				    const TestFunction &M,
				    const TestFunction &N) const
  {
    //Zeit
    A(0,0) += rho*(*oldH)[0].m()* M.m()*N.m();
    A(1,1) += rho*(*oldH)[0].m()* M.m()*N.m();
    
    // b[1] += -DT*Tref*(*oldH)[0].m()*f*N.m()*Q[0].m();
    // b[0] +=  DT*Tref*(*oldH)[0].m()*f*N.m()*Q[1].m();
    A(1,0) += -DT*Tref*(*oldH)[0].m()*f*N.m()*M.m();
    A(0,1) +=  DT*Tref*(*oldH)[0].m()*f*N.m()*M.m();
    
    
    double WZ = rhow/rho * Cdw * Lref;
   
   
    double vd_x = (*V)[0].m() - uwx;
    double vd_y = (*V)[1].m() - uwy;

    double vd = sqrt(vd_x*vd_x+vd_y*vd_y+1.e-8);
    
    
    double vd_0 = 1.0/(2.0*vd) * 2.0 * vd_x * N.m();
    double vd_1 = 1.0/(2.0*vd) * 2.0 * vd_y * N.m();
  
    
    //Correolisterm
    
    /*b[0] += DT*WZ * vd_0 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * Q[0].m();
      b[1] += DT*WZ * vd_1 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * Q[0].m();
      b[0] += DT*WZ * vd   *(N.m() * cos(theta_w)                    ) * Q[0].m();
      b[1] += DT*WZ * vd   *(                    -N.m()*sin(theta_w) ) * Q[0].m();
 
      b[0] += DT*WZ * vd_0 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * Q[1].m();
      b[1] += DT*WZ * vd_1 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * Q[1].m();
      b[1] += DT*WZ * vd   * (N.m() * cos(theta_w)                   ) * Q[1].m();
      b[0] += DT*WZ * vd   * (                     N.m()*sin(theta_w)) * Q[1].m();
    
    
    */
    

    A(0,0) += DT*WZ * vd_0 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * M.m();
    A(1,0) += DT*WZ * vd_1 *(vd_x * cos(theta_w) - vd_y*sin(theta_w) ) * M.m();
    A(0,0) += DT*WZ * vd   *(N.m() * cos(theta_w)                    ) * M.m();
    A(1,0) += DT*WZ * vd   *(                    -N.m()*sin(theta_w) ) * M.m();
 
    A(0,1) += DT*WZ * vd_0 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * M.m();
    A(1,1) += DT*WZ * vd_1 * (vd_y * cos(theta_w) + vd_x*sin(theta_w)) * M.m();
    A(1,1) += DT*WZ * vd   * (N.m() * cos(theta_w)                   ) * M.m();
    A(0,1) += DT*WZ * vd   * (                     N.m()*sin(theta_w)) * M.m();
    
    //------------------------------------------------------
    double DELTA = delta(*V);
   
    double DELTA_0 = delta_0(*V,N);
    double DELTA_1 = delta_1(*V,N);

 
    // Exponentialfaktor
    double ef   = exp(-C*(1.0-(*oldH)[1].m()));
   
    // nabla u, nabla phi
    // nable u^T, nabla phi
    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
    
	{
          
            
          //  b[0]+= DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef * DELTA_0 * ( (*V)[i][j+1] (V*)[j][i+1] ) * Q[i][j+1];
            
	  A(0,i)+= DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef * DELTA_0 * ( (*V)[i][j+1] +(*V)[j][i+1] ) * M[j+1];
           
           
	  //b[1]+= DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef * DELTA_1 * ( (V*)[i][j+1] + (V*)[j][i+1]) * Q[i][j+1];
	  A(1,i)+= DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef * DELTA_1 * ( (*V)[i][j+1] + (*V)[j][i+1]) * M[j+1];
          

           
	  // b[i] += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * N[j+1] * Q[i][j+1];
           
	  A(i,i) += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * N[j+1] * M[j+1];
          
          
	  // b[j] += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * N[i+1] * Q[i][j+1];
            
	  A(j,i) += DT*MZ * (*oldH)[0].m() * pow(ellipse,-2.0) * ef / DELTA * N[i+1] *M[j+1];


     
          //b[j]+= DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * N[j+1]*Q[i][i+1];
	  A(j,i)+= DT*MZ *(*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef / DELTA * N[j+1]*M[i+1];
          
       
	  // b[0]+=DT*MZ *  (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_0 * V[j][j+1]*Q[i][i+1];
          A(0,i)+=DT*MZ *(*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_0 * (*V)[j][j+1]*M[i+1];
          
          //b[1]+= DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_1 * V[j][j+1]*Q[i][i+1];
          A(1,i)+= DT*MZ * (*oldH)[0].m() * (1.0-pow(ellipse,-2.0)) * ef * DELTA_1 * (*V)[j][j+1]*M[i+1];
          
        } 
  }

  
 
  


 
}

