#include "loop.h"
#include "gascoignemesh2d.h"
#include <time.h>
#include  "backup.h"
#include  "stdmultilevelsolver.h"
#include  "sparseblockmatrix.h"
#include  "fmatrixblock.h"
#include  "compose_name.h"
#include  "dwrfem.h"
#include "pi.h"
using namespace Gascoigne;
using namespace std;
 


double TIME, DT,DTM, CHIGAUSS;


string Loop::SolvePrimalSingle(VectorInterface& u, VectorInterface& f, string name)
{
  
  GetMultiLevelSolver()->GetSolver()->Zero(f);

  // Rechte Seite mit Gauss-Regel
  double wgt = 0.5;
  double x1 = TIME-0.5*DT - 0.5*sqrt(1.0/3.0) * DT;
  double x2 = TIME-0.5*DT + 0.5*sqrt(1.0/3.0) * DT;
  double SAVETIME = TIME;
  TIME = x1;
  GetMultiLevelSolver()->GetSolver()->Rhs(f,wgt * DT);
  // GetMultiLevelSolver()->GetSolver()->Rhs(f,wgt );
  TIME = x2;
  GetMultiLevelSolver()->GetSolver()->Rhs(f,wgt * DT);
  //GetMultiLevelSolver()->GetSolver()->Rhs(f,wgt);
  TIME = SAVETIME;
  
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
  GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(u);
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);

  string status = GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());

  Output(u,name);

  return status;
}


void Loop::SolvePrimalProblem(vector<GlobalVector> &Utotal, nvector<double>& Jtotal, VectorInterface& u, VectorInterface& oldu, VectorInterface& f)
{
  for (_iter=1; _iter<=_niter; _iter++)
    {
      TIME += DT;
      
      cout << "\n Zeitschritt " << _iter << " " << TIME-DT << " -> " 
	   << TIME<< "  [" << DT << "]" << endl;
      GetMultiLevelSolver()->GetSolver()->Equ(oldu,1.0, u);
      
      GetMultiLevelSolver()->AddNodeVector("oldu", oldu);
      string res = SolvePrimalSingle(u,f,"Results/u");
      Utotal.push_back(GetMultiLevelSolver()->GetSolver()->GetGV(u));
      GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,_iter);
      // Funktional
      nvector<double> functionals = Functionals(u,f);
      Jtotal.push_back(functionals[0]);
      GetMultiLevelSolver()->DeleteNodeVector("oldu");
     
      
    }
}


void Loop::BoxInt(GlobalVector& avg,const vector<GlobalVector>& U, int start, int stopp)
{
  assert(start<U.size());  assert(stopp<U.size());
  avg.zero();
  for (int l=start ; l<=stopp;++l)
    avg.add(DT,U[l]);
  avg*=1/DTM;
}
void Loop::TrapezInt(GlobalVector& avg, const vector<GlobalVector>& U, int start, int stopp)
{
  assert(start<U.size());  assert(stopp<U.size());
  avg.zero();
  avg.add(0.5*DT,U[start]);
  for (int l=start+1 ; l<=stopp-1;++l)
    avg.add(DT,U[l]);
  avg.add(0.5*DT,U[stopp]);
  avg*=1/DTM;
}

void Loop::SolveDualProblem(VectorInterface& z, VectorInterface& f,double DTM)
{
  GetMultiLevelSolver()->SetProblem("dp");
 

  DoubleVector test1;
  string res_dual = Solve(z,f,"Results/z");
};


class DWRInitialRhs : public virtual DomainRightHandSide
{
  std::string GetName() const {return "DWRInitialRightHandSide";}
  int GetNcomp() const {return 1; }
  
  mutable FemFunction* U;
  
  void SetFemData(FemData& q) const 
  {
    assert(q.find("U") != q.end() );
    U = &q["U"]; 
  }
  
  
  void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
  {
    double x = v.x();
    double y = v.y();
    double t = 0.0;
    double u0 =  0.1e1 / (0.1e1 + 0.50e2 * pow(x - 0.1e1 / 0.2e1 - cos(0.2e1 * M_PI* t) / 0.4e1, 0.2e1) + 0.50e2 * pow(y - 0.1e1 / 0.2e1 - sin(0.2e1 * M_PI * t) / 0.4e1, 0.2e1));
    
    b[0] += (- (*U)[0].m()) * N.m();
  }
  
};
class DWRMassRhs : public virtual DomainRightHandSide
{
  std::string GetName() const {return "DWRMassRhs";}
  int GetNcomp() const {return 1; }

  mutable FemFunction* U;
  
  void SetFemData(FemData& q) const 
  {
    assert(q.find("U") != q.end() );
    U= &q["U"]; 
  }

  
  void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
  {
    b[0] += (*U)[0].m() * N.m();
  }
};

class DWRNonLin : public virtual DomainRightHandSide
{
  std::string GetName() const {return "DWRMassRhs";}
  int GetNcomp() const {return 1; }

  mutable FemFunction* U;
  
  void SetFemData(FemData& q) const 
  {
    assert(q.find("U") != q.end() );
    U= &q["U"]; 
  }

  
  void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
  {

    b[0] += ( (*U)[0].x() * (*U)[0].x() + (*U)[0].y() *(*U)[0].y())*N.m();
  }
};



class DWRNonLinMatrix : public virtual DomainRightHandSide
{
  std::string GetName() const {return "DWRMassRhs";}
  int GetNcomp() const {return 1; }


  mutable FemFunction* Pu;
  mutable FemFunction* W;
  
  
  void SetFemData(FemData& q) const 
  {
    
    assert(q.find("Pu") != q.end() );
    Pu= &q["Pu"];

    assert(q.find("W") != q.end() );
    W= &q["W"]; 
    
  }

  
  void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
  {

    b[0]+=2* ( (*Pu)[0].x()* (*W)[0].x()+ (*Pu)[0].y()*(*W)[0].y())*N.m();
    
  }
};

// (u0-uk[0], z(0)-i_h z(0))
void Loop::EstimateInitial(DoubleVector& eta, const GlobalVector& U, GlobalVector &Z,
			   VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f)
{
  cout << "Fehlerschaetzer: Anfangsfehler"<< endl;
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);
  DWRInitialRhs DWRIRHS;
  
  MyS->Zero(f);
    
  MyS->GetGV(u) = U;
  MyS->GetGV(z) = Z;
 

  // Std Disc
  MyS->HNAverage(u);
  MyS->HNAverage(z);
  
  MyS->AddNodeVector("U",u);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), DWRIRHS, -1.0);
  MyS->DeleteNodeVector("U");
  MyS->HNDistribute(f);

  // DWR Disc
  DwrFemQ1Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);

  MyS->AddNodeVector("U",u);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), DWRIRHS, 1.0);
  MyS->DeleteNodeVector("U");
  MyS->HNDistribute(f);


  // Std Disc
  MyS->SetDiscretization(*saveD);

  MyS->SetBoundaryVectorZero(f);

  const GlobalVector& F = MyS->GetGV(f);
  MyS->GetDiscretization()->HNAverage(Z);
  // Auswerten, kein Filtern
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*Z(i,c);
	}  
    }
  cout<<eta.sum()<<"EAT0"<<endl;
}



// EstimateDWRprim(eta1, m, Pu_kM[m],Utotal[m], Utotal[m-1], Ztotal[m], u,oldu,z,f);

void Loop::EstimateDWRprim(DoubleVector& eta, int m, const GlobalVector& Pu_kM, vector<GlobalVector>& U,
			   GlobalVector& Z,VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f)
{
  GetMultiLevelSolver()->SetProblem("LaplaceT");
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);


 
  ////////////////////////// Rechte Seite mit Gauss-Regel; TT ist mitte der Intervalle
  // Std-Diskretisierung

  MyS->Zero(f);
  int I_m = _niter / _M;

  double SAVETIME = TIME;

  for (int l=0;l<I_m; ++l)
    {
      TIME = DT * ((m-1)*I_m + l) + 0.5 * DT -  0.5*sqrt(1.0/3.0) * DT;
      MyS->Rhs(f,-0.5 * DT);
      //  MyS->Rhs(f,-0.5 );
      TIME = DT * ((m-1)*I_m + l) + 0.5 * DT +  0.5*sqrt(1.0/3.0) * DT;
      MyS->Rhs(f,-0.5 * DT);
      //  MyS->Rhs(f,-0.5);
    }
  
  MyS->HNDistribute(f);
  
  // DWR - Disc
  DwrFemQ1Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);
  
  for (int l=0;l<I_m; ++l)
    {
      TIME = DT * ((m-1)*I_m + l) + 0.5 * DT -  0.5*sqrt(1.0/3.0) * DT;
      MyS->Rhs(f, 0.5 * DT);
   
      TIME = DT * ((m-1)*I_m + l) + 0.5 * DT +  0.5*sqrt(1.0/3.0) * DT;
      MyS->Rhs(f, 0.5 * DT);
      
    }

  TIME=SAVETIME;
  MyS->HNDistribute(f);
  //////////////////////////// Form

  MyS->SetDiscretization(*saveD);
  for (int l=1;l<=I_m; l++)
    {
      // std disc
      MyS->GetGV(u)    = U[(m-1)*I_m+l];
      MyS->GetGV(oldu) = U[(m-1)*I_m -1+l];
      MyS->HNAverage(u);
      MyS->HNAverage(oldu);
      MyS->AddNodeVector("oldu", oldu);
      MyS->Form(f,u,1.0);
      MyS->DeleteNodeVector("oldu");
      MyS->HNDistribute(f);

    }
  MyS->SetDiscretization(DWRFEM,true);
  
  for (int l=1;l<=I_m; l++)
    {
      // dwr disc
      MyS->GetGV(u)    = U[(m-1)*I_m+l];
      MyS->GetGV(oldu) = U[(m-1)*I_m -1+l];
      MyS->HNAverage(u);
      MyS->HNAverage(oldu);
      MyS->AddNodeVector("oldu", oldu);
      MyS->Form(f,u,-1.0);
      MyS->DeleteNodeVector("oldu");
      MyS->HNDistribute(f);

    }
  MyS->SetDiscretization(*saveD);
  
  
  
  MyS->SetBoundaryVectorZero(f);
  MyS->Visu("Results/resi_prim",f,m);
  //////////////////////////  Auswerten, kein Filtern
  MyS->GetDiscretization()->HNAverage(Z);
  GlobalVector& F = MyS->GetGV(f);
  assert(eta.size() == Z.n());
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*Z(i,c);
	}
    }
  cout<<eta.sum()<<"ETA1 "<<endl;
}
 

void Loop::EstimateDWRdual(DoubleVector& eta, int m, const GlobalVector& Pu_kM,  GlobalVector& Pu_M,
			   const GlobalVector& OLDZ, GlobalVector& Z,VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f)
  
{

  GetMultiLevelSolver()->SetProblem("dp");
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

  
  ////////////////////////// Rechte Seite mit Gauss-Regel; TT ist mitte der Intervalle
  // Std-Diskretisierung
  MyS->Zero(f);
  MyS->Rhs(f,-1);
  MyS->HNDistribute(f);
 
   
  // DWR - Disc
  DwrFemQ1Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);
  MyS->Rhs(f,1);
  MyS->HNDistribute(f);
  MyS->SetDiscretization(*saveD);
 
  
  //Standart
  MyS->GetGV(z)    =Z;
  MyS->GetGV(u)    = Pu_kM;
  MyS->GetGV(oldu) = OLDZ;
  MyS->HNAverage(z);
  MyS->HNAverage(u);
  MyS->HNAverage(oldu);
  MyS->AddNodeVector("Pu", u);
  MyS->AddNodeVector("oldz", oldu);
  MyS->Form(f,z,1.0);
  MyS->DeleteNodeVector("Pu");
  MyS->DeleteNodeVector("oldz");
  MyS->HNDistribute(f);
 
  
  // dwr disc

  MyS->SetDiscretization(DWRFEM,true);
  MyS->GetGV(z)    =Z;
  MyS->GetGV(u)    = Pu_kM;
  MyS->GetGV(oldu) = OLDZ;
  MyS->HNAverage(u);
  MyS->HNAverage(oldu);
  MyS->HNAverage(z);
  MyS->AddNodeVector("Pu", u);
  MyS->AddNodeVector("oldz", oldu);
  MyS->Form(f,z,-1.0);
  MyS->DeleteNodeVector("Pu");
  MyS->DeleteNodeVector("oldz");
  MyS->HNDistribute(f);
  
  MyS->SetDiscretization(*saveD);
 
  MyS->SetBoundaryVectorZero(f);
  MyS->Visu("Results/resi",f,m);

  //////////////////////////  Auswerteun, kein Filtern
  MyS->GetDiscretization()->HNAverage(Pu_M);
  GlobalVector& F = MyS->GetGV(f);
  
  assert(eta.size() == Pu_M.n());
  for(int i=0; i<Pu_M.n(); i++)
    {
      for (int c=0; c<Pu_M.ncomp(); c++)
	{
	  eta[i] += F(i,c)*(Pu_M(i,c));
	}
      
    }
 
  cout<<eta.sum()<<"ETA11 "<<endl;
}
  

  
void Loop::EstimateNonU(DoubleVector& eta, int m,
		        vector<GlobalVector>& U, GlobalVector& Z,
			VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f)
{
  
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

  DWRNonLin NoLin ;

  MyS->Zero(f);

  double SAVETIME = TIME;
  // Std Disc
  MyS->HNAverage(u);

  int I_m = _niter / _M;
  //mit Box


  for (int l=1;l<=I_m; ++l)
    {
    
      MyS->GetGV(u) = U[(m-1)*I_m + l];
      MyS->AddNodeVector("U",u);
      MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, -DT);
      MyS->DeleteNodeVector("U");
     
    }

  
  // DWR Disc
  DwrFemQ1Q22d  DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);

  // Mit bOXregel 
 
  
  for (int l=1;l<=I_m; ++l)
    {

      MyS->GetGV(u) = U[(m-1)*I_m + l];
      MyS->AddNodeVector("U",u);
      MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, DT);
      MyS->DeleteNodeVector("U");
     
    }

  
 
  MyS->HNDistribute(f);

  TIME=SAVETIME;


  // Std Disc
  MyS->SetDiscretization(*saveD);

  MyS->SetBoundaryVectorZero(f);

  const GlobalVector& F = MyS->GetGV(f);
  MyS->GetDiscretization()->HNAverage(Z);
  
  // Auswerten, kein Filtern
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*Z(i,c);
	}  
    }   
}



void Loop::EstimateNonPu(DoubleVector& eta, 
			 const GlobalVector& U, GlobalVector& Z,
			 VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f)
{
  
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

  DWRNonLin NoLin;

  MyS->Zero(f);
  MyS->GetGV(u) = U;
  
  // Std Disc
  MyS->HNAverage(u);
  
  MyS->AddNodeVector("U",u);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, -DTM);
  MyS->DeleteNodeVector("U");
  MyS->HNDistribute(f);

  // DWR Disc
   DwrFemQ1Q22d   DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);

  MyS->AddNodeVector("U",u);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, DTM);
  MyS->DeleteNodeVector("U");
  MyS->HNDistribute(f);


  // Std Disc
  MyS->SetDiscretization(*saveD);

  MyS->SetBoundaryVectorZero(f);

  const GlobalVector& F = MyS->GetGV(f);
  MyS->GetDiscretization()->HNAverage(Z);
  
  // Auswerten, kein Filtern
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*Z(i,c);
	}  
    }   
}




// (Pmu[m]- Pmu[m-1], z-i_h z [m]) - (u[t_m]- u[t_{m-1}], z-i_h z [m] )




void Loop::EstimateAvg(DoubleVector& eta, GlobalVector& Pu, const GlobalVector &Puold,
		       const GlobalVector& U, const GlobalVector &Uold, GlobalVector& Z,
		       VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f)
{
  cout << "Fehlerschaetzer: Mittelfehler"<< endl;
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);
  
  DWRMassRhs DWRMass;
  MyS->Zero(f);

  ///Ich denke hier ist das vorzeichen falsch 
  /*    
	MyS->GetGV(u) = Pu;
	MyS->GetGV(u).add(-1.0,Puold);
	MyS->GetGV(u).add(-1.0,U);
	MyS->GetGV(u).add( 1.0,Uold);
  */

  MyS->GetGV(u) = U;
  MyS->GetGV(u).add(-1.0,Uold);
  MyS->GetGV(u).add(-1.0,Pu);
  MyS->GetGV(u).add( 1.0,Puold);
 


  // Std Disc
  MyS->HNAverage(u);
  
  MyS->AddNodeVector("U",u);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), DWRMass, -1.0);
  MyS->DeleteNodeVector("U");
  MyS->HNDistribute(f);

  // DWR Disc
  DwrFemQ1Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);

  MyS->AddNodeVector("U",u);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), DWRMass, 1.0);
  MyS->DeleteNodeVector("U");
  MyS->HNDistribute(f);


  // Std Disc
  MyS->SetDiscretization(*saveD);

  MyS->SetBoundaryVectorZero(f);

  const GlobalVector& F = MyS->GetGV(f);
  MyS->GetDiscretization()->HNAverage(Z);
  
  // Auswerten, kein Filtern
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*Z(i,c);
	}  
    }

}









// (Pmu[m]- Pmuk[m], i_h z[m] ) - (Pmu[m-1]- Pmuk[m-1], i_h z[m] )
//     - (u[m]- uk[m], i_h z[m] ) +  (u[m-1]- uk[m-1], i_h z[m] )

void Loop::EstimateRest(DoubleVector& eta, int m,
			const GlobalVector& Pu, const GlobalVector &Puold,
			const GlobalVector& Puk, const GlobalVector &Pukold,
			const GlobalVector& U, const GlobalVector &Uold, GlobalVector& Z,
			VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f)
{
  cout << "Fehlerschaetzer: Restfehler"<< endl;
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);
  int I_m = _niter / _M;
  
  MyS->Zero(f);
  MyS->GetGV(z) = Z;
  MyS->HNAverage(z);
  DWRMassRhs DWRMass;

  //Bleibt gleich 
  GlobalVector DU = Pu;
  DU.add(-1.0,Puold);
  DU.add(-1.0,U);
  DU.add( 1.0,Uold);


    
  MyS->AddNodeVector("U",z);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), DWRMass, -DTM);
  MyS->DeleteNodeVector("U");
      
   


  MyS->HNDistribute(f);
  MyS->GetDiscretization()->HNAverage(Z);
  const GlobalVector& F = MyS->GetGV(f);
  MyS->SetBoundaryVectorZero(f);

 
  // Auswerten, kein Filtern
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*DU(i,c);
	}  
    }   

    

  // DWR Disc
  Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);
  MyS->Zero(f);

  MyS->AddNodeVector("U",z);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), DWRMass, DTM);
  MyS->DeleteNodeVector("U");
      
   


  MyS->HNDistribute(f);
  MyS->SetBoundaryVectorZero(f);

  
  MyS->GetDiscretization()->HNAverage(Z);
  
  // Auswerten, kein Filtern
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*DU(i,c);
	}  
    }   


  /*
    GlobalVector DUk = Puk; DUk *= -1.0; // DUk = -Puk geht nicht
    DUk.add( 1.0,Pukold);
    DUk.add( 1.0,U);
    DUk.add(-1.0,Uold);
  
  
    // Std Disc
    MyS->GetGV(z) = Z;
    MyS->HNAverage(z);
    MyS->AddNodeVector("U",z);
    MyS->Zero(f);
    MyS->GetDiscretization()->Rhs(MyS->GetGV(f), DWRMass, 1.0);
    MyS->DeleteNodeVector("U");
    MyS->HNDistribute(f);
    MyS->SetBoundaryVectorZero(f);



  

    const GlobalVector& F = MyS->GetGV(f);
    // Auswerten, kein Filtern

    for(int i=0; i<Z.n(); i++)
    {
    for (int c=0; c<Z.ncomp(); c++)
    {
    eta[i] += F(i,c)*DUk(i,c);
    }  
    }
  
    // DWR Disc
    DwrFemQ1Q22d DWRFEM;
    DWRFEM.BasicInit(this->_paramfile);
    DiscretizationInterface* saveD = MyS->GetDiscretization();
    MyS->SetDiscretization(DWRFEM,true);

    MyS->AddNodeVector("U",z);
    MyS->Zero(f); // hier muss f auf null gesetzt werden, da wir zweimal auswerten
    MyS->GetDiscretization()->Rhs(MyS->GetGV(f), DWRMass, 1.0);
    MyS->DeleteNodeVector("U");
    MyS->HNDistribute(f);
    MyS->SetBoundaryVectorZero(f);

    // Auswerten, kein Filtern
    for(int i=0; i<Z.n(); i++)
    {
    for (int c=0; c<Z.ncomp(); c++)
    {
    ///Aenderung zu DU_K macht keinen sinn mit Trapezregel...DU
    eta[i] += F(i,c)*DU(i,c);
    }  
    }
  */
  
  // Std Disc
  MyS->SetDiscretization(*saveD);
}

void Loop::EstimateNonMeanPu(DoubleVector& eta, int m,
			     GlobalVector& Pu, vector<GlobalVector>& U, GlobalVector& Z,
			     VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f)
{
  
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);



  DWRNonLinMatrix NoLinM;

 

  MyS->Zero(f);
  MyS->GetGV(u) = Pu;

    
  int I_m = _niter / _M;
  //Intregartion mit Trapezregel von K'(u)  mit Q1

  MyS->GetGV(oldu) = Pu;
  MyS->GetGV(oldu).add(-1.0,U[(m-1)*I_m]);
  MyS->AddNodeVector("Pu",u);
  MyS->AddNodeVector("W",oldu);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, -DT*0.5);
  MyS->DeleteNodeVector("Pu");
  MyS->DeleteNodeVector("W");


	 
  MyS->GetGV(oldu) = Pu;
  MyS->GetGV(oldu).add(-1.0,U[(m)*I_m]);
  MyS->AddNodeVector("Pu",u);
  MyS->AddNodeVector("W",oldu);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, -DT*0.5);
  MyS->DeleteNodeVector("Pu");
  MyS->DeleteNodeVector("W");
	 

 
  for (int l=1;l<I_m; ++l)
    {
      MyS->GetGV(oldu) = Pu;
      MyS->GetGV(oldu).add(-1.0,U[(m-1)*I_m + l]);
      MyS->AddNodeVector("Pu",u);       MyS->AddNodeVector("W",oldu);
      MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, -1*DT);
      MyS->DeleteNodeVector("Pu");
      MyS->DeleteNodeVector("W");
    }
 


  MyS->HNDistribute(f);

  MyS->SetBoundaryVectorZero(f);

  const GlobalVector& F = MyS->GetGV(f);
  MyS->GetDiscretization()->HNAverage(Z);
 
  
  // Auswerten, kein Filtern
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*Z(i,c);
	}  
    }

  // DWR Disc
  Q22d DWRFEM;
  MyS->Zero(f);
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);

  //Integration mi Trapez und Q2
 
  MyS->GetGV(oldu) = Pu;
  MyS->GetGV(oldu).add(-1.0,U[(m-1)*I_m]);
  MyS->AddNodeVector("Pu",u);
  MyS->AddNodeVector("W",oldu);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM,0.5*DT);
  MyS->DeleteNodeVector("Pu");
  MyS->DeleteNodeVector("W");


	 
  MyS->GetGV(oldu) = Pu;
  MyS->GetGV(oldu).add(-1.0,U[(m)*I_m]);
  MyS->AddNodeVector("Pu",u);
  MyS->AddNodeVector("W",oldu);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, 0.5*DT);
  MyS->DeleteNodeVector("Pu");
  MyS->DeleteNodeVector("W");
	 

 
  for (int l=1;l<I_m; ++l)
    {
      MyS->GetGV(oldu) = Pu;
      MyS->GetGV(oldu).add(-1.0,U[(m-1)*I_m + l]);
      MyS->AddNodeVector("Pu",u);       MyS->AddNodeVector("W",oldu);
      MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, 1*DT);
      MyS->DeleteNodeVector("Pu");
      MyS->DeleteNodeVector("W");

    }


  

  MyS->HNDistribute(f);

 
  MyS->SetBoundaryVectorZero(f);
  MyS->GetDiscretization()->HNAverage(Z);
  // Auswerten, kein Filtern
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*Z(i,c);
	}  
    }
  MyS->SetDiscretization(*saveD);
}


 

   


void Loop::EstimateNonMeanU(DoubleVector& eta, int m,
			    GlobalVector& Pu,   vector<GlobalVector>& U, GlobalVector& Z,
			    VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f)
{
  
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

  DWRNonLinMatrix NoLinM;


  cout<<"Ohoh"<<endl;
  MyS->Zero(f);
  int I_m = _niter / _M;

  MyS->GetGV(u) = U[(m-1)*I_m];
  MyS->GetGV(oldu) = Pu;
  MyS->GetGV(oldu).add(-1.0,U[(m-1)*I_m]);
  MyS->AddNodeVector("Pu",u);
  MyS->AddNodeVector("W",oldu);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, -0.5*DT);
  MyS->DeleteNodeVector("Pu");
  MyS->DeleteNodeVector("W");

  MyS->GetGV(u) = U[(m)*I_m];
  MyS->GetGV(oldu) = Pu;
  MyS->GetGV(oldu).add(-1.0,U[(m)*I_m]);
  MyS->AddNodeVector("Pu",u);
  MyS->AddNodeVector("W",oldu);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, -0.5*DT);
  MyS->DeleteNodeVector("Pu");
  MyS->DeleteNodeVector("W");

  
    

  //Intregartion mit Trapezregel von K'(u)  mit Q1
  for (int l=1;l<I_m; ++l)
    {

      MyS->GetGV(u) = U[(m-1)*I_m + l];
      MyS->GetGV(oldu) = Pu;
      MyS->GetGV(oldu).add(-1.0,U[(m-1)*I_m + l]);
      MyS->AddNodeVector("Pu",u);
      MyS->AddNodeVector("W",oldu);
      MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, -1*DT);
      MyS->DeleteNodeVector("Pu");
      MyS->DeleteNodeVector("W");
    }
 
  MyS->HNDistribute(f);
  
  MyS->SetBoundaryVectorZero(f);

  const GlobalVector& F = MyS->GetGV(f);
  MyS->GetDiscretization()->HNAverage(Z);
     
  // Auswerten, kein Filtern
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*Z(i,c);
	}  
    }  

  

  // DWR Disc
  Q22d DWRFEM;
  MyS->Zero(f);
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);

  //Integration mi Trapez und Q2
  MyS->GetGV(u) = U[(m-1)*I_m];
  MyS->GetGV(oldu) = Pu;
  MyS->GetGV(oldu).add(-1.0,U[(m-1)*I_m]);
  MyS->AddNodeVector("Pu",u);
  MyS->AddNodeVector("W",oldu);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, 0.5*DT);
  MyS->DeleteNodeVector("Pu");
  MyS->DeleteNodeVector("W");

  MyS->GetGV(u) = U[(m)*I_m];
  MyS->GetGV(oldu) = Pu;
  MyS->GetGV(oldu).add(-1.0,U[(m)*I_m]);
  MyS->AddNodeVector("Pu",u);
  MyS->AddNodeVector("W",oldu);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, DT*0.5);
  MyS->DeleteNodeVector("Pu");
  MyS->DeleteNodeVector("W");

  //Intregartion mit Trapezregel von K'(u)  mit Q1
  for (int l=1;l<I_m; ++l)
    {

      MyS->GetGV(u) = U[(m-1)*I_m + l];
      MyS->GetGV(oldu) = Pu;
      MyS->GetGV(oldu).add(-1.0,U[(m-1)*I_m + l]);
      MyS->AddNodeVector("Pu",u);
      MyS->AddNodeVector("W",oldu);
      MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, DT);
      MyS->DeleteNodeVector("Pu");
      MyS->DeleteNodeVector("W");
    }

  MyS->HNDistribute(f);

  MyS->SetBoundaryVectorZero(f);

 
   
  MyS->GetDiscretization()->HNAverage(Z);
  MyS->GetDiscretization()->HNAverage(Pu);

  
  // Auswerten, kein Filtern
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*Z(i,c);
	}  
    }
  MyS->SetDiscretization(*saveD);
}


void Loop::EstimateDualError(DoubleVector& eta,
			     DoubleVector& eta0,
			     DoubleVector& eta1,
			     DoubleVector& eta11,
			     DoubleVector& eta2,
			     DoubleVector& eta22,
			     DoubleVector& eta23,
			     DoubleVector& eta3,
			     DoubleVector& eta4,
			     DoubleVector& eta5,
			     vector<GlobalVector>& Utotal,
			     vector<GlobalVector>& Ztotal,
			     vector<GlobalVector>& Pu_kM,
			     vector<GlobalVector>& Pu_M,
			     VectorInterface& u,
			     VectorInterface& oldu,
			     VectorInterface& z,
			     VectorInterface& f
			     )
{
  assert(Pu_kM.size()==_M+1);
  assert(Pu_M.size()==_M+1);
  assert(Ztotal.size()==_M+2);
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

 
  
  eta.zero();
  eta0.zero();  eta1.zero();eta11.zero();  eta2.zero();eta22.zero();eta23.zero(); eta3.zero();  eta4.zero(),eta5.zero();  
  assert(eta.size() == MyS->GetMesh()->nnodes());


  // Teil 0 - Anfangsdaten
  EstimateInitial(eta0,Utotal[0], Ztotal[0], u,oldu,z,f);
  
  for (int m=1;m<=_M;++m)
    {
      cout << "m="<< m << "\\Fehlerschaetzer DWR" << endl;

      ///////////// Teil 1 - DWR-Anteil     
      EstimateDWRprim(eta1, m, Pu_kM[m],Utotal,Ztotal[m], u,oldu,z,f);
      cout<<eta1.sum()<<"ETA1 "<<endl;

      /// Teil 1.1 duales residuum

    
      //  EstimateDWRdual(eta11,m,Pu_kM[m], Pu_kM[m], Ztotal[m+1], Ztotal[m],u, oldu,z,f);

      EstimateDWRdual(eta11,m,Pu_kM[m], Pu_kM[m], Ztotal[m+1], Ztotal[m],u, oldu,z,f);
      
    
      ///////////// Teil 2 - mittelungsanteil
      int I_m = _niter / _M;
      EstimateAvg(eta2, Pu_M[m], Pu_M[m-1],  Utotal[m*I_m], Utotal[(m-1)*I_m], Ztotal[m], u,oldu,z,f);
      cout<<eta2.sum()<<"ETA2 "<<endl;

      ///Teil 2.2 Nichtlinearer Teil in u K(u)(z-ihz)
      EstimateNonU(eta22,m, Utotal, Ztotal[m], u,oldu,z,f);
      cout<<eta22.sum()<<"ETA22 "<<endl;

      ///Teil 2.2 Nichtlinearer Teil in Pu K(Pu)(z-ihz)
      EstimateNonPu(eta23, Pu_kM[m], Ztotal[m], u,oldu,z,f);
      cout<<eta23.sum()<<"ETA23 "<<endl;
    
     
      // Neu
      //////////// Teil 3 - anderer teil
      EstimateRest(eta3, m,Pu_M[m], Pu_M[m-1], Pu_kM[m], Pu_kM[m-1], Utotal[m*I_m], Utotal[(m-1)*I_m], Ztotal[m], u, oldu, z,f);

	 
      //Teil 4 Fehler vom geittelten Problem zu Pu
      EstimateNonMeanPu(eta4, m, Pu_M[m],Utotal, Ztotal[m], u,oldu,z,f);
      cout<<eta4.sum()<<"ETA4"<<endl;

	
      EstimateNonMeanU(eta5, m, Pu_M[m],Utotal, Ztotal[m], u,oldu,z,f);
      cout<<eta5.sum()<<"ETA5"<<endl;
    }

  eta.add(1.0,eta0);
  eta.add(0.5,eta11);
  eta.add(0.5,eta1);
  eta.add(1.0,eta2);
  eta.add(1.0,eta22);
  eta.add(-1.0,eta23);
  eta22.add(-1.0,eta23);
  cout<<eta22.sum()<<"combi"<<endl;
  eta.add(1.0,eta3);
  eta.add(0.5,eta4);
  eta.add(0.5,eta5);
  
}

void Loop::run(const std::string& problemlabel)
{
  double endtime;

  DataFormatHandler DFH;
  DFH.insert("dt",   &DT, 0.);
  DFH.insert("time", &TIME, 0.);
  DFH.insert("endtime", &endtime, 0.);
  DFH.insert("M", &_M, 0);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile, "Loop");
  assert(DT>0.0);
  assert(_M>0);
  
  // Anzahl der n-Schritte
  _niter = endtime / DT;
  // Intervallgroesse makro
  DTM = endtime / _M;
  cout << "N="<< _niter << "\t M=" << _M << "\t dt=" << DT << "\t dtm=" << DTM << endl;
  assert( fabs(endtime - _M*DTM)<1.e-8);
  assert( fabs(endtime - _niter * DT)<1.e-8);
  
  for (int ADAITER=0;ADAITER<10;++ADAITER)
    {
      TIME=0.0;
 
      // vectors for solution and right hand side
      VectorInterface u("u"), f("f"), oldu("oldu"),z("z"),oldz("oldz");

      
      PrintMeshInformation();
      //initialize problem, solver and vectors 
      GetMultiLevelSolver()->ReInit("LaplaceT");
      GetMultiLevelSolver()->ReInitVector(u);
      GetMultiLevelSolver()->ReInitVector(oldu);
      GetMultiLevelSolver()->ReInitVector(f);
      GetMultiLevelSolver()->ReInitVector(z);
      GetMultiLevelSolver()->ReInitVector(oldz);

      StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
    
      GetMultiLevelSolver()->GetSolver()->OutputSettings();

      GetMultiLevelSolver()->SetProblem("LaplaceT");

      // Speichern der primalen Loesung u in ALLEN schritten! und alle Funktionale
      vector<GlobalVector> Utotal;
      nvector<double>      Jtotal;
      InitSolution(u); 
      GetMultiLevelSolver()->Equ(oldu,1.0,u);
      Utotal.push_back(MyS->GetGV(u));
      GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,0);
      
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;


  
      // Funktional
      nvector<double> functionals = Functionals(u,f);
      Jtotal.push_back(functionals[0]);
      
      // primale Probleme loesen
      SolvePrimalProblem(Utotal,Jtotal,u,oldu,f);
      assert(Jtotal.size() == _niter+1);
      assert(Utotal.size() == _niter+1);

      // Mitteln
    
      int I_m= _niter/_M;
     
      assert(_M*I_m == _niter);
      
      vector<GlobalVector> Pu_kM(_M+1,GlobalVector(MyS->GetGV(u).ncomp(), MyS->GetGV(u).n()));
      vector<GlobalVector> Pu_M (_M+1,GlobalVector(MyS->GetGV(u).ncomp(), MyS->GetGV(u).n()));
      Pu_kM[0]=Utotal[0];
      Pu_M[0]=Utotal[0];


      // Intregartion und Mittelung
      for (int m=1;m<=_M;++m)
	{
	   TrapezInt(Pu_kM[m], Utotal, (m-1)*I_m, m*I_m );
	   // BoxInt(Pu_kM[m],  Utotal,   (m-1)*I_m+1, m*I_m );
	  // BoxInt(Pu_M[m],  Utotal,   (m-1)*I_m+1, m*I_m );
	  GetMultiLevelSolver()->GetSolver()->GetGV(oldu)=Pu_kM[m];
	  GetMultiLevelSolver()->GetSolver()->Visu("Results/PU_KM",oldu,m);
	  TrapezInt(Pu_M[m], Utotal, (m-1)*I_m, m*I_m );
	  GetMultiLevelSolver()->GetSolver()->GetGV(f)=Pu_M[m];
	  // GetMultiLevelSolver()->GetSolver()->Visu("Results/PU_M",f,m);
	}
      
      // Integral ueber funktional berechnen mit Boxregel (deswegen ohne ersten eintrag)
      // double J = DT * Jtotal.sum() - DT * Jtotal[0];
      // Integral mit Trapezregel
       double J = DT * Jtotal.sum() - DT/2 * Jtotal[0]-DT/2 * Jtotal[_M+1];
      cout << "Integral ueber J = "<< J << endl;

      int nnodes = MyS->GetMesh()->nnodes();
      DoubleVector eta(nnodes,0.0), eta1(nnodes,0.0),eta11(nnodes,0.0),eta2(nnodes,0.0),eta22(nnodes,0.0),eta23(nnodes,0.0),eta3(nnodes,0.0),eta4(nnodes,0.0),eta5(nnodes,0.0),eta0(nnodes,0.0),zeta(nnodes,0.0);
  
      // Duales Problem loesen. Alle dualen Schritte speichern. 
      vector<GlobalVector> Ztotal(_M+2, GlobalVector(MyS->GetGV(u).ncomp(), MyS->GetGV(u).n()));
      Ztotal[_M+1].zero();
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      for (int m=_M;m>=1;--m)
	{
	  cout << "Duales Problem loesen"<< endl;
	  GetMultiLevelSolver()->GetSolver()->GetGV(oldz)=Ztotal[m+1];
	  GetMultiLevelSolver()->GetSolver()->GetGV(u)=Pu_kM[m];
	  GetMultiLevelSolver()->AddNodeVector("oldz", oldz);
	  GetMultiLevelSolver()->AddNodeVector("Pu", u); 
	  SolveDualProblem(z,f,DTM);
	  GetMultiLevelSolver()->GetSolver()->Visu("Results/z",z,1000*ADAITER+m);
	  
	  GetMultiLevelSolver()->DeleteNodeVector("oldz");
	  GetMultiLevelSolver()->DeleteNodeVector("Pu");
	  Ztotal[m] = MyS->GetGV(z);
	}
  
      Ztotal[0] = Ztotal[1]; 
     
    

      
      EstimateDualError(eta,eta0,eta1,eta11,eta2,eta22,eta23,eta3,eta4,eta5, Utotal, Ztotal, Pu_kM, Pu_M,u,oldu,z,f);

      
            
      //GetMultiLevelSolver()->GetSolver()->Visu("Results/dz",z,ADAITER);
     
      // Fehlerschaetzer auswerten
      //cout <<"ETTA"<< eta.sum() << endl;
      // 0.034519  neue Beispiel//100 // 0.0345224 ein spaeter angefangen.
      // a=0.0429182 ///1000
      // I=fabs(eta4.sum())/fabs(( 0.0429135-J));       //160
      // I=fabs(eta4.sum())/fabs((0.0429079-J)); // 80
      // I=fabs(eta4.sum())/fabs((0.0428962-J)); //40
     
   
      this->EtaVisu("Results/eta",ADAITER,eta);

      stringstream str;
      str << "CHECK.txt";
      ofstream OUTF(str.str().c_str(),ios::app);
      OUTF.precision(10);

      //  OUTF << GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes() << " " << DT << " " << DTM  << " " << J << " " << eta.sum() << " " << eta1.sum() << " " << eta2.sum() << " " << eta3.sum() <<" " << eta4.sum() << endl;


 


      OUTF << GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes() << " " << DTM  << " " << J <<  " " << eta.sum()  <<" " << eta0.sum()<<" " << eta1.sum()<<" " << eta11.sum()<<" " << eta2.sum() <<" " << eta22.sum()<<" "<< eta3.sum()<< " "<< eta4.sum()<< " "<< eta5.sum()<<endl;
      
      
      // Gitter verfeinern

      IntVector refcells, coarsecells;
      assert(eta.size()==GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes());
      double mean = eta.norm_l1()/(1.0*eta.size());
      for (int i=0;i<eta.size();++i)
	//	if (fabs(eta[i])>1.5*mean)
	refcells.push_back(i);
      GetMeshAgent()->refine_nodes(refcells);


      PrintMeshInformation();
     
      // GetMultiLevelSolver()->GetSolver()->Visu("Results/neuu",u,ADAITER);
      //GetMultiLevelSolver()->GetSolver()->Visu("Results/PU",oldu,ADAITER);
    
      
    }
  
 
}






