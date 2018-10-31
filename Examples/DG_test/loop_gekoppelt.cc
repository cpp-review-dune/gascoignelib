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
bool   PRIMALPROBLEM,FIRSTDUAL, LASTDUAL;


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

  TIME = x2;
  GetMultiLevelSolver()->GetSolver()->Rhs(f,wgt * DT);

  TIME = SAVETIME;

  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
  GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(u);
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);

  string status = GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());

  Output(u,name);

  return status;
}

string Loop::SolveTransportSingle(VectorInterface& h, VectorInterface& f, string name)
{
  GetMultiLevelSolver()->SetProblem("Transport");
  PRIMALPROBLEM = false;

  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

  GetMultiLevelSolver()->GetSolver()->Zero(f);
  // Rechte Seite mit Gauss-Regel
  double wgt = 0.5;
  double x1 = TIME-0.5*DT - 0.5*sqrt(1.0/3.0) * DT;
  double x2 = TIME-0.5*DT + 0.5*sqrt(1.0/3.0) * DT;
  double SAVETIME = TIME;
  TIME = x1;
  GetMultiLevelSolver()->GetSolver()->Rhs(f,wgt * DT);

  TIME = x2;
  GetMultiLevelSolver()->GetSolver()->Rhs(f,wgt * DT);

  TIME = SAVETIME;
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
  GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(h);
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(h);
  string status = GetMultiLevelSolver()->Solve(h,f,GetSolverInfos()->GetNLInfo());
  Output(h,name);
 PRIMALPROBLEM = false;  
  return status;
}





void Loop::SolvePrimalProblem(vector<GlobalVector> &Utotal,vector<GlobalVector> &Htotal,  nvector<double>& Jtotal, VectorInterface& u, VectorInterface& oldu,VectorInterface& h, VectorInterface& oldh, VectorInterface& f, int ADAITER)
{
 PRIMALPROBLEM = false;
  for (_iter=1; _iter<=_niter; _iter++)
    {
      TIME += DT;
      
      cout << "\n Zeitschritt " << _iter << " " << TIME-DT << " -> " 
	   << TIME<< "  [" << DT << "]" << endl;
       
      GetMultiLevelSolver()->GetSolver()->Equ(oldu,1.0, u);
      GetMultiLevelSolver()->GetSolver()->Equ(oldh,1.0, h);
     
/*
      // Momentengleichung
      GetMultiLevelSolver()->SetProblem("LaplaceT");
      PRIMALPROBLEM = true;
      
      GetMultiLevelSolver()->AddNodeVector("oldu", oldu);
      string res = SolvePrimalSingle(u,f,"Results/u");
      Utotal.push_back(GetMultiLevelSolver()->GetSolver()->GetGV(u));
       
      GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,_iter+ADAITER*1000);
      GetMultiLevelSolver()->DeleteNodeVector("oldu");

*/
      // Transportproblem
      GetMultiLevelSolver()->SetProblem("Transport");
      //PRIMALPROBLEM = true;
      
      GetMultiLevelSolver()->AddNodeVector("oldh", oldh);
     // GetMultiLevelSolver()->AddNodeVector("V", u);
      //GetMultiLevelSolver()->AddNodeVector("oldV", oldu);
      string resi = SolveTransportSingle(h,f,"Results/h");
      
      Htotal.push_back(GetMultiLevelSolver()->GetSolver()->GetGV(h));
      Utotal.push_back(GetMultiLevelSolver()->GetSolver()->GetGV(u));
      
       // Funktional
      nvector<double> functionals = Functionals(h,f);
      Jtotal.push_back(functionals[0]);

      GetMultiLevelSolver()->GetSolver()->Visu("Results/h",h,_iter+ADAITER*1000);      
      GetMultiLevelSolver()->GetSolver()->Visu("Results/oldh",oldh,_iter+ADAITER*1000);   
      
      GetMultiLevelSolver()->DeleteNodeVector("oldh");
     // GetMultiLevelSolver()->DeleteNodeVector("V");
     // GetMultiLevelSolver()->DeleteNodeVector("oldV");
    }
}

void Loop::SolveDualProblem(vector<GlobalVector>& Ztotal, vector<GlobalVector>& Wtotal,VectorInterface& f, VectorInterface& u,  VectorInterface& oldu, VectorInterface& newu,VectorInterface& h,  VectorInterface& oldh, VectorInterface& newh, VectorInterface& z,  VectorInterface& oldz,VectorInterface& w,  VectorInterface& oldw,const vector<GlobalVector>& Pu_k,vector<GlobalVector>& Htotal,int ADAITER)
{
  for (int m=_M;m>=0;--m)
    {
      FIRSTDUAL = (m==_M); 
      LASTDUAL  = (m==0);
      
      GetMultiLevelSolver()->SetProblem("Transport_Dual");
      
      string resi = SolveDualTransportSingle(Wtotal,Htotal,f, u,oldu,newu,w,oldw,Pu_k,m,"Results/w");
      GetMultiLevelSolver()->GetSolver()->Visu("Results/w",w,1000*ADAITER+m);
      Wtotal[m] = GetMultiLevelSolver()->GetSolver()->GetGV(w);  
    /*
      GetMultiLevelSolver()->SetProblem("dp"); 
      
      PRIMALPROBLEM = false;
       
      if(!FIRSTDUAL)
	{
	  GetMultiLevelSolver()->GetSolver()->GetGV(newu)=Pu_k[m+1];
	  GetMultiLevelSolver()->AddNodeVector("u3",newu);
	  GetMultiLevelSolver()->GetSolver()->GetGV(oldw)=Wtotal[m+1];
	  GetMultiLevelSolver()->AddNodeVector("oldW",oldw);
	  GetMultiLevelSolver()->GetSolver()->GetGV(newh)=Htotal[m+1];
	  GetMultiLevelSolver()->AddNodeVector("newH",newh);
	}
	
      if (!LASTDUAL)
	{
	  GetMultiLevelSolver()->GetSolver()->GetGV(oldu)=Pu_k[m-1];
	  GetMultiLevelSolver()->AddNodeVector("u1",oldu);
	  GetMultiLevelSolver()->GetSolver()->GetGV(w)=Wtotal[m];
	  GetMultiLevelSolver()->AddNodeVector("w",w);
	  GetMultiLevelSolver()->GetSolver()->GetGV(oldh)=Htotal[m-1];
	  GetMultiLevelSolver()->AddNodeVector("oldh",oldh);
	}
      
      GetMultiLevelSolver()->GetSolver()->GetGV(oldz)=Ztotal[m+1];
      GetMultiLevelSolver()->AddNodeVector("oldz",oldz);
      GetMultiLevelSolver()->GetSolver()->GetGV(u)=Pu_k[m];
      GetMultiLevelSolver()->AddNodeVector("u2",u);
      GetMultiLevelSolver()->GetSolver()->GetGV(h)=Htotal[m];
      GetMultiLevelSolver()->AddNodeVector("h",h);
   
      GetMultiLevelSolver()->GetSolver()->Zero(f);
      if (FIRSTDUAL || LASTDUAL)
	GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM/2);
      else
	GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM);
      
      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
      GetMultiLevelSolver()->Zero(z);
      GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(z);
      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(z);
      
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      string status = GetMultiLevelSolver()->Solve(z,f,GetSolverInfos()->GetNLInfo());
      GetMultiLevelSolver()->GetSolver()->Visu("Results/z",z,1000*ADAITER+m);
    
      Ztotal[m] = GetMultiLevelSolver()->GetSolver()->GetGV(z);
      
      
      GetMultiLevelSolver()->DeleteNodeVector("oldz");
      GetMultiLevelSolver()->DeleteNodeVector("u2");
      GetMultiLevelSolver()->DeleteNodeVector("h");
      
      if (!LASTDUAL){
	GetMultiLevelSolver()->DeleteNodeVector("u1");
	GetMultiLevelSolver()->DeleteNodeVector("w");
	GetMultiLevelSolver()->DeleteNodeVector("oldh");
      }
      if (!FIRSTDUAL){
	GetMultiLevelSolver()->DeleteNodeVector("u3");
	GetMultiLevelSolver()->DeleteNodeVector("oldW");   
	GetMultiLevelSolver()->DeleteNodeVector("newH");
      }
      
     */
    }
   
  
  PRIMALPROBLEM = false;

}


string Loop::SolveDualTransportSingle(vector<GlobalVector>& Wtotal,vector<GlobalVector>& Htotal,VectorInterface& f, VectorInterface& u,  VectorInterface& oldu, VectorInterface& newu,VectorInterface& w,VectorInterface& oldw,const vector<GlobalVector>& Pu_k,int m, string name)

{
  GetMultiLevelSolver()->SetProblem("Transport_Dual");
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

  PRIMALPROBLEM = false;

  cout<<m<<"Time_trans_dual"<<endl;
    
  FIRSTDUAL = (m==_M); 
  LASTDUAL  = (m==0);


  GetMultiLevelSolver()->GetSolver()->GetGV(u)=Pu_k[m];
  GetMultiLevelSolver()->AddNodeVector("u2",u);
    
  GetMultiLevelSolver()->GetSolver()->GetGV(oldw)=Wtotal[m+1];
  GetMultiLevelSolver()->AddNodeVector("oldw",oldw);
  
    
  
  if(!FIRSTDUAL)
    {
      GetMultiLevelSolver()->GetSolver()->GetGV(newu)=Pu_k[m+1];
      GetMultiLevelSolver()->AddNodeVector("u3",newu);
    }
    
  if (!LASTDUAL)
    {
      GetMultiLevelSolver()->GetSolver()->GetGV(oldu)=Pu_k[m-1];  
      GetMultiLevelSolver()->AddNodeVector("u1",oldu); 
    
    }
    
   
  GetMultiLevelSolver()->GetSolver()->Zero(f);

  if (FIRSTDUAL||LASTDUAL)
    {
      GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM/2);
    }
  else
    {
      GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM);
    }
    
    
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
    
  GetMultiLevelSolver()->Zero(w);
  GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(w);
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(w);

  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  string status = GetMultiLevelSolver()->Solve(w,f,GetSolverInfos()->GetNLInfo());

  
  GetMultiLevelSolver()->DeleteNodeVector("oldw");
  GetMultiLevelSolver()->DeleteNodeVector("u2");

  
    
  if (!LASTDUAL)
    {
      GetMultiLevelSolver()->DeleteNodeVector("u1"); 
    }
  if (!FIRSTDUAL){
    GetMultiLevelSolver()->DeleteNodeVector("u3"); 
  }
    
  return status;
    

}




void Loop::EstimateDWRprim(DoubleVector& eta, int m, const GlobalVector& Pu_kM, vector<GlobalVector>& U,vector<GlobalVector>& H,
			   GlobalVector& Z,GlobalVector& W,VectorInterface& u, VectorInterface& oldu,VectorInterface& h, VectorInterface& oldh,VectorInterface& z,VectorInterface& f)
{
  GetMultiLevelSolver()->SetProblem("Transport");

  PRIMALPROBLEM = true;

  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);



  ////////////////////////// Rechte Seite mit Gauss-Regel; TT ist mitte der Intervalle
  // Std-Diskretisierung Gerde Transport rhs 0 

  MyS->Zero(f);
   
  int I_m = _niter / _M;
    
  

  for (int l=0;l<I_m; ++l)
    {
      TIME = DT * ((m-1)*I_m + l) + 0.5 * DT -  0.5*sqrt(1.0/3.0) * DT;
      MyS->Rhs(f,-0.5 * DT);
      TIME = DT * ((m-1)*I_m + l) + 0.5 * DT +  0.5*sqrt(1.0/3.0) * DT;
      MyS->Rhs(f,-0.5 * DT);
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
  
   

  MyS->HNDistribute(f);
  //////////////////////////// Form
   

  MyS->SetDiscretization(*saveD);
  for (int l=1;l<=I_m; l++)
    {
      // std disc
     // MyS->GetGV(u)    = U[(m-1)*I_m+l]; 
     // MyS->GetGV(oldu) = U[(m-1)*I_m+l-1];  
      MyS->GetGV(h)    = H[(m-1)*I_m+l];
      MyS->GetGV(oldh) = H[(m-1)*I_m-1+l];
     // MyS->HNAverage(u);
    //  MyS->HNAverage(oldu);
      MyS->HNAverage(oldh);
      MyS->HNAverage(h);
      MyS->AddNodeVector("oldh", oldh);
     // MyS->AddNodeVector("V", u);
     // MyS->AddNodeVector("oldV", oldu);
    
      MyS->Form(f,h,1.0);
   //   MyS->DeleteNodeVector("V");
      MyS->DeleteNodeVector("oldh");
     // MyS->DeleteNodeVector("oldV");
      MyS->HNDistribute(f);

    }
  MyS->SetDiscretization(DWRFEM,true);

  for (int l=1;l<=I_m; l++)
    {
      // dwr disc
   //   MyS->GetGV(u)    = U[(m-1)*I_m+l];
   //   MyS->GetGV(oldu) = U[(m-1)*I_m -1+l];
      MyS->GetGV(h)    = H[(m-1)*I_m+l];
      MyS->GetGV(oldh) = H[ (m-1)*I_m-1+l];
  //    MyS->HNAverage(u);
  //    MyS->HNAverage(oldu);
      MyS->HNAverage(h);
      MyS->AddNodeVector("oldh", oldh);
  //    MyS->AddNodeVector("V", u);
   //   MyS->AddNodeVector("oldV", oldu);
    
      MyS->Form(f,h,-1.0);
      MyS->DeleteNodeVector("oldh");
      MyS->DeleteNodeVector("V");
   //   MyS->DeleteNodeVector("oldV");
      MyS->HNDistribute(f);
    }
  MyS->SetDiscretization(*saveD);



  MyS->SetBoundaryVectorZero(f);
 
  //////////////////////////  Auswerten, kein Filtern
  MyS->GetDiscretization()->HNAverage(Z);
  GlobalVector& F = MyS->GetGV(f);
  assert(eta.size() == W.n());
  MyS->Visu("Results/resi_prim",f,m);
    
  for(int i=0; i<W.n(); i++)
    {
      eta[i] += F(i,0)*W(i,0);
    } 
  cout<<eta.sum()<<"ETA1_T "<<endl;    
}


void Loop::EstimateDWRprimBurger(DoubleVector& eta, int m, const GlobalVector& Pu_kM, vector<GlobalVector>& U,vector<GlobalVector>& H,
				 GlobalVector& Z,GlobalVector& W,VectorInterface& u, VectorInterface& oldu,VectorInterface& h, VectorInterface& oldh,VectorInterface& z,VectorInterface& f)
{
    
    
  GetMultiLevelSolver()->SetProblem("LaplaceT");
  PRIMALPROBLEM = true;

  
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);


  ////////////////////////// Rechte Seite mit Gauss-Regel; TT ist mitte der Intervalle
  // Std-Diskretisierung

  MyS->Zero(f);
  int I_m = _niter / _M;
  for (int l=0;l<I_m; ++l)
    {
      TIME = DT * ((m-1)*I_m + l) + 0.5 * DT -  0.5*sqrt(1.0/3.0) * DT;
      MyS->Rhs(f,-0.5 * DT);
   
      TIME = DT * ((m-1)*I_m + l) + 0.5 * DT +  0.5*sqrt(1.0/3.0) * DT;
      MyS->Rhs(f,-0.5 * DT);
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
  
  MyS->Visu("Results/resi_prim_Burger",f,m);
  
  for(int i=0; i<Z.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  eta[i] += F(i,c)*Z(i,c);
	}
    }
  cout<<eta.sum()<<"ETA1 "<<endl;
}  




 
      

 
void Loop::EstimateDWRdual(DoubleVector& eta, int m, vector<GlobalVector>&Pu_kM,  GlobalVector& Pu_M,
			   const GlobalVector& OLDZ, GlobalVector& Z,vector< GlobalVector>& Htotal,vector<GlobalVector>& Wtotal,VectorInterface& u, VectorInterface& oldu,VectorInterface& newu,VectorInterface& h, VectorInterface& oldh,VectorInterface& z,VectorInterface& oldz,VectorInterface& w,VectorInterface& oldw,VectorInterface& f)
  
{
  PRIMALPROBLEM = false;

  FIRSTDUAL = (m==_M); 
  LASTDUAL  = (m==0);

   
  GetMultiLevelSolver()->SetProblem("Transport_Dual");

  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

  GetMultiLevelSolver()->GetSolver()->Zero(f);

// RHS  
  if (FIRSTDUAL || LASTDUAL)
    GetMultiLevelSolver()->GetSolver()->Rhs(f,-DTM/2.0);
  else
    GetMultiLevelSolver()->GetSolver()->Rhs(f,-DTM);

  // FORM
  GetMultiLevelSolver()->GetSolver()->GetGV(oldw)=Wtotal[m+1];
  GetMultiLevelSolver()->AddNodeVector("oldw",oldw);
  MyS->GetGV(w) =Wtotal[m];
  MyS->Form(f,w,1.0);
  GetMultiLevelSolver()->DeleteNodeVector("oldw"); 
  MyS->HNDistribute(f);



  // DWR - Disc
  DwrFemQ1Q22d DWRFEM1;
  DWRFEM1.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM1,true);

  // RHS
  if (FIRSTDUAL || LASTDUAL)
    GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM/2.);
  else
    GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM);

  // FORM
  GetMultiLevelSolver()->GetSolver()->GetGV(oldw)=Wtotal[m+1];
  GetMultiLevelSolver()->AddNodeVector("oldw",oldw);
  MyS->GetGV(w) =Wtotal[m];
  MyS->Form(f,w,-1.0);
  GetMultiLevelSolver()->DeleteNodeVector("oldw"); 
  MyS->HNDistribute(f);


  MyS->SetDiscretization(*saveD);
 
  
  MyS->SetBoundaryVectorZero(f);
  MyS->Visu("Results/resi",f,m);

  GlobalVector& F1 = MyS->GetGV(f);
  assert(eta.size() == Htotal[m].n());
  for(int i=0; i<Htotal[m].n(); i++)
    {
      eta[i] += F1(i,0)*(Htotal[m](i,0));
    }
      
   
  cout<<eta.sum()<<"ETA11 "<<endl;

  
  PRIMALPROBLEM = true;
}



void Loop::EstimateDWRdualBurger(DoubleVector& eta, int m, vector<GlobalVector>&Pu_kM,  GlobalVector& Pu_M,
				 const GlobalVector& OLDZ, GlobalVector& Z, vector<GlobalVector>& Htotal,vector<GlobalVector>& Wtotal,VectorInterface& u, VectorInterface& oldu,VectorInterface& newu,VectorInterface& h, VectorInterface& oldh,VectorInterface& newh,VectorInterface& z,VectorInterface& oldz,VectorInterface& w,VectorInterface& oldw,VectorInterface& f)
  
{
  PRIMALPROBLEM = false;
   
  FIRSTDUAL = (m==_M); 
  LASTDUAL  = (m==0);

  GetMultiLevelSolver()->SetProblem("dp");
  
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

  GetMultiLevelSolver()->GetSolver()->Zero(f);
  if (FIRSTDUAL || LASTDUAL)
    GetMultiLevelSolver()->GetSolver()->Rhs(f,-DTM/2);
  else
    GetMultiLevelSolver()->GetSolver()->Rhs(f,-DTM);
  ////////////////////////// Rechte Seite mit Gauss-Regel; TT ist mitte der Intervalle
  // Std-Diskretisierung
 
  MyS->HNDistribute(f);
   
  // DWR - Disc
  DwrFemQ1Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);
  if (FIRSTDUAL || LASTDUAL)
    GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM/2);
  else
    GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM);

  MyS->HNDistribute(f);
  MyS->SetDiscretization(*saveD);
 
  
  //Standart
  MyS->GetGV(z)    =Z;
  MyS->GetGV(u)    = Pu_kM[m];
  MyS->GetGV(oldz) = OLDZ;
  MyS->GetGV(h)    = Htotal[m];

  if(!FIRSTDUAL)
    {
      GetMultiLevelSolver()->GetSolver()->GetGV(newu)=Pu_kM[m+1];
      GetMultiLevelSolver()->AddNodeVector("u3",newu);
      MyS->HNAverage(newu);
      
        
      GetMultiLevelSolver()->GetSolver()->GetGV(oldw)=Wtotal[m+1];
      GetMultiLevelSolver()->AddNodeVector("oldW",oldw);
    
      GetMultiLevelSolver()->GetSolver()->GetGV(newh)=Htotal[m+1];
      GetMultiLevelSolver()->AddNodeVector("newH",newh);
      
    }
  if (!LASTDUAL)
    {
      GetMultiLevelSolver()->GetSolver()->GetGV(oldu)=Pu_kM[m-1];
      GetMultiLevelSolver()->AddNodeVector("u1",oldu);
      MyS->HNAverage(oldu);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(w)=Wtotal[m];
      GetMultiLevelSolver()->AddNodeVector("w",w);
  
      GetMultiLevelSolver()->GetSolver()->GetGV(oldh)=Htotal[m-1];
      GetMultiLevelSolver()->AddNodeVector("oldh",oldh);
    
      
   
      
    }
  
  MyS->HNAverage(z);
  MyS->HNAverage(u);
  MyS->HNAverage(oldz);
  MyS->AddNodeVector("u2", u);
  MyS->AddNodeVector("oldz", oldz);
  MyS->AddNodeVector("h", h);
  MyS->Form(f,z,1.0);
  MyS->DeleteNodeVector("u2");
  MyS->DeleteNodeVector("oldz");
  MyS->DeleteNodeVector("h");
  
  if(!FIRSTDUAL)
    { GetMultiLevelSolver()->DeleteNodeVector("u3");
      GetMultiLevelSolver()->DeleteNodeVector("oldW");   
      GetMultiLevelSolver()->DeleteNodeVector("newH");
      
    }
  if(!LASTDUAL)
    {
      GetMultiLevelSolver()->DeleteNodeVector("u1");
      GetMultiLevelSolver()->DeleteNodeVector("w");
      GetMultiLevelSolver()->DeleteNodeVector("oldh");
    }
    
    
  cout<<"Standart"<<endl;
  MyS->HNDistribute(f);
 
  
  // dwr disc

  MyS->SetDiscretization(DWRFEM,true);
  MyS->GetGV(z)    =Z;
  MyS->GetGV(u)    = Pu_kM[m];
  MyS->GetGV(oldz) = OLDZ;
  MyS->GetGV(h)    = Htotal[m];
  if(!FIRSTDUAL)
    {
      GetMultiLevelSolver()->GetSolver()->GetGV(newu)=Pu_kM[m+1];
      GetMultiLevelSolver()->AddNodeVector("u3",newu);
      GetMultiLevelSolver()->GetSolver()->GetGV(oldw)=Wtotal[m+1];
      GetMultiLevelSolver()->AddNodeVector("oldW",oldw);
    
      GetMultiLevelSolver()->GetSolver()->GetGV(newh)=Htotal[m+1];
      GetMultiLevelSolver()->AddNodeVector("newH",newh);
      

   
      MyS->HNAverage(newu);
    }
  if (!LASTDUAL)
    {
      GetMultiLevelSolver()->GetSolver()->GetGV(oldu)=Pu_kM[m-1];
      GetMultiLevelSolver()->AddNodeVector("u1",oldu);
      MyS->HNAverage(oldu);
      
      
      GetMultiLevelSolver()->GetSolver()->GetGV(w)=Wtotal[m];
      GetMultiLevelSolver()->AddNodeVector("w",w);
  
      GetMultiLevelSolver()->GetSolver()->GetGV(oldh)=Htotal[m-1];
      GetMultiLevelSolver()->AddNodeVector("oldh",oldh);
    
      
    }
  
  MyS->HNAverage(z);
  MyS->HNAverage(u);
  MyS->HNAverage(oldz);
  MyS->AddNodeVector("h", h);
  MyS->AddNodeVector("u2", u);
  MyS->AddNodeVector("oldz", oldz);
  MyS->Form(f,z,-1.0);
  MyS->DeleteNodeVector("u2");
  MyS->DeleteNodeVector("h");
  MyS->DeleteNodeVector("oldz");
  if(!FIRSTDUAL)
    {
      GetMultiLevelSolver()->DeleteNodeVector("u3");
      GetMultiLevelSolver()->DeleteNodeVector("oldW");   
      GetMultiLevelSolver()->DeleteNodeVector("newH");
      
      
    }
  if(!LASTDUAL)
    {
      GetMultiLevelSolver()->DeleteNodeVector("u1");
      GetMultiLevelSolver()->DeleteNodeVector("u1");
      GetMultiLevelSolver()->DeleteNodeVector("w");
      GetMultiLevelSolver()->DeleteNodeVector("oldh");
      
    }

  
  MyS->HNDistribute(f);
  
  MyS->SetDiscretization(*saveD);
 
  MyS->SetBoundaryVectorZero(f);
  MyS->Visu("Results/resi",f,m);

  //////////////////////////  Auswerteun, kein Filtern
  MyS->GetDiscretization()->HNAverage(Pu_M);
  GlobalVector& F = MyS->GetGV(f);
  MyS->Visu("Results/resi_dual_Burger",f,m);
  assert(eta.size() == Pu_M.n());
  for(int i=0; i<Pu_M.n(); i++)
    {
      for (int c=0; c<Pu_M.ncomp(); c++)
	{
	  eta[i] += F(i,c)*(Pu_M(i,c));
	}
    }
 
  cout<<eta.sum()<<"ETA_dual_Burger "<<endl;
  PRIMALPROBLEM = true;
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
			     vector<GlobalVector>& Htotal,
			     vector<GlobalVector>& Wtotal,
			     vector<GlobalVector>& Pu_kM,
			     vector<GlobalVector>& Pu_M,
			     VectorInterface& u,
			     VectorInterface& oldu,
			     VectorInterface& newu,
			     VectorInterface& h,
			     VectorInterface& oldh,
			     VectorInterface& newh,
			     VectorInterface& z,
			     VectorInterface& oldz,
			     VectorInterface& w,
			     VectorInterface& oldw,
			     VectorInterface& f
			     )
{
  assert(Pu_kM.size()==_M+1);
  assert(Pu_M.size()==_M+1);
  assert(Ztotal.size()==_M+2);
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

 
  
  eta.zero();
  eta0.zero();  eta1.zero(); eta11.zero();  eta2.zero();eta22.zero();eta23.zero(); eta3.zero();  eta4.zero(),eta5.zero();  
  assert(eta.size() == MyS->GetMesh()->nnodes());


  // Teil 0 - Anfangsdaten
  //EstimateInitial(eta0,Utotal[0], Ztotal[0], u,oldu,z,f);
  // erste Intervall dual  
  //EstimateDWRprimBurger(eta22,0,Pu_kM, Pu_kM[0], Ztotal[1], Ztotal[0],u, oldu,newu,z, oldz,f);

  for (int m=1;m<=_M;++m)
    {
      cout << "m="<< m << "\\Fehlerschaetzer DWR" << endl;

      ///////////// Teil 1 - DWR-Anteil     
      EstimateDWRprim(eta1, m, Pu_kM[m],Utotal,Htotal,Ztotal[m],Wtotal[m],u,oldu,h, oldh,z,f);
      
      
      EstimateDWRprimBurger(eta2, m, Pu_kM[m],Utotal,Htotal,Ztotal[m],Wtotal[m],u,oldu,h, oldh,z,f);
  
      cout<<eta1.sum()<<"ETA1 "<<endl;

      /// Teil 1.1 duales residuum
    
      
      EstimateDWRdual(eta11,m,Pu_kM, Pu_kM[m], Ztotal[m+1], Ztotal[m],Htotal,Wtotal, u, oldu,newu,h, oldh ,z, oldz,w, oldw,f);
      EstimateDWRdualBurger(eta3,m,Pu_kM, Pu_kM[m], Ztotal[m+1], Ztotal[m],Htotal,Wtotal,u, oldu,newu,h, oldh ,newh,z, oldz,w,oldw,f);

    }
  

  eta.add(0.5,eta11);

  eta.add(0.5,eta1);
  
 // eta.add(0.5,eta2);

 // eta.add(0.5,eta3);

  
  
  cout<<eta.sum()<<"ETA"<<endl;
  
}




void Loop::run(const std::string& problemlabel)
{
  PRIMALPROBLEM = true;
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

  Extrapolator Extra;
  
  for (int ADAITER=0;ADAITER<8;++ADAITER)
    {
      TIME=0.0;
 
      // vectors for solution and right hand side
      VectorInterface u("u"),newu("newu"),f("f"), oldu("oldu"),z("z"),oldz("oldz"),h("h"),oldh("oldh"),newh("newh"),w("w"), oldw("oldw");

      
      PrintMeshInformation();
      //initialize problem, solver and vectors 
      GetMultiLevelSolver()->ReInit("LaplaceT");
      GetMultiLevelSolver()->ReInitVector(u);
      GetMultiLevelSolver()->ReInitVector(oldu);
      GetMultiLevelSolver()->ReInitVector(newu);
      GetMultiLevelSolver()->ReInitVector(f);
      GetMultiLevelSolver()->ReInitVector(z);
      GetMultiLevelSolver()->ReInitVector(oldz);
      
    
      GetMultiLevelSolver()->ReInitVector(h);
      GetMultiLevelSolver()->ReInitVector(oldh);
      GetMultiLevelSolver()->ReInitVector(newh);
    
      GetMultiLevelSolver()->ReInitVector(w);
      GetMultiLevelSolver()->ReInitVector(oldw);

      StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
    
      GetMultiLevelSolver()->GetSolver()->OutputSettings();
    
      // Speichern der primalen Loesung u in ALLEN schritten!
      vector<GlobalVector> Utotal, Htotal;


      GetMultiLevelSolver()->SetProblem("Transport");
      InitSolution(h);
      GetMultiLevelSolver()->Equ(oldh,1.0,h);
      Htotal.push_back(MyS->GetGV(h));
      GetMultiLevelSolver()->GetSolver()->Visu("Results/h",h,0);
      GetMultiLevelSolver()->GetSolver()->Visu("Results/oldh",oldh,0);
    
      // Speichern der primalen Loesung u in ALLEN schritten! und alle Funktionale
   
      nvector<double>      Jtotal;
    
       
  
      // Funktional Transport 
      nvector<double> functionals = Functionals(h,f);
      Jtotal.push_back(functionals[0]);
    
      
      GetMultiLevelSolver()->GetSolver()->OutputSettings();

      GetMultiLevelSolver()->SetProblem("LaplaceT");
      InitSolution(u);
      GetMultiLevelSolver()->Equ(oldu,1.0,u);
      Utotal.push_back(MyS->GetGV(u));
      GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,0);
    
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  

  
      // primale Probleme loesen
      SolvePrimalProblem(Utotal,Htotal,Jtotal,u,oldu,h,oldh,f,ADAITER);
      assert(Jtotal.size() == _niter+1);
      assert(Utotal.size() == _niter+1);

      vector<GlobalVector> Pu_kM(_M+1,GlobalVector(MyS->GetGV(u).ncomp(), MyS->GetGV(u).n()));
      vector<GlobalVector> Pu_M (_M+1,GlobalVector(MyS->GetGV(u).ncomp(), MyS->GetGV(u).n()));
      Pu_M[0]=Utotal[0];
      Pu_kM[0]=Utotal[0];
           
      // Intregartion und Mittelung
      for (int m=1;m<=_M;++m)
	{   
	  Pu_kM[m] = Utotal[m];
	  Pu_M[m]  = Utotal[m];
	}
	
      // Integral mit Trapezregel
      assert(Jtotal.size() == _niter+1);
      nvector<double> J(1);
      J[0] = DT * Jtotal.sum() - DT/2. * Jtotal[0]-DT/2. * Jtotal[_niter];
      Extra.NewValues(J);
      cout << "Integral ueber J = "<< J[0] << endl;
      Extra.Print();
      
      // Duales Problem loesen. Alle dualen Schritte speichern. 
      vector<GlobalVector> Ztotal(_M+2, GlobalVector(MyS->GetGV(u).ncomp(), MyS->GetGV(u).n()));
      Ztotal[_M+1].zero();
      vector<GlobalVector> Wtotal(_M+2, GlobalVector(MyS->GetGV(u).ncomp(), MyS->GetGV(u).n()));
      Wtotal[_M+1].zero();
            
      
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      
      SolveDualProblem(Ztotal,Wtotal,f,u, oldu,newu,h, oldh,newh,z, oldz,w, oldw,Pu_kM,Htotal,ADAITER);

      // Fehlerschatzer
      int nnodes = MyS->GetMesh()->nnodes();
      DoubleVector eta(nnodes,0.0), eta0(nnodes,0.0),eta1(nnodes,0.0),eta11(nnodes,0.0),eta2(nnodes,0.0),eta22(nnodes,0.0),eta23(nnodes,0.0),eta3(nnodes,0.0),eta4(nnodes,0.0),eta5(nnodes,0.0);

      EstimateDualError(eta,eta0,eta1,eta11,eta2,eta22,eta23,eta3,eta4,eta5, Utotal, Ztotal,Htotal,Wtotal, Pu_kM, Pu_M,u,oldu,newu,h,oldh,newh ,z,oldz,w,oldw,f);
     
   
      this->EtaVisu("Results/eta",ADAITER,eta);
      
      stringstream str;
      str << "Split.txt";
      ofstream OUTF(str.str().c_str(),ios::app);
      OUTF.precision(10);

      OUTF << GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes() << " " << DTM  << " " << J <<  " " << eta.sum()  <<" " << eta0.sum()<<" " << eta1.sum()<<" " << eta11.sum()<<" " << eta2.sum() <<" " << eta22.sum()<<" "<< eta3.sum()<< " "<< eta4.sum()<< " "<< eta5.sum()<<endl;
      
      GetMultiLevelSolver()->GetSolver()->Visu("Results/neuu",u,ADAITER);
      // Gitter verfeinern

      IntVector refcells, coarsecells;
      assert(eta.size()==GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes());
      for (int i=0;i<eta.size();++i)
	refcells.push_back(i);
      GetMeshAgent()->refine_nodes(refcells);
   
      PrintMeshInformation();
    }
  
 
}






