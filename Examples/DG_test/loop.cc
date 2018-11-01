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
 


double TIME, DT,DELTAMIN;
bool  Jump;

double   weight_re=0.,  weight_li=1.0;  

class Primal : public virtual DomainRightHandSide
{
  std::string GetName() const {return "DWRMassRhs";}
  int GetNcomp() const {return 2; }
  
  mutable FemFunction* H;

  void SetFemData(FemData& q) const 
  {
    assert(q.find("H") != q.end() );
    H= &q["H"]; 
  }
  
  void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
  {
    b[0] +=  (*H)[0].m() * N.m();
    b[1] +=  (*H)[1].m() * N.m();
  }
};

class Dual : public virtual DomainRightHandSide
{
  std::string GetName() const {return "DWRMassRhs";}
  int GetNcomp() const {return 2; }
  
  mutable FemFunction* Q;
  mutable FemFunction* V;
  double _split,rho;

public:
  
  Dual(){abort();}
  Dual(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("split",   &_split,-1);
     DFH.insert("rho",   &rho,0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
    assert(_split>=0);
    assert(rho>0);
  }
  
  void SetFemData(FemData& q) const 
  {
    assert(q.find("Q") != q.end() );
    Q= &q["Q"]; 
    assert(q.find("V") != q.end() );
    V= &q["V"]; 
  }

  
  void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
  {
    b[0] +=  rho * (*Q)[0].m() * (*V)[0].m() * N.m();
    b[0] +=  rho * (*Q)[1].m() * (*V)[1].m() * N.m();
    b[1] +=  rho * (*Q)[0].m() * (*V)[0].m() * N.m();
    b[1] +=  rho * (*Q)[1].m() * (*V)[1].m() * N.m();
  }
};

class BurgerDual : public virtual DomainRightHandSide
{
  std::string GetName() const {return "DWRMassRhs";}
  int GetNcomp() const {return 2; }
  
  mutable FemFunction* Q;
  mutable FemFunction* H;

  double _split,rho;
public:

  BurgerDual(){abort();}
  BurgerDual(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("split",   &_split,-1);
    DFH.insert("rho",   &rho,0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
    assert(_split>=0);
    assert(rho>0);
  }

  void SetFemData(FemData& q) const 
  {
    assert(q.find("Q") != q.end() );
    Q= &q["Q"]; 
    assert(q.find("H") != q.end() );
    H= &q["H"]; 
  }

  
  void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
  {
    b[0] +=  rho * (*Q)[0].m() * (*H)[0].m() * N.m();
    b[1] +=  rho * (*Q)[1].m() * (*H)[0].m() * N.m();
  }
};




class SplitPrimalBurger : public virtual DomainRightHandSide
{
  std::string GetName() const {return "DWRMassRhs";}
  int GetNcomp() const {return 2; }
  
  mutable FemFunction* H;
  mutable FemFunction* V;

  double _split,rho;

public:
  SplitPrimalBurger(){abort();}
  SplitPrimalBurger(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("split",   &_split,-1);
    DFH.insert("rho",   &rho,0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
    assert(_split>=0);
    assert(rho>0);
  }

  void SetFemData(FemData& q) const 
  {
    assert(q.find("H") != q.end() );
    H= &q["H"]; 
    assert(q.find("V") != q.end() );
    V= &q["V"]; 
  }
  
  void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
  {
    b[0] += (*V)[0].m() * (rho*(*H)[0].m()) * N.m();
    b[1] += (*V)[1].m() * (rho*(*H)[0].m()) * N.m();
  }
};




class Split : public virtual DomainRightHandSide
{
  double _split,rho;
public:
  std::string GetName() const {return "DWRMassRhs";}
  int GetNcomp() const {return 2; }

  Split(){abort();}
  Split(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("split",   &_split,-1);
     DFH.insert("rho",   &rho,0.0);
    FileScanner FS(DFH);
    FS.NoComplain();
    FS.readfile(pf, "Equation");
    assert(_split>=0);
    assert(rho>0);
  }

  mutable FemFunction* V;
  mutable FemFunction* Q;
  
  void SetFemData(FemData& q) const 
  {
    assert(q.find("V") != q.end() );
    V= &q["V"]; 
    assert(q.find("Q") != q.end() );
    Q= &q["Q"]; 
  }

  
  void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
  {
    b[0] += N.m() * rho * (*V)[0].m() * (*Q)[0].m();
    b[0] += N.m() * rho*  (*V)[1].m() * (*Q)[1].m();
  }
};




string Loop::SolvePrimalSingle(VectorInterface& u, VectorInterface& f, string name)
{
    
  Jump=true;  
  GetMultiLevelSolver()->GetSolver()->Zero(f);

  GetMultiLevelSolver()->GetSolver()->Rhs(f, DT);
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
  GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(u);
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
  
  //  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  cout << endl << "*** Primal Transport" << endl;
  string status = GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());

  Output(u,name);

  return status;
}

string Loop::SolvePrimalBurgerSingle(VectorInterface& v, VectorInterface& f, string name)
{
  GetMultiLevelSolver()->GetSolver()->Zero(f);
  GetMultiLevelSolver()->GetSolver()->Rhs(f,DT);

  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
  GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(v);
  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(v);

  //  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  cout << endl << "*** Primal Burger" << endl;
  string status = GetMultiLevelSolver()->Solve(v,f,GetSolverInfos()->GetNLInfo());

  Output(v,name);

  return status;
}




void Loop::SolvePrimalProblem(vector<GlobalVector> &Htotal,vector<GlobalVector> &Vtotal,nvector<double>& Jtotal,
			      VectorInterface& h, VectorInterface& oldh,
			      VectorInterface& v, VectorInterface& oldv, VectorInterface& f, int ADAITER)
{
  TIME = 0;
  for (_iter=1; _iter<=_niter; _iter++)
    {
      TIME += DT;
      
      cout << "\n Zeitschritt " << _iter << " " << TIME-DT << " -> " 
	   << TIME<< "  [" << DT << "]" << endl;
      GetMultiLevelSolver()->GetSolver()->Equ(oldv,1.0, v); 
      GetMultiLevelSolver()->GetSolver()->Equ(oldh,1.0, h); 
      
      // Burgergleichung
      GetMultiLevelSolver()->SetProblem("Burger");
      
      Jump=true;
      
      GetMultiLevelSolver()->AddNodeVector("oldV", oldv);
      GetMultiLevelSolver()->AddNodeVector("oldH", oldh);

      if (_iter==1) GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      string resi = SolvePrimalBurgerSingle(v,f,"Results/v");
      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      Vtotal[_iter]=GetMultiLevelSolver()->GetSolver()->GetGV(v);
      
     // Vtotal.push_back(GetMultiLevelSolver()->GetSolver()->GetGV(v));
       GetMultiLevelSolver()->GetSolver()->Visu("Results/v",v,_iter+ADAITER*1000);
      
      GetMultiLevelSolver()->DeleteNodeVector("oldV");
      GetMultiLevelSolver()->DeleteNodeVector("oldH");
      // Funktional
     // nvector<double> functionals = Functionals(v,f);
      //Jtotal.push_back(functionals[0]);
     
      // Transportgleichung
      GetMultiLevelSolver()->SetProblem("LaplaceT");
      Jump=true;
     
      GetMultiLevelSolver()->AddNodeVector("oldH", oldh);
      GetMultiLevelSolver()->AddNodeVector("V", v);
      if (_iter==1) GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      string res = SolvePrimalSingle(h,f,"Results/h");

      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       Htotal[_iter]=GetMultiLevelSolver()->GetSolver()->GetGV(h);
      //  (if TIME>4.5){
     // Htotal.push_back(GetMultiLevelSolver()->GetSolver()->GetGV(u));
      GetMultiLevelSolver()->GetSolver()->Visu("Results/h",h,_iter+ADAITER*1000);
      GetMultiLevelSolver()->DeleteNodeVector("oldH");
      GetMultiLevelSolver()->DeleteNodeVector("V");

      // Funktional
      nvector<double> functionals = Functionals(h,f);
      Jtotal.push_back(functionals[0]);
    }
}

void Loop::SolveDualProblem(vector<GlobalVector>& Ztotal,vector<GlobalVector>& Qtotal,
			    vector<GlobalVector>& Vtotal,vector<GlobalVector>& Htotal,
			    VectorInterface& f,
			    VectorInterface& z,  VectorInterface& nextz,
			    VectorInterface& q,  VectorInterface& nextq, 
			    VectorInterface& v,VectorInterface& oldv,  VectorInterface& nextv, 
			    VectorInterface& h,VectorInterface& oldh,  VectorInterface& nexth,
			    int ADAITER)
{
  for (int m=_M;m>=1;--m)
    {

      /// Duales Transport
      GetMultiLevelSolver()->SetProblem("dp"); 
      Jump=true;

      GetMultiLevelSolver()->GetSolver()->Zero(f);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(h)=Htotal[m];
      GetMultiLevelSolver()->AddNodeVector("H",h);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(oldh)=Htotal[m-1];
      GetMultiLevelSolver()->AddNodeVector("oldH",oldh);
      
      GetMultiLevelSolver()->GetSolver()->Rhs(f,DT);
      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
      GetMultiLevelSolver()->GetSolver()->Visu("Results/f",f,1000*ADAITER+m);
      GetMultiLevelSolver()->Zero(z);
      GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(z);
      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(z);
      
      //      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      
      GetMultiLevelSolver()->GetSolver()->GetGV(nextv)=Vtotal[m+1];
      GetMultiLevelSolver()->AddNodeVector("nextV",nextv);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(nextq)=Qtotal[m+1];
      GetMultiLevelSolver()->AddNodeVector("nextQ",nextq);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(nextz)=Ztotal[m+1];
      GetMultiLevelSolver()->AddNodeVector("nextZ",nextz);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(v)=Vtotal[m];
      GetMultiLevelSolver()->AddNodeVector("V",v);
      

      cout << endl << "*** Dual Transport" << endl;
      if (m==_M)
	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      string status1 = GetMultiLevelSolver()->Solve(z,f,GetSolverInfos()->GetNLInfo());

      Ztotal[m] = GetMultiLevelSolver()->GetSolver()->GetGV(z);
      
        GetMultiLevelSolver()->GetSolver()->Visu("Results/z",z,1000*ADAITER+m);
      
      GetMultiLevelSolver()->DeleteNodeVector("nextZ"); 
      GetMultiLevelSolver()->DeleteNodeVector("V");
      GetMultiLevelSolver()->DeleteNodeVector("nextQ"); 
      GetMultiLevelSolver()->DeleteNodeVector("nextV");
      GetMultiLevelSolver()->DeleteNodeVector("H");
     GetMultiLevelSolver()->DeleteNodeVector("oldH");
      

      ////////////// BURGER DUAL
      GetMultiLevelSolver()->SetProblem("Burger_dual"); 
      
      Jump=true;
      
      GetMultiLevelSolver()->GetSolver()->GetGV(nextq)=Qtotal[m+1];
      GetMultiLevelSolver()->AddNodeVector("nextQ",nextq);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(v)=Vtotal[m];
      GetMultiLevelSolver()->AddNodeVector("V",v);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(h)=Htotal[m];
      GetMultiLevelSolver()->AddNodeVector("H",h);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(oldh)=Htotal[m-1];
      GetMultiLevelSolver()->AddNodeVector("oldH",oldh);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(z)=Ztotal[m];
      GetMultiLevelSolver()->AddNodeVector("Z",z);
    
      GetMultiLevelSolver()->GetSolver()->Zero(f);
      GetMultiLevelSolver()->GetSolver()->Rhs(f,DT);
      
      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
      
      GetMultiLevelSolver()->Zero(q);
      GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(q);
      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(q);
      //      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

      //      GetMultiLevelSolver()->GetSolver()->Visu("Results/qf",f,1000*ADAITER+m);
      cout << endl << "*** Dual Burger" << endl;
      if (m==_M)
	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      string status = GetMultiLevelSolver()->Solve(q,f,GetSolverInfos()->GetNLInfo());
      Qtotal[m] = GetMultiLevelSolver()->GetSolver()->GetGV(q);
    
        GetMultiLevelSolver()->GetSolver()->Visu("Results/q",q,1000*ADAITER+m);
      GetMultiLevelSolver()->DeleteNodeVector("nextQ");
      GetMultiLevelSolver()->DeleteNodeVector("Z");
      GetMultiLevelSolver()->DeleteNodeVector("V");
      GetMultiLevelSolver()->DeleteNodeVector("H");
      GetMultiLevelSolver()->DeleteNodeVector("oldH");

    }
    Ztotal[0] = Ztotal[1];
    Qtotal[0] = Qtotal[1];
    
   GetMultiLevelSolver()->GetSolver()->GetGV(q)=Qtotal[0]; 
    
        GetMultiLevelSolver()->GetSolver()->Visu("Results/q",q,1000*ADAITER); 
}

void Loop::ETAProduct(DoubleVector& eta, const GlobalVector& F, const GlobalVector& Z, double wgt, int MC)
{
  assert(eta.size()==F.n());
  assert(eta.size()==Z.n());
  assert(F.ncomp() == Z.ncomp());
  
  for(int i=0; i<Z.n(); i++)
    for (int c=0;c<MC;++c)
      eta[i] += wgt * F(i,c) * Z(i,c);
}


double Loop::DWRResidual(VectorInterface& f, VectorInterface& u,
			 string U1, VectorInterface& G1, string U2, VectorInterface& G2, string U3, VectorInterface& G3,
			 string U4, VectorInterface& G4, string U5, VectorInterface& G5, string U6, VectorInterface& G6,
			 string U7, VectorInterface& G7, string U8, VectorInterface& G8, string U9, VectorInterface& G9,
			 double rhswgt, double formwgt,
			 DoubleVector& eta, 
			 GlobalVector& Z, GlobalVector& OLDZ,  GlobalVector& NEXTZ,
			 int MC)
{
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);
  
  if (U1!="") MyS->AddNodeVector(U1,G1);if (U2!="") MyS->AddNodeVector(U2,G2);if (U3!="") MyS->AddNodeVector(U3,G3);
  if (U4!="") MyS->AddNodeVector(U4,G4);if (U5!="") MyS->AddNodeVector(U5,G5);if (U6!="") MyS->AddNodeVector(U6,G6);
  if (U7!="") MyS->AddNodeVector(U7,G7);if (U8!="") MyS->AddNodeVector(U8,G8);if (U9!="") MyS->AddNodeVector(U9,G9);
  MyS->Zero(f);
  MyS->Rhs(f,   rhswgt);
  MyS->Form(f,u,formwgt);
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
  ETAProduct(eta, MyS->GetGV(f), Z,1.0, MC);

  
  if (U1!="") MyS->DeleteNodeVector(U1);if (U2!="") MyS->DeleteNodeVector(U2);if (U3!="") MyS->DeleteNodeVector(U3);
  if (U4!="") MyS->DeleteNodeVector(U4);if (U5!="") MyS->DeleteNodeVector(U5);if (U6!="") MyS->DeleteNodeVector(U6);
  if (U7!="") MyS->DeleteNodeVector(U7);if (U8!="") MyS->DeleteNodeVector(U8);if (U9!="") MyS->DeleteNodeVector(U9);
 
 
  DwrFemQ1Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);
  if (U1!="") MyS->AddNodeVector(U1,G1);if (U2!="") MyS->AddNodeVector(U2,G2);if (U3!="") MyS->AddNodeVector(U3,G3);
  if (U4!="") MyS->AddNodeVector(U4,G4);if (U5!="") MyS->AddNodeVector(U5,G5);if (U6!="") MyS->AddNodeVector(U6,G6);
  if (U7!="") MyS->AddNodeVector(U7,G7);if (U8!="") MyS->AddNodeVector(U8,G8);if (U9!="") MyS->AddNodeVector(U9,G9);
  MyS->Zero(f);
  MyS->Rhs(f,   -rhswgt);
  MyS->Form(f,u,-formwgt);
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
   ETAProduct(eta, MyS->GetGV(f), Z,1.0, MC);
  
  /*
  ETAProduct(eta, MyS->GetGV(f), Z, 0.5*( weight_li+weight_re), MC);
  ETAProduct(eta, MyS->GetGV(f),OLDZ,0.5*weight_li, MC);
  ETAProduct(eta, MyS->GetGV(f),NEXTZ,0.5*weight_re, MC);
 
 */
 
  if (U1!="") MyS->DeleteNodeVector(U1);if (U2!="") MyS->DeleteNodeVector(U2);if (U3!="") MyS->DeleteNodeVector(U3);
  if (U4!="") MyS->DeleteNodeVector(U4);if (U5!="") MyS->DeleteNodeVector(U5);if (U6!="") MyS->DeleteNodeVector(U6);
  if (U7!="") MyS->DeleteNodeVector(U7);if (U8!="") MyS->DeleteNodeVector(U8);if (U9!="") MyS->DeleteNodeVector(U9);
  

  
  
MyS->SetDiscretization(*saveD);
  

  return eta.sum();
}	



double Loop::Splitting_Form(VectorInterface& f, VectorInterface& u,
			 string U1, VectorInterface& G1, string U2, VectorInterface& G2, string U3, VectorInterface& G3,
			 string U4, VectorInterface& G4, string U5, VectorInterface& G5, string U6, VectorInterface& G6,
			 string U7, VectorInterface& G7, string U8, VectorInterface& G8, string U9, VectorInterface& G9,
			 double formwgt_H, double formwgt,
			 DoubleVector& eta, 
			 GlobalVector& Z, GlobalVector& OLDZ,  GlobalVector& NEXTZ,
			 int MC)
{
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);
  
  if (U1!="") MyS->AddNodeVector(U1,G1);if (U2!="") MyS->AddNodeVector(U2,G2);if (U3!="") MyS->AddNodeVector(U3,G3);
  if (U4!="") MyS->AddNodeVector(U4,G4);if (U5!="") MyS->AddNodeVector(U5,G5);if (U6!="") MyS->AddNodeVector(U6,G6);
  if (U7!="") MyS->AddNodeVector(U7,G7);if (U8!="") MyS->AddNodeVector(U8,G8);if (U9!="") MyS->AddNodeVector(U9,G9);
  MyS->Zero(f);
  MyS->Form(f,u,formwgt);
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
  ETAProduct(eta, MyS->GetGV(f), Z,1.0, MC);
  
  if (U1!="") MyS->DeleteNodeVector(U1);if (U2!="") MyS->DeleteNodeVector(U2);if (U3!="") MyS->DeleteNodeVector(U3);
  if (U4!="") MyS->DeleteNodeVector(U4);if (U5!="") MyS->DeleteNodeVector(U5);if (U6!="") MyS->DeleteNodeVector(U6);
  if (U7!="") MyS->DeleteNodeVector(U7);if (U8!="") MyS->DeleteNodeVector(U8);if (U9!="") MyS->DeleteNodeVector(U9);
 
 
  DwrFemQ1Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  //MyS->SetDiscretization(DWRFEM,true);
  

  if (U1!="") MyS->AddNodeVector(U1,G1);if (U2!="") MyS->AddNodeVector(U2,G2);if (U3!="") MyS->AddNodeVector(U3,G3);
  if (U4!="") MyS->AddNodeVector(U4,G4);if (U5!="") MyS->AddNodeVector(U5,G5);if (U6!="") MyS->AddNodeVector(U6,G6);
  if (U7!="") MyS->AddNodeVector(U7,G7);if (U8!="") MyS->AddNodeVector(U8,G8);if (U9!="") MyS->AddNodeVector(U9,G9);
  MyS->Zero(f);
  MyS->Form(f,u,formwgt_H);
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
  ETAProduct(eta, MyS->GetGV(f), Z, 0.5*( weight_li+weight_re), MC);
  ETAProduct(eta, MyS->GetGV(f),OLDZ,0.5*weight_li, MC);
  ETAProduct(eta, MyS->GetGV(f),NEXTZ,0.5*weight_re, MC);
  
  if (U1!="") MyS->DeleteNodeVector(U1);if (U2!="") MyS->DeleteNodeVector(U2);if (U3!="") MyS->DeleteNodeVector(U3);
  if (U4!="") MyS->DeleteNodeVector(U4);if (U5!="") MyS->DeleteNodeVector(U5);if (U6!="") MyS->DeleteNodeVector(U6);
  if (U7!="") MyS->DeleteNodeVector(U7);if (U8!="") MyS->DeleteNodeVector(U8);if (U9!="") MyS->DeleteNodeVector(U9);
  

  
MyS->SetDiscretization(*saveD);
  

  return eta.sum();
}		    




double Loop::EstimateDWRprim(DoubleVector& eta, int m,
			     vector<GlobalVector>& H,vector<GlobalVector>& V,vector<GlobalVector>& Z,
			     VectorInterface& h, VectorInterface& oldh,VectorInterface& v, VectorInterface& oldv,VectorInterface& f)
{    
  GetMultiLevelSolver()->SetProblem("LaplaceT");
  
  Jump=true;
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

  ////////////////////////// Rechte Seite mit Gauss-Regel; TT ist mitte der Intervalle
  // Std-Diskretisierung

  MyS->Zero(f);
  
  // std disc
  MyS->GetGV(h)    = H[m];
  MyS->GetGV(oldh) = H[m-1];
  MyS->GetGV(v)    = V[m];
  
  if (m==_M)
  cout << DWRResidual(f,h,
		      "oldH",oldh,
		      "V",v,
		      "",v,
		      "",v, "",v,"",v,
		      "",v, "",v,"",v,
		      -DT,1.0,eta,Z[m], Z[m-1],Z[m], 2) << "\t";
  else  
    cout << DWRResidual(f,h,
		      "oldH",oldh,
		      "V",v,
		      "",v,
		      "",v, "",v,"",v,
		      "",v, "",v,"",v,
			-DT,1.0,eta,Z[m], Z[m-1],Z[m+1], 2) << "\t";
         
/*
  Primal Primal;
  MyS->GetGV(h)      = H[m];
  MyS->GetGV(h).add(-1.0,H[m-1]);
  MyS->HNAverage(h);
  MyS->AddNodeVector("H", h);
  
  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Primal,1.0);
  MyS->DeleteNodeVector("H");
  
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
  ETAProduct(eta, MyS->GetGV(f), Z[m],  1.0, 2);


    DwrFemQ1Q22d DWRFEM;
    DWRFEM.BasicInit(this->_paramfile);
    DiscretizationInterface* saveD = MyS->GetDiscretization();
   // MyS->SetDiscretization(DWRFEM,true);
  
  
  MyS->AddNodeVector("H", h);
  
  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Primal, -1.0);
  MyS->DeleteNodeVector("H");
  
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
 
  
// Hohe Ordnung
  ETAProduct(eta, MyS->GetGV(f), Z[m-1],  weight_li, 2);
 ETAProduct(eta, MyS->GetGV(f), Z[m],  weight_re, 2);
  
  
   MyS->SetDiscretization(*saveD);
  */ 

  return eta.sum();
}  

double Loop::EstimateDWRdual(DoubleVector& eta, int m,
			     vector<GlobalVector>& Htotal, vector<GlobalVector>& Vtotal, vector<GlobalVector>& Ztotal, vector<GlobalVector>& Qtotal,
			     VectorInterface& z, VectorInterface& nextz,
			     VectorInterface& v,  VectorInterface& nextv, 
			     VectorInterface& q, VectorInterface& nextq,
                 VectorInterface& h, VectorInterface& oldh,
			     VectorInterface& f)
{
  GetMultiLevelSolver()->SetProblem("dp");
  Jump=true;
  
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

  GetMultiLevelSolver()->GetSolver()->Zero(f);
  
  //RHS
  double rhswgt = -DT;
  
  MyS->GetGV(v)    =Vtotal[m];
  MyS->GetGV(nextv)=Vtotal[m+1];
  MyS->GetGV(z)=    Ztotal[m];
  MyS->GetGV(nextz)=Ztotal[m+1];
  MyS->GetGV(nextq)=Qtotal[m+1];
    MyS->GetGV(q)=Htotal[m];
    MyS->GetGV(oldh)=Htotal[m-1];
   
  
  if(m==_M)
    cout << DWRResidual(f,z,
			"V",v,
			"Z",z,
			"nextV",nextv,
			"nextZ",nextz,
			"nextQ",nextq,
			"H",q, "oldH",oldh,"",z,"",z,
			rhswgt,1.0,
			eta,Htotal[m],Htotal[m-1],Htotal[m],2) << "\t";
      
  else
    cout << DWRResidual(f,z,
			"V",v,
			"Z",z,
			"nextV",nextv,
			"nextZ",nextz,
			"nextQ",nextq,
			"H",q, "oldH",oldh,"",z,"",z,
			rhswgt,1.0,
			eta,Htotal[m],Htotal[m-1],Htotal[m+1],2) << "\t";
 /*
  Primal Primal; 
  MyS->GetGV(z)      = Ztotal[m];
  MyS->HNAverage(z);
  MyS->AddNodeVector("H", z);
  
  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Primal, 1.0);
     
  MyS->DeleteNodeVector("H");

  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
  ETAProduct(eta, MyS->GetGV(f), Htotal[m],  1.0, 2);
  ETAProduct(eta, MyS->GetGV(f), Htotal[m-1],  -1.0, 2); 
  
 
  DwrFemQ1Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
 // MyS->SetDiscretization(DWRFEM,true);


 
  MyS->AddNodeVector("H", z);
  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Primal, -1.0);
  MyS->DeleteNodeVector("H");
  
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
 
 
  
  
  // Hohe Ordnung

  if (m==_M)
    {
      ETAProduct(eta, MyS->GetGV(f), Htotal[m],   weight_re , 2);
      ETAProduct(eta, MyS->GetGV(f), Htotal[m],   weight_li-weight_re , 2);
      ETAProduct(eta, MyS->GetGV(f), Htotal[m-1],-weight_li , 2);
    }
  else
    {
      ETAProduct(eta, MyS->GetGV(f), Htotal[m+1], weight_re , 2);
      ETAProduct(eta, MyS->GetGV(f), Htotal[m],   weight_li-weight_re , 2);
      ETAProduct(eta, MyS->GetGV(f), Htotal[m-1],-weight_li , 2);
    }

   MyS->SetDiscretization(*saveD);
  
 
 
  Dual Dual(_paramfile);
  MyS->GetGV(v)        = Vtotal[m];
  MyS->GetGV(v).add(-1.0,Vtotal[m-1]);
  MyS->HNAverage(v);
  MyS->AddNodeVector("V", v);

  MyS->GetGV(q)      = Qtotal[m];
  MyS->HNAverage(q);
  MyS->AddNodeVector("Q", q);

  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Dual, 1.0);
  MyS->DeleteNodeVector("V");
  MyS->DeleteNodeVector("Q"); 
  
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
  ETAProduct(eta, MyS->GetGV(f), Htotal[m-1],  1.0, 2);
  
  

  //MyS->SetDiscretization(DWRFEM,true);
  MyS->AddNodeVector("V", v);
  MyS->AddNodeVector("Q", q);
  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Dual, -1.0);
  MyS->DeleteNodeVector("V");
  MyS->DeleteNodeVector("Q"); 
  
//   ETAProduct(eta, MyS->GetGV(f), Htotal[m-1],  1.0, 2);
   
   //Hohhe Ordnung
  
 ETAProduct(eta, MyS->GetGV(f), Htotal[m],  weight_re, 2);
 ETAProduct(eta, MyS->GetGV(f), Htotal[m-1],  weight_li, 2);

  MyS->SetDiscretization(*saveD);
*/
  
  return eta.sum();      
}

double Loop::EstimateDWRprimBurger(DoubleVector& eta, int m,
				   vector<GlobalVector>& V,vector<GlobalVector>& H, vector<GlobalVector>& Q,
				   VectorInterface& v, VectorInterface& oldv,
				   VectorInterface& h, VectorInterface& oldh,VectorInterface& f)
{
  GetMultiLevelSolver()->SetProblem("Burger");
  
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

  MyS->Zero(f);

  TIME = DT * m;
  MyS->GetGV(oldv) = V[m-1];
  MyS->GetGV(v)    = V[m];
  MyS->GetGV(oldh) = H[m-1];
    MyS->GetGV(h) = H[m];
  
  Jump=true;

  if(m==_M)
    cout << DWRResidual(f,
			v, "oldV",oldv,
			"oldH", oldh,
			"",h,
			"",h,"",h,"",h,
			"",h,"",h,"",h,
			-DT, 1.0,
			eta,Q[m],Q[m-1],Q[m],2) << "\t";
  else
    cout << DWRResidual(f,
			v, "oldV",oldv,
			"oldH",oldh,
			"",h,
			"",h,"",h,"",h,
			"",h,"",h,"",h,
			-DT, 1.0,
			eta,Q[m],Q[m-1],Q[m+1],2) << "\t";
 /*
  SplitPrimalBurger Split(_paramfile);
  MyS->GetGV(v)      = V[m];
  MyS->GetGV(v).add(-1.0,V[m-1]);
  MyS->GetGV(oldh)    = H[m-1];

  MyS->HNAverage(v);
  MyS->HNAverage(oldh);
  MyS->AddNodeVector("V", v);
  MyS->AddNodeVector("H", oldh);
  
  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Split, 1.0);
  MyS->DeleteNodeVector("H");
  MyS->DeleteNodeVector("V");
  
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);

  ETAProduct(eta, MyS->GetGV(f), Q[m],  1.0, 2);   // - richtig*


  DwrFemQ1Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
//  MyS->SetDiscretization(DWRFEM,true);
  
  
  MyS->AddNodeVector("V", v);
  MyS->AddNodeVector("H", oldh);
  
  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Split, -1.0);
  MyS->DeleteNodeVector("H");
  MyS->DeleteNodeVector("V");
  
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
  
//Hohe Ordnung
  ETAProduct(eta, MyS->GetGV(f), Q[m-1],  weight_li, 2);        
  ETAProduct(eta, MyS->GetGV(f), Q[m],  weight_re, 2);   
    
MyS->SetDiscretization(*saveD);   
    */

  return eta.sum();
}


double Loop::EstimateDWRdualBurger(DoubleVector& eta, int m, 
				   vector<GlobalVector>&Vtotal,
				   vector<GlobalVector>&Htotal,
				   vector<GlobalVector>&Ztotal,
				   vector<GlobalVector>&Qtotal,
				   VectorInterface& v, VectorInterface& oldv,VectorInterface& nextv,
				   VectorInterface& h, VectorInterface& oldh,VectorInterface& nexth,
				   VectorInterface& q,VectorInterface& nextq,
				   VectorInterface& z,VectorInterface& nextz,
				   VectorInterface& f)
  
{

  GetMultiLevelSolver()->SetProblem("Burger_dual");
  
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);

  GetMultiLevelSolver()->GetSolver()->Zero(f);
  Jump=true;
  
  double rhswgt = -DT;
  
  //Standart
  MyS->GetGV(v)    = Vtotal[m];
  MyS->GetGV(h)    = Htotal[m];
  MyS->GetGV(oldh) = Htotal[m-1];
  MyS->GetGV(z)    = Ztotal[m];
  MyS->GetGV(nextq)    = Qtotal[m+1];
  MyS->GetGV(q)    = Qtotal[m];
   
  if (m==_M)
    cout << DWRResidual(f,q, 
			"V",v,
			"H",h,
			"Z",z,
			"nextQ",nextq,
			"oldH",oldh,
			"",q,
			"",q,
			"",q,
			"",q, 
			rhswgt,1.0,
			eta, Vtotal[m],Vtotal[m-1],Vtotal[m],2) << "\t";  
 
    cout << DWRResidual(f,q, 
			"V",v,
			"H",h,
			"Z",z,
			"nextQ",nextq,
			"oldH",oldh,
			"",q,
			"",q,
			"",q,
			"",q, 
			rhswgt,1.0,
			eta, Vtotal[m],Vtotal[m-1],Vtotal[m+1],2) << "\t";
/*
 SplitPrimalBurger Split(_paramfile);
  

  MyS->GetGV(h)      = Htotal[m-1];
  MyS->GetGV(q)      = Qtotal[m];

  MyS->HNAverage(h);
  MyS->HNAverage(q);

  MyS->AddNodeVector("H", h);
  MyS->AddNodeVector("V", q);   // SplitPrimalBurger erwartet V
  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Split, 1.0);
  MyS->DeleteNodeVector("H");
  MyS->DeleteNodeVector("V");
  
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
  
  ETAProduct(eta, MyS->GetGV(f), Vtotal[m],    1.0, 2);
  ETAProduct(eta, MyS->GetGV(f), Vtotal[m-1], -1.0, 2);  
 
  
  DwrFemQ1Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();
  MyS->SetDiscretization(DWRFEM,true);
  

  
  MyS->AddNodeVector("H", h);
  MyS->AddNodeVector("V", q);
  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Split, -1.0);
  MyS->DeleteNodeVector("H");
  MyS->DeleteNodeVector("V");
  
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
  if (m == _M)
    {
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m],   weight_re, 2);
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m],   weight_li-weight_re , 2);
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m-1],-weight_li , 2);
    }
  else
    {
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m+1], weight_re, 2);
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m],   weight_li-weight_re , 2);
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m-1],-weight_li ,2);
    }
    

   MyS->SetDiscretization(*saveD);  
  
  */
  return eta.sum();
}
  
  


void Loop::Splittingerror(DoubleVector& eta, int m,
			  vector<GlobalVector>& V,vector<GlobalVector>& H, vector<GlobalVector>& Q,
                          VectorInterface& v, VectorInterface& oldv,
			  VectorInterface& h, VectorInterface& oldh,
			  VectorInterface& f)
  
{
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);
  GetMultiLevelSolver()->SetProblem("Burger");
  
  Jump=false;
  
  TIME = DT * m;
  
  
  
  BurgerDual Split(_paramfile);
  
  
  MyS->GetGV(v)      = V[m];
  MyS->GetGV(v).add(-1.0,V[m-1]);
  MyS->GetGV(h)    = H[m-1];
  MyS->GetGV(h).add(-1.0,H[m]);
  
  MyS->HNAverage(v);
  MyS->HNAverage(h);
  MyS->AddNodeVector("H", h);
  MyS->AddNodeVector("Q", v);
  
  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Split, 1.0);
  
  MyS->DeleteNodeVector("Q");
  MyS->DeleteNodeVector("H");
  
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
  ETAProduct(eta, MyS->GetGV(f), Q[m],  1, 2);
  
   cout<<"Splitting"<<eta.sum()<< "\t";
  
   
  DwrFemQ1Q22d DWRFEM;
  DWRFEM.BasicInit(this->_paramfile);
  DiscretizationInterface* saveD = MyS->GetDiscretization();

 // MyS->SetDiscretization(DWRFEM,true);
  
  MyS->HNAverage(v);
  MyS->HNAverage(h);
 
  MyS->AddNodeVector("H", h);
  MyS->AddNodeVector("Q", v);

  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Split, 1.0);
  
  MyS->DeleteNodeVector("H");
  MyS->DeleteNodeVector("Q");

  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);
  
  ETAProduct(eta, MyS->GetGV(f), Q[m-1], weight_li, 2);
  ETAProduct(eta, MyS->GetGV(f), Q[m]  , weight_re, 2);
  
  MyS->SetDiscretization(*saveD);
  MyS->Zero(f);

  TIME = DT * m;
  MyS->GetGV(oldv) = V[m-1];
  MyS->GetGV(v)    = V[m];
  MyS->GetGV(oldh) = H[m-1];
    MyS->GetGV(h) = H[m];
 
  Jump=false;
  // Zuerst Bs
  // niederige Ordnung plus hocher 
  
  
  
  if(m==_M)
    cout << Splitting_Form(f,
			v, "oldV",oldv,
			"oldH", oldh,
			"",h,
			"",h,"",h,"",h,
			"",h,"",h,"",h,
			1.0, 1.0,
			eta,Q[m],Q[m-1],Q[m],2) << "\t";
  else
    cout << Splitting_Form(f,
			v, "oldV",oldv,
			"oldH",oldh,
			"",h,
			"",h,"",h,"",h,
			"",h,"",h,"",h,
			1.0, 1.0,
			eta,Q[m],Q[m-1],Q[m+1],2) << "\t";
 
 if(m==_M)
    cout << Splitting_Form(f,
			v, "oldV",oldv,
			"oldH", h,
			"",h,
			"",h,"",h,"",h,
			"",h,"",h,"",h,
			-1.0, -1.0,
			eta,Q[m],Q[m-1],Q[m],2) << "\t";
  else
    cout << Splitting_Form(f,
			v, "oldV",oldv,
			"oldH",h,
			"",h,
			"",h,"",h,"",h,
			"",h,"",h,"",h,
			-1.0, -1.0,
			eta,Q[m],Q[m-1],Q[m+1],2) << "\t";
  
  
}


void Loop::SplittingErrorDual(DoubleVector& eta, int m,
			      vector<GlobalVector>& V,vector<GlobalVector>& H,
			      vector<GlobalVector>& Q,vector<GlobalVector>& Z,
			      VectorInterface& q, VectorInterface& nextq,
                  VectorInterface& z, VectorInterface& nextz,
			      VectorInterface& v, VectorInterface& nextv,
			      VectorInterface& h,
                  VectorInterface& oldh,VectorInterface& f)
{
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);
  GetMultiLevelSolver()->SetProblem("dp");
  
  Jump=false;
  
  Split Split(_paramfile);

  GetMultiLevelSolver()->GetSolver()->GetGV(v)=V[m];
  GetMultiLevelSolver()->GetSolver()->GetGV(v).add(-1.0,V[m-1]);
  GetMultiLevelSolver()->AddNodeVector("V",v);
  
  GetMultiLevelSolver()->GetSolver()->GetGV(q)=Q[m];
  GetMultiLevelSolver()->AddNodeVector("Q",q);
  
  MyS->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Split, -1.0);
  GetMultiLevelSolver()->DeleteNodeVector("V");
  GetMultiLevelSolver()->DeleteNodeVector("Q");
  
  MyS->SetBoundaryVectorZero(f);
  MyS->HNDistribute(f);

  ETAProduct(eta, MyS->GetGV(f), H[m],   -1, 1);
  ETAProduct(eta, MyS->GetGV(f), H[m-1],  1, 1);
  
  
  
   Jump=false;
  

  GetMultiLevelSolver()->GetSolver()->Zero(f);
  
  //RHS
  double rhswgt = -DT;
  
  MyS->GetGV(v)    =V[m];
  MyS->GetGV(nextv)=V[m+1];
  MyS->GetGV(z)=    Z[m];
  MyS->GetGV(nextz)=Z[m+1];
  MyS->GetGV(nextq)=Q[m+1];
    MyS->GetGV(oldh)=H[m-1];
    MyS->GetGV(h)=H[m];
  
  
  if(m==_M)
    cout << Splitting_Form(f,z,
			"V",v,
			"Z",z,
			"nextV",nextv,
			"nextZ",nextz,
			"nextQ",nextq,
			"",z, "oldH",oldh,"",z,"",z,
			1.0,-1.0,
			eta,H[m],H[m-1],H[m],2) << "\t";
      
  else
    cout << Splitting_Form(f,z,
			"V",v,
			"Z",z,
			"nextV",nextv,
			"nextZ",nextz,
			"nextQ",nextq,
			"",z, "oldH",oldh,"",z,"",z,
			1.0,-1.0,
			eta,H[m],H[m-1],H[m+1],2) << "\t";
            
            
         if(m==_M)
    cout << Splitting_Form(f,z,
			"V",v,
			"Z",z,
			"nextV",nextv,
			"nextZ",nextz,
			"nextQ",nextq,
			"",z, "oldH",h,"",z,"",z,
			-1.0,1.0,
			eta,H[m],H[m-1],H[m],2) << "\t";
      
  else
    cout << Splitting_Form(f,z,
			"V",v,
			"Z",z,
			"nextV",nextv,
			"nextZ",nextz,
			"nextQ",nextq,
			"",z, "oldH",h,"",z,"",z,
			-1.0,1.0,
			eta,H[m],H[m-1],H[m+1],2) << "\t";    
            
          
        
  cout<<eta.sum()<<"\t";
  
  
  
}

 
 
void Loop::SplittingErrorBurgerDual(DoubleVector& eta, int m,
				    vector<GlobalVector>& Vtotal,vector<GlobalVector>& Htotal,
				    vector<GlobalVector>& Qtotal,vector<GlobalVector>& Ztotal,
				    VectorInterface& q, VectorInterface& nextq,VectorInterface& z, VectorInterface& v,  
				    VectorInterface& h,   VectorInterface& oldh, VectorInterface& f)

{

  GetMultiLevelSolver()->SetProblem("Burger_dual");
  
  Jump=false;
  
  BurgerDual Split(_paramfile);
  
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);
  
  GetMultiLevelSolver()->GetSolver()->GetGV(q)= Qtotal[m];
  GetMultiLevelSolver()->AddNodeVector("Q",q);
  GetMultiLevelSolver()->GetSolver()->GetGV(h)= Htotal[m-1];
  MyS->GetGV(h).add(-1.0,Htotal[m]);
  GetMultiLevelSolver()->AddNodeVector("H",h);


  GetMultiLevelSolver()->GetSolver()->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Split, -1.0);

  GetMultiLevelSolver()->DeleteNodeVector("Q");
  GetMultiLevelSolver()->DeleteNodeVector("H");

  
  ETAProduct(eta, MyS->GetGV(f), Vtotal[m],     1.0, 2);
  ETAProduct(eta, MyS->GetGV(f), Vtotal[m-1],-1.0, 2);

    DwrFemQ1Q22d DWRFEM;
    DWRFEM.BasicInit(this->_paramfile);
    DiscretizationInterface* saveD = MyS->GetDiscretization();
   // MyS->SetDiscretization(DWRFEM,true);
    

  GetMultiLevelSolver()->AddNodeVector("Q",q);
  GetMultiLevelSolver()->AddNodeVector("H",h);

  GetMultiLevelSolver()->GetSolver()->Zero(f);
  MyS->GetDiscretization()->Rhs(MyS->GetGV(f), Split, 1.0);

  GetMultiLevelSolver()->DeleteNodeVector("Q");
  GetMultiLevelSolver()->DeleteNodeVector("H");

  if (m==_M)
    {
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m],    weight_li, 2);
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m],    weight_re, 2);

      ETAProduct(eta, MyS->GetGV(f), Vtotal[m-1],-weight_li, 2);
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m],  - weight_re, 2);
    }
  else
    {
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m],   weight_li, 2);
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m+1],  weight_re, 2);

      ETAProduct(eta, MyS->GetGV(f), Vtotal[m-1],- weight_li, 2);
      ETAProduct(eta, MyS->GetGV(f), Vtotal[m],  - weight_re, 2);
    }
    
   MyS->SetDiscretization(*saveD);
      
 
  
     GetMultiLevelSolver()->SetProblem("Burger_dual");
//StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
//  assert(MyS);
  GetMultiLevelSolver()->GetSolver()->Zero(f);
  Jump=false;
 
  double rhswgt = -DT;
  
  //Standart
  MyS->GetGV(v)    = Vtotal[m];
  MyS->GetGV(h)    = Htotal[m];
  MyS->GetGV(oldh) = Htotal[m-1];
  MyS->GetGV(h)=Htotal[m];


  MyS->GetGV(z)    = Ztotal[m];
  MyS->GetGV(nextq)    = Qtotal[m+1];
  MyS->GetGV(q)    = Qtotal[m];
   
  if (m==_M)
    cout << Splitting_Form(f,q, 
			"V",v,
			"H",h,
			"Z",z,
			"nextQ",nextq,
			"oldH",oldh,
			"",q,
			"",q,
			"",q,
			"",q, 
			1.0,-1.0,
			eta, Vtotal[m],Vtotal[m-1],Vtotal[m],2) << "\t";  
  else 
    cout << Splitting_Form(f,q, 
			"V",v,
			"H",h,
			"Z",z,
			"nextQ",nextq,
			"oldH",oldh,
			"",q,
			"",q,
			"",q,
			"",q, 
			1.0,-1.0,
			eta, Vtotal[m],Vtotal[m-1],Vtotal[m+1],2) << "\t";
            
            
            if (m==_M)
    cout << Splitting_Form(f,q, 
			"V",v,
			"H",h,
			"Z",z,
			"nextQ",nextq,
			"oldH",h,
			"",q,
			"",q,
			"",q,
			"",q, 
			-1.0,1.0,
			eta, Vtotal[m],Vtotal[m-1],Vtotal[m],2) << "\t";  
  else 
    cout << Splitting_Form(f,q, 
			"V",v,
			"H",h,
			"Z",z,
			"nextQ",nextq,
			"oldH",h,
			"",q,
			"",q,
			"",q,
			"",q, 
			-1.0,1.0,
			eta, Vtotal[m],Vtotal[m-1],Vtotal[m+1],2) << "\t";

            


 cout<<eta.sum()<< "\t";
  
}


 







void Loop::EstimateDualError(DoubleVector& eta,
			     DoubleVector& eta1,
			     DoubleVector& eta2,
			     DoubleVector& eta3,
			     DoubleVector& eta4,
			     DoubleVector& eta5,
			     DoubleVector& eta6,
			     DoubleVector& eta7,
			     vector<GlobalVector>& Htotal,
			     vector<GlobalVector>& Ztotal,
			     vector<GlobalVector>& Vtotal,
			     vector<GlobalVector>& Qtotal,
			     VectorInterface& h,
			     VectorInterface& oldh,
			     VectorInterface& nexth,
			     
			     VectorInterface& z,VectorInterface& nextz,
			     
			     VectorInterface& v,VectorInterface& oldv,VectorInterface& nextv,
			     
			     VectorInterface& q,VectorInterface& nextq,
			     VectorInterface& f)
{
   
  assert(Ztotal.size()==_M+2);
   
  StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(MyS);
 
  eta.zero(); eta1.zero();   eta2.zero();   eta3.zero();   eta4.zero(); eta5.zero(); eta6.zero(); eta7.zero(); 
  assert(eta.size() == MyS->GetMesh()->nnodes());

  for (int m=1;m<=_M;++m)
    {
      cout << "m="<< m << "\\Fehlerschaetzer DWR" << "\t";
      cout.precision(5);
    
     EstimateDWRprim(eta1, m,
		     Htotal,Vtotal,Ztotal,
		      h, oldh, v, oldv, f);
    
 
 
  
      EstimateDWRdual(eta2, m, Htotal, Vtotal, Ztotal, Qtotal,
		     z, nextz,
		     v, nextv, 
		     q, nextq,
            h, oldh,
		      f);

     
      EstimateDWRprimBurger(eta3, m,
			    Vtotal, Htotal, Qtotal,
			    v, oldv,
			    h, oldh, f);
  
      EstimateDWRdualBurger(eta4, m,
			    Vtotal,Htotal,Ztotal,Qtotal,
			    v, oldv,nextv,
			    h, oldh,nexth,
			    q, nextq,
			    z, nextz,f);
			   
  /*    
    Splittingerror(eta5, m, Vtotal, Htotal,Qtotal,
		    v, oldv,
		     h, oldh,f);
  
  
     SplittingErrorDual(eta6, m, Vtotal,Htotal,Qtotal,Ztotal,
			 q,nextq,
             z,nextz,
            v, nextv,
			 h, oldh,f);
    
     
     
     SplittingErrorBurgerDual(eta7, m,
		       Vtotal,Htotal,Qtotal,Ztotal,
		      q,nextq,z,v,h,oldh,f);
     
  */   
     
      cout << endl;
      
    }
  eta.add(0.5,eta1);
  
  eta.add(0.5,eta2);
  
  eta.add(0.5,eta3);

  eta.add(0.5,eta4);
  
 eta.add(0.5,eta5);
  
  eta.add(0.5,eta6);

  eta.add(0.5,eta7);

  cout<<eta.sum()<<"ETA"<<endl;
  

}

void Loop::run(const std::string& problemlabel)
{
  double endtime;

  DataFormatHandler DFH;
  DFH.insert("dt",   &DT, 0.);
  DFH.insert("time", &TIME, 0.);
  DFH.insert("endtime", &endtime, 0.);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile, "Loop");
  assert(DT>0.0);

  
  
    DFH.insert("deltamin",   &DELTAMIN, 0.);
    FS.NoComplain();
    FS.readfile(_paramfile, "Equation");
    assert(DELTAMIN>0.0);
  
  


  // Anzahl der n-Schritte
  _niter = endtime / DT;
  _M = _niter;
  // Intervallgroesse makro

  assert(_M == _niter);


  cout << "N="<< _niter << "\t M=" << _M << "\t dt=" << DT <<  endl;


  assert( fabs(endtime - _niter * DT)<1.e-8);

  Extrapolator Extra;
  
  for (int ADAITER=0;ADAITER<6;++ADAITER)
    {
      TIME=0.0;
 
      // vectors for solution and right hand side
      VectorInterface
	h("h"),f("f"), oldh("oldh"),nexth("nexth"),
	v("v"), oldv("oldv"), nextv("nextv"),
	z("z"),nextz("nextz"),
	q("q"), nextq("nextq");
      
      PrintMeshInformation();
      //initialize problem, solver and vectors 
      GetMultiLevelSolver()->ReInit("LaplaceT");

      GetMultiLevelSolver()->ReInitVector(h);
      GetMultiLevelSolver()->ReInitVector(oldh);
      GetMultiLevelSolver()->ReInitVector(nexth);

      GetMultiLevelSolver()->ReInitVector(v);
      GetMultiLevelSolver()->ReInitVector(oldv);
      GetMultiLevelSolver()->ReInitVector(nextv);

      GetMultiLevelSolver()->ReInitVector(f);

      GetMultiLevelSolver()->ReInitVector(z);
      GetMultiLevelSolver()->ReInitVector(nextz);
      
      GetMultiLevelSolver()->ReInitVector(q);
      GetMultiLevelSolver()->ReInitVector(nextq);
    
      StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
    
      GetMultiLevelSolver()->GetSolver()->OutputSettings();
    
      // Speichern der primalen Loesung u in ALLEN schritten!
    //  vector<GlobalVector> Htotal,Vtotal,Ztotal,Qtotal;

      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       vector<GlobalVector> Htotal(_niter+2, GlobalVector(MyS->GetGV(h).ncomp(), MyS->GetGV(h).n()));
       Htotal[_niter+1].zero();
      
       vector<GlobalVector> Vtotal(_niter+2, GlobalVector(MyS->GetGV(v).ncomp(), MyS->GetGV(v).n()));
       Vtotal[_niter+1].zero();

       vector<GlobalVector> Ztotal(_M+2, GlobalVector(MyS->GetGV(z).ncomp(), MyS->GetGV(z).n()));
       Ztotal[_M+1].zero();
      
       vector<GlobalVector> Qtotal(_M+2, GlobalVector(MyS->GetGV(q).ncomp(), MyS->GetGV(q).n()));
       Qtotal[_M+1].zero();


      // Speichern der primalen Loesung u in ALLEN schritten! und alle Funktionale
   
      nvector<double>      Jtotal;
    
       
     
      
      GetMultiLevelSolver()->GetSolver()->OutputSettings();

      GetMultiLevelSolver()->SetProblem("LaplaceT");
      InitSolution(h);
       Htotal[0]=MyS->GetGV(h);
         GetMultiLevelSolver()->GetSolver()->Visu("Results/h",h,0);
      GetMultiLevelSolver()->Equ(oldh,1.,h);
      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          
      
      // Funktional Transport 
      nvector<double> functionals = Functionals(h,f);
      Jtotal.push_back(functionals[0]);
    
      ////// Burger PRoblem
      GetMultiLevelSolver()->SetProblem("Burger");
      InitSolution(v);
      GetMultiLevelSolver()->Equ(oldv,1.0,v);
      //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx      
            Vtotal[0]=MyS->GetGV(v);
      
      GetMultiLevelSolver()->GetSolver()->Visu("Results/v",v,0);
     // nvector<double> functionals = Functionals(v,f);
     // Jtotal.push_back(functionals[0]);
     // NICHT!!!!Jtotal.push_back(functionals[0]);

  
      
       cout<<Jtotal.size()<<"vorscheife"<<endl;
      // primale Probleme loesen

      SolvePrimalProblem(Htotal,Vtotal,Jtotal,
			 h,oldh,
			 v,oldv,f,ADAITER);
  
      
      assert(Jtotal.size() == _niter+1);
       assert(Htotal.size() == _niter+2);
       assert(Vtotal.size() == _niter+2);
       assert(Ztotal.size() == _niter+2);
       assert(Qtotal.size() == _niter+2);
       
      // Integral mit Boxregel
      assert(Jtotal.size() == _niter+1);
      nvector<double> J(1);
      J[0] = DT * Jtotal.sum();

      Extra.NewValues(J);
      cout << "Integral ueber J = "<< J[0] << " "<<"10"<<endl;
      Extra.Print();

      // Duales Problem loesen. Alle dualen Schritte speichern. 
      

      
      SolveDualProblem(Ztotal,Qtotal,Vtotal,Htotal,f,
		       z, nextz,
		       q, nextq,
		       v, oldv, nextv,
		       h, oldh, nexth,
		       ADAITER);

      // Fehlerschatzer
      int nnodes = MyS->GetMesh()->nnodes();
      DoubleVector eta(nnodes,0.0), eta1(nnodes,0.0),eta2(nnodes,0.0), eta3(nnodes,0.0),eta4(nnodes,0.0),eta5(nnodes,0.0),eta6(nnodes,0.0),eta7(nnodes,0.0);

      EstimateDualError(eta,eta1,eta2, eta3,eta4,eta5,eta6,eta7,Htotal, Ztotal,Vtotal, Qtotal,
			h,oldh,nexth,
			z,nextz,
			v,oldv,nextv,
			q, nextq,
			f);
     
   
      this->EtaVisu("Results/eta",ADAITER,eta);
      
      stringstream str;
      str << "ice_A_Ort.txt";
      ofstream OUTF(str.str().c_str(),ios::app);
      OUTF.precision(10);

           OUTF << GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes() << " " << DT << " " << J <<  " " << eta.sum()  <<" " << eta1.sum()<<" " << eta2.sum()<<" " << eta3.sum()<<" " << eta4.sum()<<" " << eta5.sum()<<" " << eta6.sum()<<" " << eta7.sum()<<endl;
      //OUTF << GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes() << " " << DT << " " << J <<  " " << eta.sum()  <<" " << eta1.sum() +  eta2.sum() + eta3.sum() +  eta4.sum()<<" " << eta5.sum()<<" " << eta6.sum()<<" " << eta7.sum()<<endl;      

      
      // Gitter verfeinern

      IntVector refcells, coarsecells;
      assert(eta.size()==GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes());
      for (int i=0;i<eta.size();++i)
	  refcells.push_back(i);
      GetMeshAgent()->refine_nodes(refcells);
   
      PrintMeshInformation();
    }
  
 
}






