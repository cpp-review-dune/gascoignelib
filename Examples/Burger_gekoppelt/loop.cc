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



double TIME, DT, CHIGAUSS,DTM1,DTM2;
bool   PRIMALPROBLEM,FIRSTDUAL, LASTDUAL;


string Loop::SolveTransportSingle(VectorInterface& h, VectorInterface& f, string name)
{
GetMultiLevelSolver()->SetProblem("Transport");
PRIMALPROBLEM = false;

GetMultiLevelSolver()->GetSolver()->Zero(f);
GetMultiLevelSolver()->GetSolver()->Rhs(f, DT);
GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(h);
GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(h);
string status = GetMultiLevelSolver()->Solve(h,f,GetSolverInfos()->GetNLInfo());
//  Output(h,name);
PRIMALPROBLEM = true;  
return status;
}



string Loop::SolvePrimalSingle(VectorInterface& u, VectorInterface& f, string name)
{
PRIMALPROBLEM = true;
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


void Loop::SolvePrimalProblem(vector<GlobalVector> &Utotal, VectorInterface& u, VectorInterface& oldu, VectorInterface& f,vector<GlobalVector> &Htotal, VectorInterface& h, VectorInterface& oldh, int ADAITER)
{
PRIMALPROBLEM = true;
for (_iter=1; _iter<=_niter; _iter++)
    {
    TIME += DT;
    
    cout << "\n Zeitschritt " << _iter << " " << TIME-DT << " -> " 
    << TIME<< "  [" << DT << "]" << endl;
    
    GetMultiLevelSolver()->GetSolver()->Equ(oldu,1.0, u);
    GetMultiLevelSolver()->GetSolver()->Equ(oldh,1.0, h);
    
    GetMultiLevelSolver()->SetProblem("Transport");
    GetMultiLevelSolver()->AddNodeVector("oldh", oldh);
    GetMultiLevelSolver()->AddNodeVector("V", oldu);
    
    string resi = SolveTransportSingle(h,f,"Results/h");
    Htotal.push_back(GetMultiLevelSolver()->GetSolver()->GetGV(h));
    GetMultiLevelSolver()->GetSolver()->Visu("Results/h",h,_iter+ADAITER*1000);      
    GetMultiLevelSolver()->DeleteNodeVector("oldh");
    GetMultiLevelSolver()->DeleteNodeVector("V");
    
    
    // Momentengleichung
    GetMultiLevelSolver()->SetProblem("LaplaceT");
    
    GetMultiLevelSolver()->AddNodeVector("oldu", oldu);
    GetMultiLevelSolver()->AddNodeVector("H", h);
    
    string res = SolvePrimalSingle(u,f,"Results/u");
    Utotal.push_back(GetMultiLevelSolver()->GetSolver()->GetGV(u));
    GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,_iter+ADAITER*1000);
    
    GetMultiLevelSolver()->DeleteNodeVector("oldu");
    GetMultiLevelSolver()->DeleteNodeVector("H");
    
    
    }
}

double Loop::CompFunctional(vector<GlobalVector> &Utotal,VectorInterface& u, VectorInterface& f)
{

double Ji=0;
for (int ii=1;ii<=_niter; ii++)

    {
    GetMultiLevelSolver()->GetSolver()->GetGV(u).equ(0.5,Utotal[ii-1], 0.5, Utotal[ii]);
    TIME = (1.0*ii-0.5)*DT;
    
    DoubleVector juh = Functionals(u,f);
    
    Ji += DT*juh[0];

    }
return Ji;


}





void Loop::MittelInt(GlobalVector& avg_old,GlobalVector& avg, const vector<GlobalVector>& U, int start, int stopp, double DTM)
{
assert(start<U.size());  assert(stopp<U.size());
avg.zero();
avg.add(1*DT,U[start]);

for (int l=start+1 ; l<=stopp-1;++l)
    avg.add(2*DT,U[l]);
avg.add(1*DT,U[stopp]);
avg*=1/DTM;
avg.add(-1,avg_old);

}

string Loop::SolveDualSingle(vector<GlobalVector>& Ztotal,vector<GlobalVector>& Wtotal,vector<GlobalVector>& Htotal,VectorInterface& f, VectorInterface& u,  VectorInterface& oldu, VectorInterface& newu, VectorInterface& z,  VectorInterface& oldz,VectorInterface& w,VectorInterface& oldw,VectorInterface& h,VectorInterface& oldh,const vector<GlobalVector>& Pu_k,int m, vector<double>& DT_M,vector<double>& T,string name)

{
GetMultiLevelSolver()->SetProblem("dp");
PRIMALPROBLEM = false;

    DTM1=0.0;
    DTM2=0.0;
    
    DTM1=DT_M[m];

    
    
    FIRSTDUAL = (m==_M); 
    LASTDUAL  = (m==0);


    if(LASTDUAL)
    {TIME=0.0;}
    else
    {TIME=T[m]-DTM1/2;}
    
    if(!FIRSTDUAL)
    {
    GetMultiLevelSolver()->GetSolver()->GetGV(newu)=Pu_k[m+1];
    GetMultiLevelSolver()->AddNodeVector("u3",newu);
    DTM2=DT_M[m+1];
    Wtotal[m+1]=1.0;
    GetMultiLevelSolver()->GetSolver()->GetGV(oldw)=Wtotal[m+1];
    GetMultiLevelSolver()->AddNodeVector("oldW",oldw);
    }
    GetMultiLevelSolver()->GetSolver()->GetGV(u)=Pu_k[m];
    GetMultiLevelSolver()->AddNodeVector("u2",u);
    if (!LASTDUAL)
    {
    Wtotal[m]=1.0; 
    GetMultiLevelSolver()->GetSolver()->GetGV(oldu)=Pu_k[m-1];
    GetMultiLevelSolver()->AddNodeVector("u1",oldu);
    GetMultiLevelSolver()->GetSolver()->GetGV(w)=Wtotal[m];
    GetMultiLevelSolver()->AddNodeVector("w",w);
    }
    
    GetMultiLevelSolver()->GetSolver()->GetGV(oldz)=Ztotal[m+1];
    GetMultiLevelSolver()->AddNodeVector("oldz",oldz);
    
    GetMultiLevelSolver()->GetSolver()->GetGV(h)=Htotal[m];
    GetMultiLevelSolver()->AddNodeVector("h",h);
    GetMultiLevelSolver()->GetSolver()->GetGV(oldh)=Htotal[m+1];
    GetMultiLevelSolver()->AddNodeVector("newH",oldh);
    
    GetMultiLevelSolver()->GetSolver()->Zero(f);

    if (FIRSTDUAL)
    {
    GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM1/2);
    }
    else
    {
    if (LASTDUAL)
    {
        GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM2/2);
    }
    else
    {
        DTM2=DT_M[m+1];  
        GetMultiLevelSolver()->GetSolver()->Rhs(f,(DTM1)/2);
        TIME=T[m+1]-DTM2/2;
        GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM2/2);
    
    }
    }
    GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
    
    GetMultiLevelSolver()->Zero(z);
    GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(z);
    GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(z);

    string status = GetMultiLevelSolver()->Solve(z,f,GetSolverInfos()->GetNLInfo());

    
    GetMultiLevelSolver()->DeleteNodeVector("oldz");
    GetMultiLevelSolver()->DeleteNodeVector("u2");
    GetMultiLevelSolver()->DeleteNodeVector("h");
    GetMultiLevelSolver()->DeleteNodeVector("newH");    
    
    if (!LASTDUAL)
    {
    GetMultiLevelSolver()->DeleteNodeVector("u1");
    GetMultiLevelSolver()->DeleteNodeVector("w");   
    }
    if (!FIRSTDUAL){
    GetMultiLevelSolver()->DeleteNodeVector("u3");
    GetMultiLevelSolver()->DeleteNodeVector("oldW");   
   
    }
    
    DTM1=0.0;
    DTM2=0.0;
    
    return status;
    
}

void Loop::SolveDualProblem(vector<GlobalVector>& Ztotal, vector<GlobalVector>& Wtotal,vector<GlobalVector>& Htotal,VectorInterface& f, VectorInterface& u,  VectorInterface& oldu, VectorInterface& newu, VectorInterface& z,  VectorInterface& oldz,VectorInterface& w,VectorInterface& oldw,VectorInterface& h,VectorInterface& oldh,const vector<GlobalVector>& Pu_k,int ADAITER, vector<double>& DT_M,vector<double>& T)
{


PRIMALPROBLEM = false;

for (int m=_M;m>=0;--m)
    {
        
    
    GetMultiLevelSolver()->SetProblem("dp");

            
    string res = SolveDualSingle(Ztotal,Wtotal,Htotal,f, u,oldu,newu,z,oldz,w,oldw,h,oldh,Pu_k,m, DT_M,T,"Results/z");
    GetMultiLevelSolver()->GetSolver()->Visu("Results/z",z,1000*ADAITER+m);
    Ztotal[m] = GetMultiLevelSolver()->GetSolver()->GetGV(z);   

    }

PRIMALPROBLEM = true;

}






class DWRMassRhs : public virtual DomainRightHandSide
{
std::string GetName() const {return "DWRMassRhs";}
int GetNcomp() const {return 2; }

mutable FemFunction* U;

void SetFemData(FemData& q) const 
{
    assert(q.find("U") != q.end() );
    U= &q["U"]; 
}


void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
{
    b[0] += (*U)[0].m() * N.m();
    b[1]  +=(*U)[1].m() * N.m();
}
};


class DWRNonLin : public virtual DomainRightHandSide
{
std::string GetName() const {return "DWRMassRhs";}
int GetNcomp() const {return 2; }

mutable FemFunction* U;

void SetFemData(FemData& q) const 
{
    assert(q.find("U") != q.end() );
    U= &q["U"]; 
}


void operator()(VectorIterator b, const TestFunction& N, const Vertex2d& v) const 
{

    b[0] += ( (*U)[0].m() * (*U)[0].x() + (*U)[1].m() *(*U)[0].y())*N.m();
    b[1] += ( (*U)[0].m() * (*U)[1].x() + (*U)[1].m() *(*U)[1].y())*N.m();
}
};


class DWRNonLinMatrix : public virtual DomainRightHandSide
{
std::string GetName() const {return "DWRMassRhs";}
int GetNcomp() const {return 2; }


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

    b[0]+=( (*Pu)[0].m()* (*W)[0].x()+ (*Pu)[1].m()*(*W)[0].y())*N.m();
    b[1]+=( (*Pu)[0].m()* (*W)[1].x()+ (*Pu)[1].m()*(*W)[1].y())*N.m();
}
};


void Loop::EstimateDWRprim(DoubleVector& eta, int m, const GlobalVector& Pu_kM, vector<GlobalVector>& U,
            GlobalVector& Z,VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f,vector<double>&T)
{
GetMultiLevelSolver()->SetProblem("LaplaceT");
PRIMALPROBLEM = true;


StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
assert(MyS);



////////////////////////// Rechte Seite mit Gauss-Regel; TT ist mitte der Intervalle
// Std-Diskretisierung

MyS->Zero(f);
int I_m = (T[m]-T[m-1])/DT+1.e-10;

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
int start=T[m-1]/DT+1.e-10;

MyS->SetDiscretization(*saveD);
for (int l=1;l<=I_m; l++)
    {
    // std disc
    MyS->GetGV(u)    = U[start+l];
    MyS->GetGV(oldu) = U[start-1+l];
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
    MyS->GetGV(u)    = U[start+l];
    MyS->GetGV(oldu) = U[start -1+l];
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

void Loop::EstimateDWRdual(DoubleVector& eta, int m, vector<GlobalVector>&Pu_kM,  GlobalVector& Pu_M,
            const GlobalVector& OLDZ, GlobalVector& Z,VectorInterface& u, VectorInterface& oldu,VectorInterface& newu,VectorInterface& z,VectorInterface& oldz,VectorInterface& f, vector<double>& DT_M, vector<double>&T)

{
PRIMALPROBLEM = false;

DTM1=0.0;
DTM2=0.0;

    DTM1=DT_M[m];
    
FIRSTDUAL = (m==_M); 
LASTDUAL  = (m==0);

GetMultiLevelSolver()->SetProblem("dp");

StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
assert(MyS);

GetMultiLevelSolver()->GetSolver()->Zero(f);

TIME=T[m]-DT_M[m]/2;
if (FIRSTDUAL || LASTDUAL)
    { if(m==_M)
    GetMultiLevelSolver()->GetSolver()->Rhs(f,-DTM1/2);
    else {
    TIME=0.0;
            DTM2=DT_M[m+1];  
    GetMultiLevelSolver()->GetSolver()->Rhs(f,-DTM2/2);
    }
}
    
else{

    DTM2=DT_M[m+1];
    GetMultiLevelSolver()->GetSolver()->Rhs(f,-(DTM1)/2);
    TIME=T[m+1]-DTM2/2;
    GetMultiLevelSolver()->GetSolver()->Rhs(f,-(DTM2)/2);   
}////////////////////////// Rechte Seite mit Gauss-Regel; TT ist mitte der Intervalle
// Std-Diskretisierung


MyS->HNDistribute(f);

// DWR - Disc
DwrFemQ1Q22d DWRFEM;
DWRFEM.BasicInit(this->_paramfile);
DiscretizationInterface* saveD = MyS->GetDiscretization();
MyS->SetDiscretization(DWRFEM,true);

TIME=T[m]-DT_M[m]/2;
if (FIRSTDUAL || LASTDUAL)
{
    if(m==_M)
    GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM1/2);
    else{
        DTM2=DT_M[m+1];
    TIME=0.0;;
GetMultiLevelSolver()->GetSolver()->Rhs(f,DTM2/2);
        }
    
}
else
    {
    DTM2=DT_M[m+1];   
    GetMultiLevelSolver()->GetSolver()->Rhs(f,(DTM1)/2);
    TIME=T[m+1]-DTM2/2;
    GetMultiLevelSolver()->GetSolver()->Rhs(f,(DTM2)/2);
    }
MyS->HNDistribute(f);
MyS->SetDiscretization(*saveD);
//Standart
MyS->GetGV(z)    =Z;
MyS->GetGV(u)    = Pu_kM[m];
MyS->GetGV(oldz) = OLDZ;

if(!FIRSTDUAL)
    {
    GetMultiLevelSolver()->GetSolver()->GetGV(newu)=Pu_kM[m+1];
    GetMultiLevelSolver()->AddNodeVector("u3",newu);
    MyS->HNAverage(newu);
    DTM2=DT_M[m+1];
    
    }
if (!LASTDUAL)
    {
    GetMultiLevelSolver()->GetSolver()->GetGV(oldu)=Pu_kM[m-1];
    GetMultiLevelSolver()->AddNodeVector("u1",oldu);
    MyS->HNAverage(oldu);
    }

MyS->HNAverage(z);
MyS->HNAverage(u);
MyS->HNAverage(oldz);
MyS->AddNodeVector("u2", u);
MyS->AddNodeVector("oldz", oldz);
MyS->Form(f,z,1.0);
MyS->DeleteNodeVector("u2");
MyS->DeleteNodeVector("oldz");
if(!FIRSTDUAL)
    {    
    GetMultiLevelSolver()->DeleteNodeVector("u3");
    }
if(!LASTDUAL)
    {
    GetMultiLevelSolver()->DeleteNodeVector("u1");
    }

MyS->HNDistribute(f);


// dwr disc


// DTM2=DTM1;



MyS->SetDiscretization(DWRFEM,true);
MyS->GetGV(z)    =Z;
MyS->GetGV(u)    = Pu_kM[m];
MyS->GetGV(oldz) = OLDZ;

if(!FIRSTDUAL)
    {
    GetMultiLevelSolver()->GetSolver()->GetGV(newu)=Pu_kM[m+1];
    GetMultiLevelSolver()->AddNodeVector("u3",newu);
    MyS->HNAverage(newu);
    
    DTM2=DT_M[m+1];
    
    
    }
if (!LASTDUAL)
    {
    GetMultiLevelSolver()->GetSolver()->GetGV(oldu)=Pu_kM[m-1];
    GetMultiLevelSolver()->AddNodeVector("u1",oldu);
    MyS->HNAverage(oldu);
    }

MyS->HNAverage(z);
MyS->HNAverage(u);
MyS->HNAverage(oldz);
MyS->AddNodeVector("u2", u);
MyS->AddNodeVector("oldz", oldz);
MyS->Form(f,z,-1.0);
MyS->DeleteNodeVector("u2");
MyS->DeleteNodeVector("oldz");
if(!FIRSTDUAL)
    {
    GetMultiLevelSolver()->DeleteNodeVector("u3");
    }
if(!LASTDUAL)
    {
    GetMultiLevelSolver()->DeleteNodeVector("u1");
    }


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
    
        DTM1=0.0;
        DTM2=0.0;
cout<<eta.sum()<<"ETA11 "<<endl;
PRIMALPROBLEM = true;
}

void Loop::EstimateNonU(DoubleVector& eta,
                vector<GlobalVector>& U, GlobalVector& Z,
            VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f,int start, int stopp)
{

StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
assert(MyS);

DWRNonLin NoLin ;

MyS->Zero(f);

double SAVETIME = TIME;
// Std Disc



int I_m=stopp-start;


//Mit Gauss
double gx = 0.5-0.5*sqrt(1.0/3.0);
for (int l=0;l<I_m; ++l)
    {
    // u im linken gauss-punkt
    MyS->GetGV(u).equ((1.0-gx), U[start+ l],gx,  U[start + l+1]);
    MyS->HNAverage(u);
    MyS->AddNodeVector("U",u); 
    MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, -DT*0.5);
    MyS->DeleteNodeVector("U");

    // u im rechten gauss-punkt
    MyS->GetGV(u).equ(gx, U[start + l],(1.0-gx),  U[start + l+1]);
    MyS->AddNodeVector("U",u);
    MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, -DT*0.5);
    MyS->DeleteNodeVector("U");

    }





// DWR Disc
DwrFemQ1Q22d  DWRFEM;
DWRFEM.BasicInit(this->_paramfile);
DiscretizationInterface* saveD = MyS->GetDiscretization();
MyS->SetDiscretization(DWRFEM,true);

// Mit Gauss


for (int l=0;l<I_m; ++l)
    {
    MyS->GetGV(u).equ((1.0-gx), U[start+ l],gx,  U[start + l+1]);
    MyS->HNAverage(u);
    MyS->AddNodeVector("U",u); 
    MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, DT*0.5);
    MyS->DeleteNodeVector("U");

    MyS->GetGV(u).equ(gx, U[start + l],(1.0-gx),  U[start + l+1]);
    MyS->AddNodeVector("U",u);
    MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, DT*0.5);
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
    cout<<eta.sum()<<"ETA5"<<endl;
}

void Loop::EstimateNonPu(DoubleVector& eta, 
            vector<GlobalVector>& U, GlobalVector& Z,
            VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f, int m,double DTM_PU)
{

StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
assert(MyS);

DWRNonLin NoLin;

MyS->Zero(f);
MyS->GetGV(u) = U[m];
MyS->GetGV(oldu) = U[m-1];


// Std Disc
MyS->HNAverage(u);

MyS->AddNodeVector("U",u); 
MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, -DTM_PU/2);
MyS->DeleteNodeVector("U");

MyS->AddNodeVector("U",oldu); 
MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, -DTM_PU/2);
MyS->DeleteNodeVector("U");

MyS->HNDistribute(f);

// DWR Disc
DwrFemQ1Q22d   DWRFEM;
DWRFEM.BasicInit(this->_paramfile);
DiscretizationInterface* saveD = MyS->GetDiscretization();
MyS->SetDiscretization(DWRFEM,true);

MyS->HNAverage(u);
MyS->AddNodeVector("U",u); 
MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, DTM_PU/2);
MyS->HNDistribute(f);
MyS->DeleteNodeVector("U");
MyS->HNAverage(u);
MyS->AddNodeVector("U",oldu); 
MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLin, DTM_PU/2);
MyS->HNDistribute(f);
MyS->DeleteNodeVector("U");
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

void Loop::EstimateRest(DoubleVector& eta, int m,
            const GlobalVector& Pu, const GlobalVector &Puold,
            const GlobalVector& Puk, const GlobalVector &Pukold,
            const GlobalVector& U, const GlobalVector &Uold, GlobalVector& Z,
            VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f)
{
cout << "Fehlerschaetzer: Restfehler"<< endl;
StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
assert(MyS);

MyS->Zero(f);
MyS->GetGV(z) = Z;
MyS->HNAverage(z);
DWRMassRhs DWRMass;

GlobalVector DUk = Puk;
DUk.add(-1.0,Pukold);
DUk.add(-1.0,U);
DUk.add( 1.0,Uold);


    
MyS->AddNodeVector("U",z);
MyS->GetDiscretization()->Rhs(MyS->GetGV(f), DWRMass, -1);
MyS->DeleteNodeVector("U");
    



MyS->HNDistribute(f);
MyS->GetDiscretization()->HNAverage(Z);
MyS->SetBoundaryVectorZero(f);


// Auswerten, kein Filtern
for(int i=0; i<Z.n(); i++)
    {
    for (int c=0; c<Z.ncomp(); c++)
    { 
    eta[i] += MyS->GetGV(f)(i,c) * DUk(i,c);
    }  
    }   

    

// DWR Disc
Q22d DWRFEM;
DWRFEM.BasicInit(this->_paramfile);
DiscretizationInterface* saveD = MyS->GetDiscretization();
MyS->SetDiscretization(DWRFEM,true);
MyS->Zero(f);

//Bleibt gleich 
GlobalVector DU = Pu;
DU.add(-1.0,Puold);
DU.add(-1.0,U);
DU.add( 1.0,Uold);

MyS->AddNodeVector("U",z);
MyS->GetDiscretization()->Rhs(MyS->GetGV(f), DWRMass, 1);
MyS->DeleteNodeVector("U");
    



MyS->HNDistribute(f);
MyS->SetBoundaryVectorZero(f);


MyS->GetDiscretization()->HNAverage(Z);

// Auswerten, kein Filtern
for(int i=0; i<Z.n(); i++)
    {
    for (int c=0; c<Z.ncomp(); c++)
    {
    eta[i] += MyS->GetGV(f)(i,c)*DU(i,c);
    }  
    }   


MyS->SetDiscretization(*saveD);
}

void Loop::EstimateNonMeanU(DoubleVector& eta, int m,
                GlobalVector& Pu, GlobalVector& Pu_k,  vector<GlobalVector>& U, vector<GlobalVector>& U_2, GlobalVector& Z,
                VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f, int start, int stopp)

{
    

    
    //-----------------------------------------------------------------

StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
assert(MyS);

DWRNonLin NoLinM;

MyS->Zero(f);


int I_m = stopp-start;


double gx = 0.5-0.5*sqrt(1.0/3.0);
// Mit gauss Q1

for (int l=0;l<I_m; ++l)
    {
    // MyS->GetGV(u) = U[(m-1)*I_m + l];
    MyS->GetGV(u).equ((1.0-gx), U[start+ l],gx,  U[start + l+1]);
    // MyS->GetGV(oldu).equ(1,Pu_k,-1*(1.0-gx), U[(m-1)*I_m + l],-gx,  U[(m-1)*I_m + l+1]);
    MyS->AddNodeVector("U",u);       
    //  MyS->AddNodeVector("W",oldu);
    //  TIME = DT * ((m-1)*I_m + l) + 0.5 * DT -  0.5*sqrt(1.0/3.0) * DT;
    MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, -DT*0.5);
    // MyS->DeleteNodeVector("Pu");

    
    MyS->DeleteNodeVector("U");
    
    MyS->GetGV(u).equ(gx, U[start + l],(1.0-gx),  U[start + l+1]);
    // MyS->GetGV(oldu).equ(1,Pu_k,-gx, U[(m-1)*I_m + l],-1*(1.0-gx),  U[(m-1)*I_m + l+1]);
    MyS->AddNodeVector("U",u);       
    //   MyS->AddNodeVector("W",oldu);
    //  TIME = DT * ((m-1)*I_m + l) + 0.5 * DT +  0.5*sqrt(1.0/3.0) * DT;
    MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, -DT*0.5);
    MyS->DeleteNodeVector("U");
    //   MyS->DeleteNodeVector("W");

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

//Gauss Q22


for (int l=0;l<I_m; ++l)
    {
    MyS->GetGV(u).equ((1.0-gx), U[start + l],gx,  U[start + l+1]);
    // MyS->GetGV(oldu).equ(1,Pu_k,-1*(1.0-gx), U[(m-1)*I_m + l],-gx,  U[(m-1)*I_m + l+1]);
    MyS->AddNodeVector("U",u);       
    //MyS->AddNodeVector("W",oldu);
    
    // TIME = DT * ((m-1)*I_m + l) + 0.5 * DT -  0.5*sqrt(1.0/3.0) * DT;
    MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, DT*0.5);
    MyS->DeleteNodeVector("U");
    //MyS->DeleteNodeVector("W");

    MyS->GetGV(u).equ(gx, U[start + l],(1.0-gx),  U[start + l+1]);
    //MyS->GetGV(oldu).equ(1,Pu_k,-gx, U[(m-1)*I_m + l],-1*(1.0-gx),  U[(m-1)*I_m + l+1]);
    MyS->AddNodeVector("U",u);       
    // MyS->AddNodeVector("W",oldu);
    //  TIME = DT * ((m-1)*I_m + l) + 0.5 * DT +  0.5*sqrt(1.0/3.0) * DT;
    MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, DT*0.5);
    MyS->DeleteNodeVector("U");
    //  MyS->DeleteNodeVector("W");

    }

MyS->HNDistribute(f);

MyS->SetBoundaryVectorZero(f);

const GlobalVector& F1 = MyS->GetGV(f);

MyS->GetDiscretization()->HNAverage(Z);
MyS->GetDiscretization()->HNAverage(Pu);


// Auswerten, kein Filtern
for(int i=0; i<Z.n(); i++)
    {
    for (int c=0; c<Z.ncomp(); c++)
    {
    eta[i] += F1(i,c)*Z(i,c);
    }  
    }
MyS->SetDiscretization(*saveD);

}


void Loop::EstimateNonMeanPu(DoubleVector& eta, int m,
                GlobalVector& Pu,  vector<GlobalVector>& Pu_k,vector<GlobalVector>& U, GlobalVector& Z,
                VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f, double DTM_U)
{

StdSolver* MyS = dynamic_cast<StdSolver*> (GetMultiLevelSolver()->GetSolver());
assert(MyS);



// DWRNonLinMatrix NoLinM;

DWRNonLin NoLinM;



MyS->Zero(f);
MyS->GetGV(u) = Pu_k[m];
MyS->GetGV(oldu) = Pu_k[m-1];


MyS->HNAverage(u);

MyS->AddNodeVector("U",u); 
MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, -DTM_U/2);
MyS->DeleteNodeVector("U");

MyS->AddNodeVector("U",oldu); 
MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, -DTM_U/2);
MyS->DeleteNodeVector("U");



MyS->HNDistribute(f);

MyS->SetBoundaryVectorZero(f);

MyS->GetDiscretization()->HNAverage(Z);


// Auswerten, kein Filtern
for(int i=0; i<Z.n(); i++)
    {
    for (int c=0; c<Z.ncomp(); c++)
    {
    eta[i] += MyS->GetGV(f)(i,c)*Z(i,c);
    }  
    }

// DWR Disc
Q22d DWRFEM;
MyS->Zero(f);  
MyS->GetGV(u) = Pu;

DWRFEM.BasicInit(this->_paramfile);
DiscretizationInterface* saveD = MyS->GetDiscretization();
MyS->SetDiscretization(DWRFEM,true);

MyS->HNAverage(u);
MyS->AddNodeVector("U",u); 
MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, DTM_U/2);
MyS->HNDistribute(f);
MyS->DeleteNodeVector("U");

MyS->HNAverage(u);
MyS->AddNodeVector("U",oldu); 
MyS->GetDiscretization()->Rhs(MyS->GetGV(f), NoLinM, DTM_U/2);
MyS->HNDistribute(f);
MyS->DeleteNodeVector("U");



MyS->SetBoundaryVectorZero(f);
MyS->GetDiscretization()->HNAverage(Z);
// Auswerten, kein Filtern
for(int i=0; i<Z.n(); i++)
    {
    for (int c=0; c<Z.ncomp(); c++)
    {
    eta[i] +=MyS->GetGV(f)(i,c)*Z(i,c);
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
                vector<GlobalVector>& U_2,
                vector<GlobalVector>& Ztotal,
                vector<GlobalVector>& Pu_kM,
                vector<GlobalVector>& Pu_M,
                VectorInterface& u,
                VectorInterface& oldu,
                VectorInterface& newu,
                VectorInterface& z,
                VectorInterface& oldz,
                VectorInterface& f,
                vector<double>& DT_M,
                vector<double>& T,
        vector<double>& eta_time
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
// EstimateInitial(eta0,Utotal[0], Ztotal[0], u,oldu,z,f);

for (int m=1;m<=_M;++m)
    {
    cout << "m="<< m << "\\Fehlerschaetzer DWR" << endl;
    
    

    eta.zero();
    eta1.zero();
    eta11.zero(),eta2.zero(),eta22.zero(),eta23.zero(),eta3.zero(),eta4.zero(),eta5.zero();

    ///////////// Teil 1 - DWR-Anteil     
    EstimateDWRprim(eta1, m, Pu_kM[m],Utotal,Ztotal[m], u,oldu,z,f,T);
    cout<<eta1.sum()<<"ETA1 "<<endl;

    /// Teil 1.1 duales residuum

    
    EstimateDWRdual(eta11,m,Pu_kM, Pu_kM[m], Ztotal[m+1], Ztotal[m],u, oldu,newu,z, oldz,f,DT_M,T);
    
    //Zeitterme

    int start=T[m-1]/DT+1.e-10;
    int stoppi=T[m]/DT+1.e-10;	
    
        EstimateAvg(eta2, Pu_M[m], Pu_M[m-1],  Utotal[stoppi], Utotal[start], Ztotal[m], u,oldu,z,f);
    cout<<eta2.sum()<<"ETA2 "<<endl;

    
        EstimateRest(eta3, m,Pu_M[m], Pu_M[m-1], Pu_kM[m], Pu_kM[m-1], Utotal[stoppi], Utotal[start], Ztotal[m], u, oldu, z,f);
    cout<<eta3.sum()<<"ETA3 "<<endl;	

    // Nichtlinearitaet

    EstimateNonU(eta22,Utotal, Ztotal[m], u,oldu,z,f,T[m-1]/DT+1.e-10,T[m]/DT+1.e-10);
    cout<<eta22.sum()<<"ETA22 "<<endl;

    EstimateNonPu(eta23, Pu_kM, Ztotal[m], u,oldu,z,f,m,DT_M[m]);
    cout<<eta23.sum()<<"ETA23 "<<endl;

    EstimateNonMeanU(eta5, m, Pu_M[m],Pu_kM[m],Utotal,U_2, Ztotal[m], u,oldu,z,f,T[m-1]/DT+1.e-10,T[m]/DT+1.e-10);
    cout<<eta5.sum()<<"ETA5"<<endl;

    //Teil 4 Fehler vom geittelten Problem zu Pu
        EstimateNonMeanPu(eta4, m, Pu_M[m],Pu_kM,Utotal, Ztotal[m], u,oldu,z,f,DT_M[m]);
    cout<<eta4.sum()<<"ETA4"<<endl;

    eta_time[m-1]=0.5*eta1.sum()+0.5*eta11.sum()+eta2.sum()+eta22.sum()-eta23.sum()+eta3.sum()+eta4.sum()+eta5.sum();
    

    
    }

// eta.add(1.0,eta0);
eta.add(0.5,eta11);

eta.add(0.5,eta1);
eta.add(1.0,eta2);
eta.add(1.0,eta22);
eta.add(-1.0,eta23);
eta22.add(-1.0,eta23);
cout<<eta22.sum()<<"combi"<<endl;
eta.add(1.0,eta3);
eta.add(1.0,eta4);
eta.add(-1.0,eta5);



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
double DTM = endtime / _M;
cout << "N="<< _niter << "\t M=" << _M << "\t dt=" << DT << "\t dtm=" << DTM << endl;

assert( fabs(endtime - _M*DTM)<1.e-8);
assert( fabs(endtime - _niter * DT)<1.e-8);

Extrapolator Extra;

for (int ADAITER=0;ADAITER<8;++ADAITER)
    {
    TIME=0.0;
    DTM1=0.0;
    DTM2=0.0;

    // vectors for solution and right hand side
    VectorInterface u("u"), 
    newu("newu"), f("f"), 
    oldu("oldu"),z("z"),oldz("oldz"),h("h"),oldh("oldh"),w("w"), oldw("oldw");

    
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
    

    GetMultiLevelSolver()->SetProblem("LaplaceT");


    InitSolution(u);
    
    GetMultiLevelSolver()->Equ(oldu,1.0,u);
    Utotal.push_back(MyS->GetGV(u));
    GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,0);
    
    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    
    
    // primale Probleme loesen
    SolvePrimalProblem(Utotal,u,oldu,f,Htotal, h, oldh, ADAITER);
    
    assert(Utotal.size() == _niter+1);



    double Ji=CompFunctional(Utotal,u,f);
    nvector<double> J(1,Ji);
    Extra.NewValues(J);
    cout << "Integral ueber J = "<< Ji << endl;
    Extra.Print();
    
    

    // Mitteln
    
    
    assert(_niter%2==0);
    int K=_niter/2;
    

    vector<GlobalVector> U_2(K+1,GlobalVector(MyS->GetGV(u).ncomp(), MyS->GetGV(u).n()));
    vector<GlobalVector> Pu_kM(_M+1,GlobalVector(MyS->GetGV(u).ncomp(), MyS->GetGV(u).n()));
    vector<GlobalVector> Pu_M (_M+1,GlobalVector(MyS->GetGV(u).ncomp(), MyS->GetGV(u).n()));
    Pu_M[0]=Utotal[0];
    Pu_kM[0]=Utotal[0];
    U_2[0]=Utotal[0];
    
    
        vector<double> DT_M(_M+1);
        vector<double> T(_M+1);
        
    T[0]=0;
    T[1]=2.5;
    T[2]=3.0;
    T[3]=3.5;
    T[4]=4.0;
    
    
    double checki=0.0;
    for( int k=1;k<=_M;k++)
    {
        DT_M[k]=T[k]-T[k-1];
    checki+=DT_M[k];
    }
        assert( checki== endtime);
    
    

    
    // Intregartion und Mittelung
    for (int m=1;m<=_M;++m)
    {  
    //void TrapezInt(GlobalVector& avg, const vector<GlobalVector>& U, int start, int stopp);
    MittelInt(Pu_kM[m-1],Pu_kM[m], Utotal, T[m-1]/DT+1.e-10, T[m]/DT+1.e-10,DT_M[m]);
    MittelInt(Pu_M[m-1],Pu_M[m], Utotal, T[m-1]/DT+1.e-10, T[m]/DT+1.e-10,DT_M[m]);

    GetMultiLevelSolver()->GetSolver()->GetGV(oldu)=Pu_M[m];
    nvector<double> functionals_P = Functionals(oldu,f);
    GetMultiLevelSolver()->GetSolver()->Visu("Results/PU_kM",oldu,m);
    }
    
    
    // Duales Problem loesen. Alle dualen Schritte speichern. 
    vector<GlobalVector> Ztotal(_M+2, GlobalVector(MyS->GetGV(u).ncomp(), MyS->GetGV(u).n()));
    Ztotal[_M+1].zero();
    
    cout<<"dualsingel"<<endl;
    vector<GlobalVector> Wtotal(_M+2, GlobalVector(MyS->GetGV(h).ncomp(), MyS->GetGV(h).n()));
        cout<<"dualsingel"<<endl;
    Wtotal[_M+1].zero();
    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    cout<<"dualsingel"<<endl;
    SolveDualProblem(Ztotal,Wtotal,Htotal,f,u, oldu,newu,z, oldz,w,oldw,h,oldh,Pu_kM,ADAITER,DT_M,T);

    abort();
    // Fehlerschatzer
    int nnodes = MyS->GetMesh()->nnodes();
    DoubleVector eta(nnodes,0.0), eta1(nnodes,0.0),eta11(nnodes,0.0),eta2(nnodes,0.0),eta22(nnodes,0.0),eta23(nnodes,0.0),eta3(nnodes,0.0),eta4(nnodes,0.0),eta5(nnodes,0.0),eta0(nnodes,0.0),zeta(nnodes,0.0);

    vector<double> eta_time(_M,0.0);

    EstimateDualError(eta,eta0,eta1,eta11,eta2,eta22,eta23,eta3,eta4,eta5, Utotal,U_2, Ztotal, Pu_kM, Pu_M,u,oldu,newu,z,oldz,f,DT_M,T,eta_time);
    

    this->EtaVisu("Results/eta",ADAITER,eta);

    stringstream str;
    str << "eta_3_3.5.txt";
    ofstream OUTF(str.str().c_str(),ios::app);
    OUTF.precision(10);

    // OUTF << GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes() << " " << DTM  << " " << J <<  " " << eta.sum()  <<" " << eta0.sum()<<" " << eta1.sum()<<" " << eta11.sum()<<" " << eta2.sum() <<" " << eta22.sum()<<" "<< eta3.sum()<< " "<< eta4.sum()<< " "<< eta5.sum()<<endl;
    OUTF<< GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes() << " " << eta_time[0]<<" "<<eta_time[1]<<" "<< eta_time[2]<<" "<< eta_time[3]<<endl;
    
    GetMultiLevelSolver()->GetSolver()->Visu("Results/neuu",u,ADAITER);
    // Gitter verfeinern

    IntVector refcells, coarsecells;
    assert(eta.size()==GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes());
    double mean = eta.norm_l1()/(1.0*eta.size());
    for (int i=0;i<eta.size();++i)
    //	if (fabs(eta[i])>1.5*mean)
    refcells.push_back(i);
    GetMeshAgent()->refine_nodes(refcells);


    PrintMeshInformation();
    //  GetMultiLevelSolver()->GetSolver()->Visu("Results/neuu",u,ADAITER);
    //GetMultiLevelSolver()->GetSolver()->Visu("Results/PU",oldu,ADAITER);
    } 
    


}






