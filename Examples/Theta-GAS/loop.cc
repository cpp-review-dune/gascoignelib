#include "loop.h"
#include "stdtimesolver.h"
#include  "filescanner.h"
#include  "solvers.h"


#include <fstream>
#include <string>
#include <iostream>
using namespace std;
using namespace Gascoigne;

double SHIFT = 0.0;


double __GLOBAL_TIME,__GLOBAL_DT;



//double exact_time = 0.4493560614; // 1.5
//double exact_time = 0.4505536071;  // 2.0
//double exact_time =  0.4471761764; //0.1435160481; // 4.0


/////////// so klappts, beta wechselt einmal
//double exact_time = 0.4931128801; // extra p = 1.5
double exact_time;

//0.9477662297; // extra p = 1.5



///// Knotenbasispolynome vom Grad 2 in Stuetzstellen x0, x1, x2
double PQ0(double x,double x0,double x1,double x2)  { return (x-x1)*(x-x2)/(x0-x1)/(x0-x2); }
double PQ1(double x,double x0,double x1,double x2)  { return (x-x0)*(x-x2)/(x1-x0)/(x1-x2); }
double PQ2(double x,double x0,double x1,double x2)  { return (x-x0)*(x-x1)/(x2-x0)/(x2-x1); }

double Q(double x,double x0,double x1,double x2,double w0,double w1,double w2)
{
  return PQ0(x,x0,x1,x2) * w0 +PQ1(x,x0,x1,x2) * w1 +PQ2(x,x0,x1,x2) * w2;
}

double L(double x,double x0,double x1, double w0,double w1)
{
  return (x-x1)/(x0-x1) * w0 + (x-x0)/(x1-x0) * w1;
}

// {
//   double g0 = t1 + (t2-t1)*(0.5+sqrt(1.0/12.0));
//   double g1 = t2 - (t2-t1)*(0.5+sqrt(1.0/12.0));
  
//   double Q0 = Q(g0, t0,t1,t2, z0,z2,z4);
//   double Q1 = Q(g1, t0,t1,t2, z0,z2,z4);

//   double L0 = L(g0, t1,t2, z2,z3);
//   double L1 = L(g1, t1,t2, z2,z3);
    
//   double i = sqrt(0.5 * (t2-t1) * (pow(Q0-L0,2.0) + pow(Q1-L1,2.0)));
// }




///// ableitungen davon.
double DPQ0(double x,double x0,double x1,double x2) { return ((x-x2)+(x-x1))/(x0-x1)/(x0-x2); }
double DPQ1(double x,double x0,double x1,double x2) { return ((x-x2)+(x-x0))/(x1-x0)/(x1-x2); }
double DPQ2(double x,double x0,double x1,double x2) { return ((x-x0)+(x-x1))/(x2-x0)/(x2-x1); }

///// Knotenbasispolynome vom Grad 2 in Stuetzstellen x0, x1, x2
double P3Q0(double x,double x0,double x1,double x2, double x3)  { return (x-x1)*(x-x2)*(x-x3)/(x0-x1)/(x0-x2)/(x0-x3); }
double P3Q1(double x,double x0,double x1,double x2, double x3)  { return (x-x0)*(x-x2)*(x-x3)/(x1-x0)/(x1-x2)/(x1-x3); }
double P3Q2(double x,double x0,double x1,double x2, double x3)  { return (x-x0)*(x-x1)*(x-x3)/(x2-x0)/(x2-x1)/(x2-x3); }
double P3Q3(double x,double x0,double x1,double x2, double x3)  { return (x-x0)*(x-x1)*(x-x2)/(x3-x0)/(x3-x1)/(x3-x2); }
///// ableitungen davon.
double DP3Q0(double x,double x0,double x1,double x2, double x3)  { return ((x-x2)*(x-x3)+(x-x1)*(x-x3)+(x-x1)*(x-x2))/(x0-x1)/(x0-x2)/(x0-x3); }
double DP3Q1(double x,double x0,double x1,double x2, double x3)  { return ((x-x2)*(x-x3)+(x-x0)*(x-x3)+(x-x0)*(x-x2))/(x1-x0)/(x1-x2)/(x1-x3); }
double DP3Q2(double x,double x0,double x1,double x2, double x3)  { return ((x-x1)*(x-x3)+(x-x0)*(x-x3)+(x-x0)*(x-x1))/(x2-x0)/(x2-x1)/(x2-x3); }
double DP3Q3(double x,double x0,double x1,double x2, double x3)  { return ((x-x1)*(x-x2)+(x-x0)*(x-x2)+(x-x0)*(x-x1))/(x3-x0)/(x3-x1)/(x3-x2); }

void InterpolateLinear(GlobalVector& U, const GlobalVector& U0, const GlobalVector& U1, double x, double x0, double x1) 
{
  U.ncomp()=U0.ncomp(); U.resize(U0.n());
  U.zero(); 
  U.add( (x-x1)/(x0-x1), U0);  
  U.add( (x-x0)/(x1-x0), U1);
}


void InterpolateQuadratic(GlobalVector& U, const GlobalVector& U0, const GlobalVector& U1, const GlobalVector& U2, double x, double x0, double x1, double x2) 
{
  U.ncomp()=U0.ncomp(); U.resize(U0.n());
  U.zero(); 
  U.add(PQ0(x,x0,x1,x2), U0);  
  U.add(PQ1(x,x0,x1,x2), U1);  
  U.add(PQ2(x,x0,x1,x2), U2); 
}
void InterpolateQuadraticDerivative(GlobalVector& U, const GlobalVector& U0, const GlobalVector& U1, const GlobalVector& U2, double x, double x0, double x1, double x2) 
{
  U.ncomp()=U0.ncomp(); U.resize(U0.n());
  U.zero(); 
  U.add(DPQ0(x,x0,x1,x2), U0);  
  U.add(DPQ1(x,x0,x1,x2), U1);  
  U.add(DPQ2(x,x0,x1,x2), U2); 
}
void InterpolateCubic(GlobalVector& U, const GlobalVector& U0, const GlobalVector& U1, const GlobalVector& U2, const GlobalVector& U3, double x, double x0, double x1, double x2, double x3) 
{
  U.ncomp()=U0.ncomp(); U.resize(U0.n());
  U.zero(); 
  U.add(P3Q0(x,x0,x1,x2,x3), U0);  
  U.add(P3Q1(x,x0,x1,x2,x3), U1);  
  U.add(P3Q2(x,x0,x1,x2,x3), U2); 
  U.add(P3Q3(x,x0,x1,x2,x3), U3); 
}
void InterpolateCubicDerivative(GlobalVector& U, const GlobalVector& U0, const GlobalVector& U1, const GlobalVector& U2, const GlobalVector& U3, double x, double x0, double x1, double x2, double x3) 
{
  U.ncomp()=U0.ncomp(); U.resize(U0.n());
  U.zero(); 
  U.add(DP3Q0(x,x0,x1,x2,x3), U0);  
  U.add(DP3Q1(x,x0,x1,x2,x3), U1);  
  U.add(DP3Q2(x,x0,x1,x2,x3), U2); 
  U.add(DP3Q3(x,x0,x1,x2,x3), U3); 
}

double ThetaLoop::GetEndFunctional(VectorInterface& f,VectorInterface& u)
{
  SetTimeData(__Tstop,0,0);
  assert(__primal_solution.size()==_niter+1);
  GetMultiLevelSolver()->GetSolver()->GetGV(u) = __primal_solution[_niter];
  return GetMultiLevelSolver()->ComputeFunctional(f,u,"end");
}

double ThetaLoop::GetTimeFunctional(VectorInterface& f,VectorInterface& u)
{
  assert(__TIME.size()==_niter+1);

  double GX[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};
  double res = 0.0;
  for (int i=0;i<_niter;++i)
    {
      double gx[2] = { __TIME[i] + __DT[i] * GX[0], __TIME[i] + __DT[i] * GX[1] };
      SetTimeData(gx[0],__DT[i], __THETA[i]);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(u) = __primal_solution[i];
      GetMultiLevelSolver()->GetSolver()->GetGV(u).sadd(GX[1], GX[0], __primal_solution[i+1]);
      double j1 = 0.5 * __DT[i] * GetMultiLevelSolver()->ComputeFunctional(f,u,"time");
      res += j1;
      SetTimeData(gx[1],__DT[i], __THETA[i]);
      
      GetMultiLevelSolver()->GetSolver()->GetGV(u) = __primal_solution[i];
      GetMultiLevelSolver()->GetSolver()->GetGV(u).sadd(GX[0], GX[1], __primal_solution[i+1]);
      double j2 = 0.5 * __DT[i] * GetMultiLevelSolver()->ComputeFunctional(f,u,"time");
      res += j2;
    }
  
  
  return res;
}


void ThetaLoop::BasicInit(const ParamFile* paramfile, const ProblemContainer* PC, const FunctionalContainer* FC)
{  
  GetMultiLevelSolverPointer() = new ThetaMultiLevelSolver();
  
  StdLoop::BasicInit(paramfile, PC, FC);
}


void ThetaLoop::SetTimeData(double time, double dt, double theta)
{
  __GLOBAL_TIME = time;
  __GLOBAL_DT   = dt;
  for(int l=0;l<GetMultiLevelSolver()->nlevels();l++)
    {
      StdTimeSolver* TS = dynamic_cast<StdTimeSolver*>(GetMultiLevelSolver()->GetSolver(l));
      assert(TS);
      TS->SetTimeData(dt,theta,time, 0,0);
    }
}

void ThetaLoop::PrimalLoop(VectorInterface& u, VectorInterface& f)
{
  SetTimeData(__TIME[0], __DT[0], __THETA[0]);
  InitSolution(u);
 
   ThetaSolver* TS = dynamic_cast<ThetaSolver*>(GetMultiLevelSolver()->GetSolver());

   ofstream PRIMAL("primal.txt");
   PRIMAL.precision(10);
  

  __primal_solution.clear();
  __primal_solution.push_back(TS->GetGV(u));

   //char s[20];
   //sprintf(s,"primal_%i.txt", _niter);
   //ofstream OUT(s);
   //OUT << __TIME[0] << "\t" << TS->GetGV(u)(4,0) << endl;
  stringstream filename;
  filename << "Results/u.reflev"<<reflevel;
    TS->Visu(filename.str(),u,0); 
    //PRIMAL << __TIME[0] << "\t" << TS->GetGV(u)(0,0) << endl;
  for (_iter=1; _iter<=_niter; _iter++)
    {
      assert(_iter<__TIME.size());
      assert(_iter-1<__DT.size());
      assert(_iter-1<__THETA.size());
      double dt    = __DT[_iter-1];
      double time  = __TIME[_iter-1];
      double theta = __THETA[_iter-1];
      cout<<"--------------------------------------------------------"<<endl;

      cout << "step " << _iter << "\t time:" << time << " -> "
	   << __TIME[_iter] << "\t dt = " << dt << "\t theta = "  << theta << endl;

      TS->Zero(f);
      SetTimeData(time,dt,theta);
      TS->RhsOld(f,u, time,dt,theta);

      SetTimeData(time+dt,dt,theta);
      TS->RhsNew(f, time+dt,dt,theta);

      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);

      
      string status = GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());
            TS->Visu(filename.str(),u,_iter);
      //      OUT << __TIME[_iter] << "\t" << TS->GetGV(u)(4,0) << endl;
      //      PRIMAL << __TIME[_iter] << "\t" << TS->GetGV(u)(0,0) << endl;
      
      __primal_solution.push_back(TS->GetGV(u));
    }
   // PRIMAL.close();
   // OUT.close();
  
}


// --------------------------------------------------

/**
 *  Das duale Problem wird auf der Basis des kontinuierlichen Theta-Verfahrens berechnet,
 *  nicht als duales Theta-Verfahren.
 *
 *  Denn, wir brauchen zur Approximation sowieso (dwr:2) Galerkin-Orthogonalitaet.
 *  Die haetten wir dann im dualen und die Approximation (dwr:3.5) muesste dann nur
 *  fuer das primale durchgefuert werden.
 *
 *  Ausserdem benoetigen wir fuer den Interpolations-fehler z-i_k z eine gute approximation
 *  an z und das ist durch die richtige Galerkin-formulierung definiert.
 *
 ***
 **  Doch anders. duales theta-schema.
 *
 **/
void ThetaLoop::AdjointLoop(VectorInterface& z, VectorInterface& f, VectorInterface& u)
{
  assert(__primal_solution.size()==_niter+1);
  __dual_solution.clear();
  __dual_solution.resize(_niter+1);

  ThetaSolver* TS = dynamic_cast<ThetaSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TS);
  ThetaMultiLevelSolver* TMLS = dynamic_cast<ThetaMultiLevelSolver*>(GetMultiLevelSolver());


  // char s[20];
  // sprintf(s,"dual_%i.txt", _niter);
  // ofstream DUAL(s);
  // DUAL.precision(10);
  
  double GX[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};
  
  stringstream filename;
  filename << "Results/z.reflev"<<reflevel;
    
  assert(__TIME.size()==_niter+1);
  for (_iter=_niter; _iter>=0; --_iter)
    {
      // Intervall     __TIME[_iter] -> __TIME[_iter-1]
      cout<<"--------------------------------------------------------"<<endl;
      cout << "step " << _iter << "\t solve adjoint at time:" << __TIME[_iter] << endl;

      // RHS
      TS->Zero(f);

      // Schritt i, testfunktion phi_i, linearisiert um u[i]
      GetMultiLevelSolver()->GetSolver()->GetGV(u) = __primal_solution[_iter];
      SetTimeData(__TIME[_niter],0,0); // ???????
      // Adjoint Solution in the two GP of this Interval
      vector<GlobalVector> GU(2);


      //////////////////////////////////////////////// FUNKTIONAL t=T
      // if (_iter==_niter) // 1-ter Schritt
      // 	TS->RhsEndFunctional(f,__primal_solution[_niter]);
      
      if (_iter<_niter) // zeitschritt  z[i+1] - k (1-theta) a' 
	{
	  GetMultiLevelSolver()->GetSolver()->GetGV(z) = __dual_solution[_iter+1];
	  TS->AdjointRhsOld(f,z,u, __TIME[_iter], __DT[_iter], __THETA[_iter]);

	  double gx[2] = { __TIME[_iter] + __DT[_iter] * GX[0], __TIME[_iter] + __DT[_iter] * GX[1] };
	  
	  GU[0] = __primal_solution[_iter]; GU[0].sadd(GX[1], GX[0], __primal_solution[_iter+1]);
	  GU[1] = __primal_solution[_iter]; GU[1].sadd(GX[0], GX[1], __primal_solution[_iter+1]);

	  SetTimeData(gx[0],0,0);
	  TS->AdjointRhsFunctional(f, GU[0], 0.5 * __DT[_iter] * GX[1]);
	  SetTimeData(gx[1],0,0);
	  TS->AdjointRhsFunctional(f, GU[1], 0.5 * __DT[_iter] * GX[0]);
	}
      if (_iter>0)
	{
	  double gx[2] = { __TIME[_iter-1] + __DT[_iter-1] * GX[0], __TIME[_iter-1] + __DT[_iter-1] * GX[1] };

	  GU[0] = __primal_solution[_iter-1]; GU[0].sadd(GX[1], GX[0], __primal_solution[_iter]);
	  GU[1] = __primal_solution[_iter-1]; GU[1].sadd(GX[0], GX[1], __primal_solution[_iter]);

	  SetTimeData(gx[0],0,0);
	  TS->AdjointRhsFunctional(f, GU[0], 0.5 * __DT[_iter-1] * GX[0]);
	  SetTimeData(gx[1],0,0);
	  TS->AdjointRhsFunctional(f, GU[1], 0.5 * __DT[_iter-1] * GX[1]);
 	}
      TS->SetBoundaryVectorZero(f);
      
      
      //// Aufbau Matrix & Loesen. Der letzte schritt ist anders, 
      //// nicht nur die Matrix sondern auch das residuum!!!
      GetMultiLevelSolver()->GetSolver()->GetGV(u) = __primal_solution[_iter];
      if (_iter>0)
	{
	  SetTimeData(__TIME[_iter], __DT[_iter-1], __THETA[_iter-1]);
      	  TS->AssembleAdjointMatrix(u, __DT[_iter-1], __THETA[_iter-1]);
	  GetMultiLevelSolver()->AddNodeVector("u",u);
	  TMLS->SolveAdjoint(z,f,GetSolverInfos()->GetLInfo());
	  GetMultiLevelSolver()->DeleteNodeVector("u");
	}      
      else
	{
	  SetTimeData(__TIME[0], __DT[0], __THETA[0]);
	  TS->AssembleAdjointMatrixFirst(u);
	  GetMultiLevelSolver()->AddNodeVector("u",u);
	  TMLS->SolveAdjointFirst(z,f,GetSolverInfos()->GetLInfo());
	  GetMultiLevelSolver()->DeleteNodeVector("u");
	}
      
      //      DUAL << __TIME[_iter] << "\t" << TS->GetGV(z)(4,0) << endl;
      

      //      TS->Visu("Results/z",z,_iter);
      __dual_solution[_iter] = TS->GetGV(z);
          TS->Visu(filename.str(),z,_iter);

    }
  //  DUAL.close();
  
}

double ThetaLoop::EstimateInt(DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f)
{   
  assert(eta.size()==_niter);
  assert(_niter+1==__primal_solution.size());
  assert(_niter+1==__dual_solution.size());

  ThetaSolver* TS = dynamic_cast<ThetaSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TS);
  
  // integration of Galerkin residual with 2-point Gauss rule
  double GX[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};
  
  // solution and adjoint solution in Gauss-Points
  vector<GlobalVector> GU; GU.resize(2);
  vector<GlobalVector> GZ; GZ.resize(2);

  double res = 0.0;
  
  for (int i=0;i<_niter;++i)
    {
      double time  = __TIME[i];
      double dt    = __DT[i];
      double theta = __THETA[i];

      // transformed GP
      double gx[2] = { time + dt * GX[0], time + dt * GX[1] };
      
      // solution  in GP 
      //////// was passiert hier bei zeitabhaengigen Randwerten???
      //////// vielleicht muss da ein settime(time + dt gx[i]) und setboundaryvalue hin
      GU[0] = __primal_solution[i]; GU[0].sadd(GX[1], GX[0], __primal_solution[i+1]);
      GU[1] = __primal_solution[i]; GU[1].sadd(GX[0], GX[1], __primal_solution[i+1]);
      
      // weights for adjoint solution in GP
      double gzw[2] = {1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * GX[0] - 1.0) ,
		       1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * GX[1] - 1.0) };

      TS->Zero(f);
      
      // (dt u, z)
      TS->GetGV(u) = __primal_solution[i+1];  
      TS->GetMassMatrix()->vmult_time(TS->GetGV(f), TS->GetGV(u), TS->GetTimePattern(),  1.0);
      TS->GetGV(u) = __primal_solution[i  ];  
      TS->GetMassMatrix()->vmult_time(TS->GetGV(f), TS->GetGV(u), TS->GetTimePattern(), -1.0);
      
      // a(u)(z)
      TS->GetGV(u) = GU[0];
      SetTimeData(gx[0], dt, theta);
      TS->Rhs(f,-0.5 * dt);
      TS->FormOnly(f,u, 0.5 * dt * gzw[0]);

      TS->GetGV(u) = GU[1];
      SetTimeData(gx[1], dt, theta);
      TS->Rhs(f,-0.5 * dt);
      TS->FormOnly(f,u, 0.5 * dt * gzw[1]);
      
      TS->GetGV(z) = __dual_solution[i+1];
      TS->SetBoundaryVectorZero(f);
      double r = TS->ScalarProduct(f,z);

      res += r;
      eta[i] += r;
    }  
  return res;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    Estimator for the FS-theta-scheme

double ThetaLoop::EstimatePrimalTheta(DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f)
{
  assert(eta.size()==_niter);
  assert(_niter+1==__primal_solution.size());
  assert(_niter+1==__dual_solution.size());
  
  ThetaSolver* TS = dynamic_cast<ThetaSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TS);
  
  // integration of Galerkin residual with 2-point Gauss rule
  double GX[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};

  // solution and adjoint solution in Gauss-Points
  vector<GlobalVector> GU; GU.resize(2);
  vector<GlobalVector> GZ; GZ.resize(2);

  double res = 0.0;
  for (int i=0;i<_niter;++i)
    {
      double time  = __TIME[i];
      double dt    = __DT[i];
      double theta = __THETA[i];
      // transformed GP
      double gx[2] = { time + dt * GX[0], time + dt * GX[1] };

      double r = 0.0;
      
      // uk in Gauss-Points
      GU[0] = __primal_solution[i]; GU[0].sadd(GX[1], GX[0], __primal_solution[i+1]);
      GU[1] = __primal_solution[i]; GU[1].sadd(GX[0], GX[1], __primal_solution[i+1]);

      
      // weights for adjoint solution in GP
      double gzw[2] = {1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * GX[0] - 1.0) ,
		       1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * GX[1] - 1.0) };
      

      //////////// approximation using big patches.
      
      int i0 = i-i%2;
      // if ((i0>0)&&(i0+1<__DT.size()))
      // 	  if (__DT[i0]!=__DT[i0+1]) --i0;
      assert(i0+1 < __dual_solution.size());
      assert(i0+1 < __DT.size());
      assert(__DT[i0]==__DT[i0+1]);
      assert(__THETA[i0]==__THETA[i0+1]);
      assert(i0+2 < __dual_solution.size());
      
      

      // InterpolateLinear(GZ[0], __dual_solution[i0], __dual_solution[i0+3], gx[0], __TIME[i0]-0.5*__DT[i0], __TIME[i0+3]-0.5*__DT[i0+3]);
      // InterpolateLinear(GZ[1], __dual_solution[i0], __dual_solution[i0+3], gx[1], __TIME[i0]-0.5*__DT[i0], __TIME[i0+3]-0.5*__DT[i0+3]);
      InterpolateLinear(GZ[0], __dual_solution[i0+1], __dual_solution[i0+2], gx[0], __TIME[i0]+1.0*__DT[i0], __TIME[i0+1]+1.0*__DT[i0+1]);
      InterpolateLinear(GZ[1], __dual_solution[i0+1], __dual_solution[i0+2], gx[1], __TIME[i0]+1.0*__DT[i0], __TIME[i0+1]+1.0*__DT[i0+1]);
      GZ[0].add(-gzw[0], __dual_solution[i+1]);
      GZ[1].add(-gzw[1], __dual_solution[i+1]);

      
      ///////// (dt u, z)
      // GP1
      TS->Zero(f);
      SetTimeData(gx[0], dt, theta);
      TS->GetGV(u) = __primal_solution[i+1]; TS->GetGV(u).add(-1.0,__primal_solution[i]);
      TS->GetMassMatrix()->vmult_time(TS->GetGV(f), TS->GetGV(u), TS->GetTimePattern(),  0.5);

      TS->Rhs(f,-0.5 * dt);
      TS->GetGV(u) = GU[0]; 
      TS->FormOnly(f,u, 0.5 * dt);
      TS->GetGV(z) = GZ[0];
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,z);

      // GP1
      TS->Zero(f);
      SetTimeData(gx[1], dt, theta);
      TS->GetGV(u) = __primal_solution[i+1]; TS->GetGV(u).add(-1.0,__primal_solution[i]);
      TS->GetMassMatrix()->vmult_time(TS->GetGV(f), TS->GetGV(u), TS->GetTimePattern(),  0.5);

      TS->Rhs(f,-0.5 * dt);
      TS->GetGV(u) = GU[1]; 
      TS->FormOnly(f,u, 0.5 * dt);
      TS->GetGV(z) = GZ[1];
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,z);

      res += r;
      eta[i] += r;
    }

  
  return res;
}

double ThetaLoop::EstimateDualTheta(DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f)
{
  //  assert(_niter%2==0);
  
  assert(eta.size()==_niter);
  assert(_niter+1==__primal_solution.size());
  assert(_niter+1==__dual_solution.size());
  
  ThetaSolver* TS = dynamic_cast<ThetaSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TS);
  
  // integration of Galerkin residual with 2-point Gauss rule
  double GX[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};

  // solution and adjoint solution in Gauss-Points
  vector<GlobalVector> GU;   GU.resize(2, __primal_solution[0]);
  vector<GlobalVector> GUH;  GUH.resize(2, __primal_solution[0]);
  vector<GlobalVector> GDUH; GDUH.resize(2, __primal_solution[0]);
  vector<GlobalVector> GZ;   GZ.resize(2, __dual_solution[0]);
  double res = 0.0;
  for (int i=0;i<_niter;++i)
    {

      double time  = __TIME[i];
      double dt    = __DT[i];
      double theta = __THETA[i];
      // transformed GP
      double gx[2] = { time + dt * GX[0], time + dt * GX[1] };

      double r = 0.0;
      
      // uk in Gauss-Points
      GU[0].equ(GX[1], __primal_solution[i], GX[0], __primal_solution[i+1]);
      GU[1].equ(GX[0], __primal_solution[i], GX[1], __primal_solution[i+1]);

      // d_t uk in Gauss-Points
      GlobalVector DU = __primal_solution[i+1];
      DU.sadd(1.0/dt, -1.0/dt, __primal_solution[i]);
      
      // weights for adjoint solution in GP
      double gzw[2] = {1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * GX[0] - 1.0),  
		       1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * GX[1] - 1.0) };
      

      // higher approximation of U and dt U in GP;
      int i0 = i-i%2;
      assert(__DT[i0]==__DT[i0+1]);
      
      
      assert(i0+2 < __primal_solution.size());
      

      if (1)
	{
	  
	  double T0 = __TIME[i0];
	  double T1 = __TIME[i0+1];
	  double T2 = __TIME[i0+2];

	  assert(fabs(T1-T0-__DT[i0])<1.e-10);
	  assert(fabs(T2-T1-__DT[i0+1])<1.e-10);
	  assert(__DT[i0]==__DT[i0+1]);
	  assert(__THETA[i0]==__THETA[i0+1]);

	  GlobalVector U0 = __primal_solution[i0];
	  GlobalVector U1 = __primal_solution[i0+1];
	  GlobalVector U2 = __primal_solution[i0+2];

	  //T1 = T0+__DT[i0]*__THETA[i0]*2.0;

	  //cout << 2.0 * __THETA[i0] << endl;
	  
	  GlobalVector M1=U0;
	  M1.add( 2.0*__THETA[i0],U1);
	  M1.add(-2.0*__THETA[i0],U0);
	  

	  InterpolateQuadratic(GUH[0], U0,U1,U2, gx[0],T0,T1,T2);
	  InterpolateQuadratic(GUH[1], U0,U1,U2, gx[1],T0,T1,T2);

	  InterpolateQuadraticDerivative(GDUH[0], U0,U1,U2,gx[0], T0,T1,T2);
	  InterpolateQuadraticDerivative(GDUH[1], U0,U1,U2,gx[1], T0,T1,T2);

	  
	}
      

      
      /////// Form & Functional
      TS->GetGV(z) = __dual_solution[i+1];
      

      // GP1
      SetTimeData(gx[0], dt, theta);
      TS->Zero(f);
      TS->AdjointRhsFunctional(f, GU[0], 0.5* __DT[i]);
      
      TS->GetGV(u) = GU[0];
      TS->AddNodeVector("u",u);
      TS->AdjointFormOnly(f,z, - 0.5 * dt * gzw[0]); // - weil unten u-u_higher steht
      TS->DeleteNodeVector("u");

      TS->GetGV(u).add(-1.0, GUH[0]);
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,u);

      // GP2
      SetTimeData(gx[1], dt, theta);
      TS->Zero(f);
      TS->GetGV(u) = GU[1]; 

      TS->AdjointRhsFunctional(f, GU[1], 0.5*  __DT[i]);
      
      TS->AddNodeVector("u",u);
      TS->AdjointFormOnly(f,z, - 0.5 * dt * gzw[1]); // - weil unten u-u_higher steht
      TS->DeleteNodeVector("u");
      TS->GetGV(u).add(-1.0, GUH[1]);
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,u);


      ///////////// time-part
      // GP1
      TS->Zero(f);
      TS->GetGV(u) = DU;
      TS->GetGV(u).add(-1.0, GDUH[0]);                // k (dt uk - dt uk_high) Minus !!!
      TS->GetMassMatrix()->vmult_time(TS->GetGV(f), TS->GetGV(u), TS->GetTimePattern(),-0.5 * dt * gzw[0]);
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,z);

      // GP2
      TS->Zero(f);
      TS->GetGV(u) = DU;
      TS->GetGV(u).add(-1.0, GDUH[1]);                // k (dt uk - dt uk_high) Minus !!!
      TS->GetMassMatrix()->vmult_time(TS->GetGV(f), TS->GetGV(u), TS->GetTimePattern(),-0.5 * dt * gzw[1]);
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,z);

      res += r;
      eta[i] += r;
    }

  return res;
}




// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    Estimator for the theta-scheme




double ThetaLoop::EstimatePrimal(DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f)
{
  assert(eta.size()==_niter);
  assert(_niter+1==__primal_solution.size());
  assert(_niter+1==__dual_solution.size());
  
  ThetaSolver* TS = dynamic_cast<ThetaSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TS);
  
  // integration of Galerkin residual with 2-point Gauss rule
  double GX[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};

  // solution and adjoint solution in Gauss-Points
  vector<GlobalVector> GU; GU.resize(2);
  vector<GlobalVector> GZ; GZ.resize(2);

  double res = 0.0;
  for (int i=0;i<_niter;++i)
    {
      double time  = __TIME[i];
      double dt    = __DT[i];
      double theta = __THETA[i];
      // transformed GP
      double gx[2] = { time + dt * GX[0], time + dt * GX[1] };

      double r = 0.0;
      
      // uk in Gauss-Points
      GU[0] = __primal_solution[i]; GU[0].sadd(GX[1], GX[0], __primal_solution[i+1]);
      GU[1] = __primal_solution[i]; GU[1].sadd(GX[0], GX[1], __primal_solution[i+1]);

      
      // weights for adjoint solution in GP
      double gzw[2] = {1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * GX[0] - 1.0) ,
		       1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * GX[1] - 1.0) };
      

      //////////// approximation using big patches.
      int i0 = i - i%3;
      assert(i0+3 < __dual_solution.size());
      // InterpolateLinear(GZ[0], __dual_solution[i0], __dual_solution[i0+3], gx[0], __TIME[i0]-0.5*__DT[i0], __TIME[i0+3]-0.5*__DT[i0+3]);
      // InterpolateLinear(GZ[1], __dual_solution[i0], __dual_solution[i0+3], gx[1], __TIME[i0]-0.5*__DT[i0], __TIME[i0+3]-0.5*__DT[i0+3]);
      InterpolateLinear(GZ[0], __dual_solution[i0], __dual_solution[i0+3], gx[0], __TIME[i0], __TIME[i0+3]);
      InterpolateLinear(GZ[1], __dual_solution[i0], __dual_solution[i0+3], gx[1], __TIME[i0], __TIME[i0+3]);
      GZ[0].add(-gzw[0], __dual_solution[i+1]);
      GZ[1].add(-gzw[1], __dual_solution[i+1]);
      
      


      ///////// (dt u, z)
      // GP1
      TS->Zero(f);
      SetTimeData(gx[0], dt, theta);
      TS->GetGV(u) = __primal_solution[i+1]; TS->GetGV(u).add(-1.0,__primal_solution[i]);
      TS->GetMassMatrix()->vmult_time(TS->GetGV(f), TS->GetGV(u), TS->GetTimePattern(),  0.5);

      TS->Rhs(f,-0.5 * dt);
      TS->GetGV(u) = GU[0]; 
      TS->FormOnly(f,u, 0.5 * dt);
      TS->GetGV(z) = GZ[0];
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,z);

      // GP1
      TS->Zero(f);
      SetTimeData(gx[1], dt, theta);
      TS->GetGV(u) = __primal_solution[i+1]; TS->GetGV(u).add(-1.0,__primal_solution[i]);
      TS->GetMassMatrix()->vmult_time(TS->GetGV(f), TS->GetGV(u), TS->GetTimePattern(),  0.5);

      TS->Rhs(f,-0.5 * dt);
      TS->GetGV(u) = GU[1]; 
      TS->FormOnly(f,u, 0.5 * dt);
      TS->GetGV(z) = GZ[1];
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,z);

      res += r;
      eta[i] += r;
    }

  
  return res;
}

double ThetaLoop::EstimateDual(DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f)
{
  //  assert(_niter%2==0);
  
  assert(eta.size()==_niter);
  assert(_niter+1==__primal_solution.size());
  assert(_niter+1==__dual_solution.size());
  
  ThetaSolver* TS = dynamic_cast<ThetaSolver*>(GetMultiLevelSolver()->GetSolver());
  assert(TS);
  
  // integration of Galerkin residual with 2-point Gauss rule
  double GX[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};

  // solution and adjoint solution in Gauss-Points
  vector<GlobalVector> GU;   GU.resize(2, __primal_solution[0]);
  vector<GlobalVector> GUH;  GUH.resize(2, __primal_solution[0]);
  vector<GlobalVector> GDUH; GDUH.resize(2, __primal_solution[0]);
  vector<GlobalVector> GZ;   GZ.resize(2, __dual_solution[0]);
  
  double res = 0.0;
  for (int i=0;i<_niter;++i)
    {
      double time  = __TIME[i];
      double dt    = __DT[i];
      double theta = __THETA[i];
      // transformed GP
      double gx[2] = { time + dt * GX[0], time + dt * GX[1] };

      double r = 0.0;
      
      // uk in Gauss-Points
      GU[0].equ(GX[1], __primal_solution[i], GX[0], __primal_solution[i+1]);
      GU[1].equ(GX[0], __primal_solution[i], GX[1], __primal_solution[i+1]);

      // d_t uk in Gauss-Points
      GlobalVector DU = __primal_solution[i+1];
      DU.sadd(1.0/dt, -1.0/dt, __primal_solution[i]);
      
      // weights for adjoint solution in GP
      double gzw[2] = {1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * GX[0] - 1.0),  
		       1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * GX[1] - 1.0) };
      

      // higher approximation of U and dt U in GP;
      //      assert(_niter%6==0);
      assert(_niter%3==0);
      int i0 = (i/3)*3;
      //      int i0 = (i/6)*6;
      if (i0<3) i0+=3;
      assert(i0+3 < __primal_solution.size());
      i0 -= 3;
      //      assert(i0+6<__primal_solution.size());
      
      InterpolateQuadratic(GUH[0], __primal_solution[i0], __primal_solution[i0+3], __primal_solution[i0+6], gx[0], __TIME[i0], __TIME[i0+3], __TIME[i0+6]);
      InterpolateQuadratic(GUH[1], __primal_solution[i0], __primal_solution[i0+3], __primal_solution[i0+6], gx[1], __TIME[i0], __TIME[i0+3], __TIME[i0+6]);
      
      InterpolateQuadraticDerivative(GDUH[0], __primal_solution[i0], __primal_solution[i0+3], __primal_solution[i0+6], gx[0], __TIME[i0], __TIME[i0+3], __TIME[i0+6]);
      InterpolateQuadraticDerivative(GDUH[1], __primal_solution[i0], __primal_solution[i0+3], __primal_solution[i0+6], gx[1], __TIME[i0], __TIME[i0+3], __TIME[i0+6]);      



      
      /////// Form & Functional
      TS->GetGV(z) = __dual_solution[i+1];
      

      // GP1
      SetTimeData(gx[0], dt, theta);
      TS->Zero(f);
      TS->AdjointRhsFunctional(f, GU[0], 0.5* __DT[i]);
      
      TS->GetGV(u) = GU[0];
      TS->AddNodeVector("u",u);
      TS->AdjointFormOnly(f,z, - 0.5 * dt * gzw[0]); // - weil unten u-u_higher steht
      TS->DeleteNodeVector("u");

      TS->GetGV(u).add(-1.0, GUH[0]);
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,u);

      // GP2
      SetTimeData(gx[1], dt, theta);
      TS->Zero(f);
      TS->GetGV(u) = GU[1]; 

      TS->AdjointRhsFunctional(f, GU[1], 0.5*  __DT[i]);
      
      TS->AddNodeVector("u",u);
      TS->AdjointFormOnly(f,z, - 0.5 * dt * gzw[1]); // - weil unten u-u_higher steht
      TS->DeleteNodeVector("u");
      TS->GetGV(u).add(-1.0, GUH[1]);
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,u);


      ///////////// time-part
      // GP1
      TS->Zero(f);
      TS->GetGV(u) = DU;
      TS->GetGV(u).add(-1.0, GDUH[0]);                // k (dt uk - dt uk_high) Minus !!!
      TS->GetMassMatrix()->vmult_time(TS->GetGV(f), TS->GetGV(u), TS->GetTimePattern(),-0.5 * dt * gzw[0]);
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,z);

      // GP2
      TS->Zero(f);
      TS->GetGV(u) = DU;
      TS->GetGV(u).add(-1.0, GDUH[1]);                // k (dt uk - dt uk_high) Minus !!!
      TS->GetMassMatrix()->vmult_time(TS->GetGV(f), TS->GetGV(u), TS->GetTimePattern(),-0.5 * dt * gzw[1]);
      TS->SetBoundaryVectorZero(f);
      r += TS->ScalarProduct(f,z);
      
      res += r;
      eta[i] += r;
    }
  return res;
}




















double ThetaLoop::Estimator(double& ei, double& ep, double& ed, DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f)
{
  ei=ep=ed=0;
  
  eta.resize(_niter);
  eta.zero();

  DoubleVector etaI(_niter), etaP(_niter), etaD(_niter);
  etaI.zero(); etaP.zero(); etaD.zero();
   

  ei = EstimateInt    (etaI,u,z,f);
  ep = EstimatePrimal (etaP,u,z,f);
  ed = EstimateDual   (etaD,u,z,f);


  
  eta=etaI; eta.add(0.5,etaP); eta.add(0.5, etaD);

  assert(eta.size()%3==0);
  ofstream ETA("eta.txt");
  double s1=0,s2=0;
  
  for (int i=0;i<eta.size()/3;++i)
    {
      double ee = eta[3*i] + eta[3*i+1] + eta[3*i+2];
      ETA << i << " " <<scientific<< ee << endl;;
      //      ETA << i << " " << eta[i] << " " << etaP[i] << " " << etaD[i] << " " << etaI[i] << endl;
      s1 += ee;
      s2 += fabs(ee);
    }
  
  ETA.close();

  return ei + 0.5 * ep + 0.5 * ed;
}


double ThetaLoop::EstimatorTheta(double& ei, double& ep, double& ed, DoubleVector& eta, VectorInterface& u, VectorInterface& z, VectorInterface& f)
{
  ei=ep=ed=0;
  
  eta.resize(_niter);
  eta.zero();

  DoubleVector etaI(_niter), etaP(_niter), etaD(_niter);
  etaI.zero(); etaP.zero(); etaD.zero();


  ei = EstimateInt         (etaI,u,z,f);

  ep = EstimatePrimalTheta (etaP,u,z,f);

  ed = EstimateDualTheta   (etaD,u,z,f);

  
  eta=etaI; eta.add(0.5,etaP); eta.add(0.5, etaD);

  ofstream ETA("eta.txt");
  double s1=0,s2=0;
  ETA << "step \t estimator \t estimator \t etaI[i] \t etaP[i] \t etaD[i]" <<  endl;;
  for (int i=0;i<eta.size();++i)
    {
      double ee = eta[i];
      ETA << i << "\t " << scientific<< ee <<"\t "<<etaI[i] + 0.5 * etaP[i] + 0.5 * etaD[i]<<"\t "<<etaI[i]<<"\t "<<etaP[i]<<"\t "<<etaD[i]<<  endl;;
      //      ETA << i << " " << eta[i] << " " << etaP[i] << " " << etaD[i] << " " << etaI[i] << endl;
      s1 += ee;
      s2 += fabs(ee);
    }
  
  ETA.close();

  return ei + 0.5 * ep + 0.5 * ed;
}


// --------------------------------------------------

void ThetaLoop::run_theta(const std::string& problemlabel, int NITER, double THETA)
{
  VectorInterface u("u"), f("f"), z("z");
  
  // Programm-Parameter
  _niter   = NITER;
  __theta0 = THETA;
  double T0 = THETA;
  /////////// Gitter
  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(z);
  GetMultiLevelSolver()->ReInitVector(f);
  
  cout << "\nMesh [l,nn,nc]: ";
  cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " " << GetMeshAgent()->ncells() << endl;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();

  ///////////
  /////////// Init Zeitgitter
  ///////////
  DataFormatHandler DFH;
  DFH.insert("Tstart",&__Tstart);
  DFH.insert("Tstop",&__Tstop);
  FileScanner FS(DFH,_paramfile,"Loop");

  assert(_niter%2==0);
  __DT.clear(); __THETA.clear(); __TIME.clear();
  __TIME.push_back(__Tstart);

  double dt = (__Tstop-__Tstart)/static_cast<double>(_niter);
  THETA = T0 + SHIFT * dt;
  //  THETA = 1.0;
  for (int i=0;i<_niter;++i)
    {
      __DT.push_back(dt);
      __THETA.push_back(THETA);
      __TIME.push_back(__TIME[i]+__DT[i]);
    }
  assert(fabs(__TIME[_niter]-__Tstop)<1.e-12);
  SetTimeData(__TIME[0], __DT[0], __THETA[0]);


  

  ///////////// Vorwaerts
  PrimalLoop(u,f);
  cout<<"--------------------------------------------------------"<<endl;
  
  ///////////// Funktionale
  SetTimeData(__Tstop, __DT[0], __THETA[0]);
  assert(__primal_solution.size()==_niter+1);
  double func_time = GetTimeFunctional(f,u);
  cout.precision(10);
  cout <<"func_time:"<< func_time << endl;

  double est=0,esti=0,estp=0,estd=0;

  if (1)
    {
      cout<<"Adjoint--------------------------------------------------------"<<endl;
      //  Duale Loesung
      AdjointLoop(z,f,u);
      cout<<"Estimator--------------------------------------------------------"<<endl;

      // //  Fehler Schaetzen...
      DoubleVector eta;
      
      est = EstimatorTheta(esti,estp,estd, eta, u,z,f);
      cout.precision(10);
    }

  double error = func_time - exact_time;

  int KK =NITER;
  


  //  cerr << "$2^{-" << log(KK)/log(2) << "}\\pi$&$";
  ofstream cerr("cerr",ios::app);
  cerr.precision(10);
  
  cerr << "KK" << KK<< "\t";
  cerr.precision(2);
  cerr << "\t error: "<<scientific <<error << "\t est:"  <<est;
  
  cerr.precision(2);
  //cerr << "$&$" << fixed <<  est/error << "$&$" << scientific << estp << "$$" << estd << "$&$" << esti << "$\\\\" << endl;
  cerr << "\t" << fixed <<"\test/error: "<< scientific<<  est/error << "\t est_primal: "  << estp << "\t est_dua: " << estd << " \t est_i: " << esti << "\\\\" << endl;
  cout<<"--------------------------------------------------------"<<endl;

  ofstream out("out",ios::app);
  out.precision(10);
  
  out << "dt: "<<__DT[0] << "\t theta: " << __THETA[0] << "\t" << "\t nnodes: " << 
    GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes() << "\t func_time-exact_time:" << func_time-exact_time << "\t func_time:"<<func_time<< endl;
  out.close();
}












void ThetaLoop::run(const std::string& problemlabel, int NITER, double THETA)
{
  
  /////////////////// Alles ohne Zeitgewicht
  //exact_time = 0.6264769363;  // p = 1.2 (check, ok.)
  //exact_time = 0.631715778;   // p = 2.0 (check, ok)
  //exact_time = 0.6285842472; // p = 4.0 (check, ok)
  //exact_time = 0.6190051696; // p = 8.0 (check, ok)



  VectorInterface u("u"), f("f"), z("z");
  
  // Programm-Parameter
  _niter   = NITER;
  __theta0 = THETA;
  
  /////////// Gitter
  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(z);
  GetMultiLevelSolver()->ReInitVector(f);
  cout << "\nMesh [l,nn,nc]: ";
  cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " " << GetMeshAgent()->ncells() << endl;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();

  ///////////
  /////////// Init Zeitgitter
  ///////////
  DataFormatHandler DFH;
  DFH.insert("Tstart",&__Tstart);
  DFH.insert("Tstop",&__Tstop);
  FileScanner FS(DFH,_paramfile,"Loop");

  // Test.     abwechselnd 0.5 dt   1.5 dt  0.5 dt usw.
  assert(_niter%2==0);
  __DT.clear(); __THETA.clear(); __TIME.clear();
  __TIME.push_back(__Tstart);

  double alpha = 1.0 - 0.5 * sqrt(2.0);
  //  THETA = (1.0-2.0*alpha)/(1.0-alpha);

  for (int i=0;i<_niter;++i)

    {
      double dt = (__Tstop-__Tstart)/static_cast<double>(_niter);

      if      (i%3==0)   __DT.push_back(dt * 3.0 * alpha);
      else if (i%3==1)   __DT.push_back(dt * 3.0 * (1.0 - 2.0 * alpha));
      else if (i%3==2)   __DT.push_back(dt * 3.0 * alpha);
      else abort();



      if      (i%3==0)   __THETA.push_back(THETA);
      else if (i%3==1)   __THETA.push_back(1.0-THETA);
      else if (i%3==2)   __THETA.push_back(THETA);
      else abort();
      
      // if (i%4==0) __THETA.push_back(1.0);
      // else 
      //      __THETA.push_back(THETA);
      
      
      __TIME.push_back(__TIME[i]+__DT[i]);
    }
  assert(fabs(__TIME[_niter]-__Tstop)<1.e-12);
  SetTimeData(__TIME[0], __DT[0], __THETA[0]);



  ///////////// Vorwaerts
  PrimalLoop(u,f);
  ///////////// Funktionale
  SetTimeData(__Tstop, __DT[0], __THETA[0]);
  assert(__primal_solution.size()==_niter+1);
  double func_time = GetTimeFunctional(f,u);
  cout.precision(10);
  
  cout << func_time << endl;


  double est=0,esti=0,estp=0,estd=0;

  if (1)
    {
      //  Duale Loesung
      AdjointLoop(z,f,u);
      
      // //  Fehler Schaetzen...
      DoubleVector eta;
      
      
      est = Estimator(esti,estp,estd, eta, u,z,f);
      cout.precision(10);
      
    }

  double error = func_time - exact_time;
  
  int KK =NITER/3;
  


  cerr << KK << "&$";
  
  cerr.precision(2);
  cerr << scientific << error << "$&$" << est;
  
  cerr.precision(2);
  cerr << "$&$" << fixed <<  est/error << "$&$" << scientific << estp << "$&$" << estd << "$&$" << esti << "$\\\\" << endl;
  
  ofstream out("out",ios::app);
  out.precision(10);
  
  out << __DT[0] << "\t" << __THETA[0] << "\t" << "\t" << 
    GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes() << "\t" << func_time-exact_time << endl;
  out.close();
}















void ThetaLoop::run_adaptive(const std::string& problemlabel, int NITER, double THETA)
{
  //  exact_time =  0.07262488678 ; // p = 8.0 (check, ok)

  //  exact_time =  0.004546656091 ; // p = 8.0 (check, ok)

  VectorInterface u("u"), f("f"), z("z");
  
  // Programm-Parameter
  _niter   = NITER;
  __theta0 = THETA;
  
  /////////// Gitter
  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(z);
  GetMultiLevelSolver()->ReInitVector(f);
  cout << "\nMesh [l,nn,nc]: ";
  cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " " << GetMeshAgent()->ncells() << endl;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();

  ///////////
  /////////// Init Zeitgitter
  ///////////
  DataFormatHandler DFH;
  DFH.insert("Tstart",&__Tstart);
  DFH.insert("Tstop",&__Tstop);
  FileScanner FS(DFH,_paramfile,"Loop");

  // Test.     abwechselnd 0.5 dt   1.5 dt  0.5 dt usw.

  __DT.clear(); __THETA.clear(); __TIME.clear();
  __TIME.push_back(__Tstart);

  double alpha = 1.0 - 0.5 * sqrt(2.0);

  for (int i=0;i<_niter;++i)
    {
      double dt = (__Tstop-__Tstart)/static_cast<double>(_niter);

      if      (i%3==0)   __DT.push_back(dt * 3.0 * alpha);
      else if (i%3==1)   __DT.push_back(dt * 3.0 * (1.0 - 2.0 * alpha));
      else if (i%3==2)   __DT.push_back(dt * 3.0 * alpha);
      else abort();

      if      (i%3==0)   __THETA.push_back(THETA);
      else if (i%3==1)   __THETA.push_back(1.0-THETA);
      else if (i%3==2)   __THETA.push_back(THETA);
      else abort();
      __TIME.push_back(__TIME[i]+__DT[i]);
    }
  assert(fabs(__TIME[_niter]-__Tstop)<1.e-12);
  SetTimeData(__TIME[0], __DT[0], __THETA[0]);



  for (int itt=0;itt<10;++itt)
    {
      _niter = __TIME.size()-1;
      
      ///////////// Vorwaerts
      PrimalLoop(u,f);
      reflevel=itt+1;
      ///////////// Funktionale
      SetTimeData(__Tstop, __DT[0], __THETA[0]);
      double func_time = GetTimeFunctional(f,u);
      
      double est=0,esti=0,estp=0,estd=0;
      
      AdjointLoop(z,f,u);
      DoubleVector eta;
      est = Estimator(esti,estp,estd, eta, u,z,f);

      double error = func_time - exact_time;
      cerr << "$" << _niter <<  "$&$";
      cerr.precision(2);
      cerr << scientific << error << "$&$" << est;
      
      cerr.precision(2);
      cerr << "$&$" << fixed <<  est/error << "$&$" << scientific << estp << "$&$" << estd << "$&$" << esti << "$\\\\" << endl;
      
      
      
      
      // Adapt Time - Mesh
      DoubleVector EE;
      // assert(eta.size()%6==0);
      // EE.resize(eta.size()/6);
      assert(eta.size()%3==0);
       EE.resize(eta.size()/3);
      EE.zero();
      
      for (int i=0;i<eta.size()/3;++i)
	{
	  for (int j=0;j<3;++j)
	    EE[i]+=(fabs(eta[3*i+j]));
	}
      char sEE[20];
      sprintf(sEE,"X/eta_%i.txt", _niter);
      
      ofstream ETA(sEE);
      for (int i=0;i<EE.size(); ++i)
	{
	  ETA << __TIME[3*i] << " " << EE[i] << endl
	      << __TIME[3*i+3] << " " << EE[i] << endl << endl << endl;
	}
      ETA.close();
      

      double mean = EE.sum() / EE.size();
      nvector<double> DTO = __DT;
      assert(EE.size()*3 == __DT.size());
      assert(EE.size()*3 == __THETA.size());
      assert(EE.size()*3 == __TIME.size()-1);
      __DT.clear();
      
      for (int i=0;i<EE.size();++i)
	{
	  if (EE[i]>=2.0 * mean)
	    {
	      __DT.push_back(DTO[3*i]*0.5);
	      __DT.push_back(DTO[3*i+1]*0.5);
	      __DT.push_back(DTO[3*i+2]*0.5);
	      // __DT.push_back(DTO[6*i+3]*0.5);
	      // __DT.push_back(DTO[6*i+4]*0.5);
	      // __DT.push_back(DTO[6*i+5]*0.5);
	      __DT.push_back(DTO[3*i]*0.5);
	      __DT.push_back(DTO[3*i+1]*0.5);
	      __DT.push_back(DTO[3*i+2]*0.5);
	      // __DT.push_back(DTO[6*i+3]*0.5);
	      // __DT.push_back(DTO[6*i+4]*0.5);
	      // __DT.push_back(DTO[6*i+5]*0.5);
	    }
	  else
	    {
	      __DT.push_back(DTO[3*i]);
	      __DT.push_back(DTO[3*i+1]);
	      __DT.push_back(DTO[3*i+2]);
	      // __DT.push_back(DTO[6*i+3]);
	      // __DT.push_back(DTO[6*i+4]);
	      // __DT.push_back(DTO[6*i+5]);
	    }
	}
      __THETA.clear();
      assert(__DT.size()%3==0);
      for (int i=0;i<__DT.size()/3;++i)
	{

	  __THETA.push_back(THETA);
	  __THETA.push_back(1.0-THETA);
	  __THETA.push_back(THETA);
	}
      __TIME.clear();
      __TIME.push_back(0);
      for (int i=0;i<__DT.size();++i)
	__TIME.push_back(__DT[i] + __TIME[i]);
      cout << __TIME[__TIME.size()-1] << endl;
      

      assert(fabs(__TIME[__TIME.size()-1]-__Tstop)<1.e-10);
      

      // cerr << __DT.size() << " " << __THETA.size() << " " << __TIME.size() << " "
      // 	   << __TIME[__TIME.size()-1] << endl;
      

      char s[20];
      sprintf(s,"X/tut_%i.txt", _niter);
      
      ofstream TUT(s);
      for (int i=0;i<__DT.size()/3;++i)
	{
	  double dd = __DT[3*i] + __DT[3*i+1] + __DT[3*i+2];
	  TUT << __TIME[3*i] << " " << dd << endl
	      << __TIME[3*i+3] << " " << dd << endl << endl << endl;
	}
      TUT.close();
      
    }
  
}





















void ThetaLoop::run_adaptive_theta2(const std::string& problemlabel, int NITER, double THETA)
{
  VectorInterface u("u"), f("f"), z("z");
  
  // Programm-Parameter
  _niter   = NITER;
  __theta0 = THETA;
  double T0 = THETA;
  
  
  /////////// Gitter
  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(z);
  GetMultiLevelSolver()->ReInitVector(f);
  cout << "\nMesh [l,nn,nc]: ";
  cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " " << GetMeshAgent()->ncells() << endl;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();

  ///////////
  /////////// Init Zeitgitter
  ///////////
  DataFormatHandler DFH;
  DFH.insert("Tstart",&__Tstart);
  DFH.insert("Tstop",&__Tstop);
  FileScanner FS(DFH,_paramfile,"Loop");

  // Test.     abwechselnd 0.5 dt   1.5 dt  0.5 dt usw.

  __DT.clear(); __THETA.clear(); __TIME.clear();
  __TIME.push_back(__Tstart);

  double alpha = 1.0/3.0;

  
  
  for (int i=0;i<_niter;++i)
    {
      double dt = (__Tstop-__Tstart)/static_cast<double>(_niter);

      if      (i%3==0)   __DT.push_back(dt * 3.0 * alpha);
      else if (i%3==1)   __DT.push_back(dt * 3.0 * alpha);
      else if (i%3==2)   __DT.push_back(dt * 3.0 * alpha);
      else abort();

      if      (i%3==0)   __THETA.push_back(T0 + SHIFT * __DT[i]);
      else if (i%3==1)   __THETA.push_back(T0 + SHIFT * __DT[i]);
      else if (i%3==2)   __THETA.push_back(T0 + SHIFT * __DT[i]);
      else abort();
      __TIME.push_back(__TIME[i]+__DT[i]);
    }
  assert(fabs(__TIME[_niter]-__Tstop)<1.e-12);
  SetTimeData(__TIME[0], __DT[0], __THETA[0]);



  for (int itt=0;itt<10;++itt)
    {
      _niter = __TIME.size()-1;
      
      ///////////// Vorwaerts 
      reflevel=itt+1;
      PrimalLoop(u,f);
      ///////////// Funktionale
      SetTimeData(__Tstop, __DT[0], __THETA[0]);
      double func_time = GetTimeFunctional(f,u);
      
      double est=0,esti=0,estp=0,estd=0;
      
      AdjointLoop(z,f,u);
      DoubleVector eta;
      est = Estimator(esti,estp,estd, eta, u,z,f);

      double error = func_time - exact_time;
      cerr << "$" << _niter <<  "$&$";
      cerr.precision(2);
      cerr << scientific << error << "$&$" << est;
      
      cerr.precision(2);
      cerr << "$&$" << fixed <<  est/error << "$&$" << scientific << estp << "$&$" << estd << "$&$" << esti << "$\\\\" << endl;
      
      
      
      
      // Adapt Time - Mesh
      DoubleVector EE;
      // assert(eta.size()%6==0);
      // EE.resize(eta.size()/6);
      assert(eta.size()%3==0);
       EE.resize(eta.size()/3);
      EE.zero();
      
      for (int i=0;i<eta.size()/3;++i)
	{
	  for (int j=0;j<3;++j)
	    EE[i]+=(fabs(eta[3*i+j]));
	}
      char sEE[20];
      sprintf(sEE,"X/eta_%i.txt", _niter);
      
      ofstream ETA(sEE);
      for (int i=0;i<EE.size(); ++i)
	{
	  ETA << __TIME[3*i] << " " << EE[i] << endl
	      << __TIME[3*i+3] << " " << EE[i] << endl << endl << endl;
	}
      ETA.close();
      

      double mean = EE.sum() / EE.size();
      nvector<double> DTO = __DT;
      assert(EE.size()*3 == __DT.size());
      assert(EE.size()*3 == __THETA.size());
      assert(EE.size()*3 == __TIME.size()-1);
      __DT.clear();
      
      for (int i=0;i<EE.size();++i)
	{
	  if (EE[i]>=mean)
	    {
	      __DT.push_back(DTO[3*i]*0.5);
	      __DT.push_back(DTO[3*i+1]*0.5);
	      __DT.push_back(DTO[3*i+2]*0.5);
	      __DT.push_back(DTO[3*i]*0.5);
	      __DT.push_back(DTO[3*i+1]*0.5);
	      __DT.push_back(DTO[3*i+2]*0.5);
	    }
	  else
	    {
	      __DT.push_back(DTO[3*i]);
	      __DT.push_back(DTO[3*i+1]);
	      __DT.push_back(DTO[3*i+2]);
	    }
	}
      __THETA.clear();
      assert(__DT.size()%3==0);
      for (int i=0;i<__DT.size();++i)
	__THETA.push_back(T0 + SHIFT * __DT[i]);
      __TIME.clear();
      __TIME.push_back(0);
      for (int i=0;i<__DT.size();++i)
	__TIME.push_back(__DT[i] + __TIME[i]);
      cout << __TIME[__TIME.size()-1] << endl;
      

      assert(fabs(__TIME[__TIME.size()-1]-__Tstop)<1.e-10);
      

      // cerr << __DT.size() << " " << __THETA.size() << " " << __TIME.size() << " "
      // 	   << __TIME[__TIME.size()-1] << endl;
      

      char s[20];
      sprintf(s,"X/tut_%i.txt", _niter);
      
      ofstream TUT(s);
      for (int i=0;i<__DT.size()/3;++i)
	{
	  double dd = __DT[3*i] + __DT[3*i+1] + __DT[3*i+2];
	  TUT << __TIME[3*i] << " " << dd << endl
	      << __TIME[3*i+3] << " " << dd << endl << endl << endl;
	}
      TUT.close();
      
    }
  
  
}















void ThetaLoop::run_adaptive_theta(const std::string& problemlabel, int NITER, double THETA)
{
  // Programm-Parameter
  _niter   = NITER;
  VectorInterface u("u"), f("f"), z("z");
  double T0 = THETA;
  
  
  /////////// Gitter
  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(z);
  GetMultiLevelSolver()->ReInitVector(f);
  cout << "\nMesh [l,nn,nc]: ";
  cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " " << GetMeshAgent()->ncells() << endl;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();

  ///////////
  /////////// Init Zeitgitter
  ///////////
  DataFormatHandler DFH;
  DFH.insert("Tstart",&__Tstart);
  DFH.insert("Tstop",&__Tstop);
  FileScanner FS(DFH,_paramfile,"Loop");


  __DT.clear(); __THETA.clear(); __TIME.clear();
  __TIME.push_back(__Tstart);

  double dt = (__Tstop-__Tstart)/static_cast<double>(_niter);
  
  THETA = T0+ SHIFT * dt;
  
  for (int i=0;i<_niter;++i)
    {
      __DT.push_back(dt);
      __THETA.push_back(T0 + SHIFT * dt);
      __TIME.push_back(__TIME[i]+__DT[i]);
    }

  assert(fabs(__TIME[_niter]-__Tstop)<1.e-12);
  SetTimeData(__TIME[0], __DT[0], __THETA[0]);

  for (int itt=0;itt<30;++itt)
    {
      _niter = __TIME.size()-1;
      
      ///////////// Vorwaerts
      reflevel=itt+1;
      PrimalLoop(u,f);
      
      ///////////// Funktionale
      SetTimeData(__Tstop, __DT[0], __THETA[0]);
      double func_time = GetTimeFunctional(f,u);
      
      double est=0,esti=0,estp=0,estd=0;
      
      AdjointLoop(z,f,u);
      DoubleVector eta;
      est = EstimatorTheta(esti,estp,estd, eta, u,z,f);

      double error = func_time - exact_time;
      ofstream cerr("cerr",ios::app);
      cerr.precision(10);
      cerr << "" << _niter <<  "$\t$";

      cerr.precision(2);
      cerr << "\t error: "<<scientific <<error << "\t est:"  <<est;
      
      cerr.precision(2);
      //cerr << "$&$" << fixed <<  est/error << "$&$" << scientific << estp << "$$" << estd << "$&$" << esti << "$\\\\" << endl;
      cerr << "\t" << fixed <<"\test/error: "<< scientific<<  est/error << "\t est_primal: "  << estp << "\t est_dua: " << estd << " \t est_i: " << esti << "\\\\" << endl;
      cout<<"--------------------------------------------------------"<<endl;

      
      // Adapt Time - Mesh
      DoubleVector EE;
      // assert(eta.size()%6==0);
      // EE.resize(eta.size()/6);
      assert(eta.size()%4==0);
      EE.resize(eta.size()/4);
      EE.zero();
      
      for (int i=0;i<EE.size();++i)
	{
	  EE[i]=0.0;
	  for (int j=0;j<4;++j)
	    EE[i] += fabs(eta[4*i+j]);
	}
      
	  
      char sEE[20];
      sprintf(sEE,"X/eta_%i.txt", _niter);
      

      
      ofstream ETA(sEE);
      for (int i=0;i<EE.size(); ++i)
	{
	  ETA << __TIME[4*i] << " " << EE[i] << endl
	      << __TIME[4*i+4] << " " << EE[i] << endl << endl << endl;
	}
      ETA.close();
      

      /////// Verfeinern
      double mean = EE.sum() / EE.size();
      nvector<double> DTO = __DT;
      assert(4*EE.size() == __DT.size());
      assert(4*EE.size() == __THETA.size());
      assert(4*EE.size() == __TIME.size()-1);
      __DT.clear();
      
      double FAC= 1.0/0.8;
      int ii=0;
      do { FAC *= 0.8;
	for (ii=0;ii<EE.size();++ii)
	  if (EE[ii]> FAC * mean) break;
      }  while (ii==EE.size());
      
      
      
      for (int i=0;i<EE.size();++i)
	{
	  double dt =  DTO[4*i];
	  assert(dt== DTO[4*i+1]);
	  assert(dt== DTO[4*i+2]);
	  assert(dt== DTO[4*i+3]);
	  if (EE[i]>=FAC * mean)
	    {
	      __DT.push_back(dt*0.5);
	      __DT.push_back(dt*0.5);
	      __DT.push_back(dt*0.5);
	      __DT.push_back(dt*0.5);
	      __DT.push_back(dt*0.5);
	      __DT.push_back(dt*0.5);
	      __DT.push_back(dt*0.5);
	      __DT.push_back(dt*0.5);
	    }
	  else 
	    {
	      __DT.push_back(dt);
	      __DT.push_back(dt);
	      __DT.push_back(dt);
	      __DT.push_back(dt);
	    }
	}
      __THETA.clear();
      for (int i=0;i<__DT.size();++i)
	__THETA.push_back(T0 + SHIFT * __DT[i]); //THETA);
      __TIME.clear();
      __TIME.push_back(0);
      for (int i=0;i<__DT.size();++i)
	__TIME.push_back(__DT[i] + __TIME[i]);
      cout << __TIME[__TIME.size()-1] << endl;
      

      assert(fabs(__TIME[__TIME.size()-1]-__Tstop)<1.e-10);
      
      char s[20];
      sprintf(s,"X/tut_%i.txt", _niter);
      
      ofstream TUT(s);
      for (int i=0;i<__DT.size()/4;++i)
	{
	  double dd = __DT[4*i];
	  TUT << __TIME[4*i] << " " << dd << endl
	      << __TIME[4*i+4] << " " << dd << endl << endl << endl;
	}
      TUT.close();
      
    }
  
  
}

void ThetaLoop::Extra(vector<double>& f)
{
  if (f.size()<3) return;
  double a2 = f[f.size()-1];
  double a1 = f[f.size()-2];
  double a0 = f[f.size()-3];
  double A = (-a1*a1+a0*a2)/ (a0+a2-2.0*a1);
  double q = -log((a1-a2)/(a0-a1)) / log(2.0);
  
  exact_time = A;
  
  cerr << "EXTRA: " << A << "\t [" << q << "]";
  
}



void  ThetaLoop::run_exact(const std::string& problemlabel, int NITER, double THETA)
{
  VectorInterface u("u"), f("f"), z("z");

  vector<double>  func_time;

  for (int II=0;II<8;++II)
    {
      __theta0 = THETA;
      
      _niter = NITER * pow(2.0,II);
      cout << "Solve with " << _niter << " steps" << endl;
      
      

      /////////// Gitter
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      GetMultiLevelSolver()->ReInit(problemlabel);
      GetMultiLevelSolver()->ReInitVector(u);
      GetMultiLevelSolver()->ReInitVector(z);
      GetMultiLevelSolver()->ReInitVector(f);
      cout << "\nMesh [l,nn,nc]: ";
      cout << GetMeshAgent()->nlevels() << " " << GetMeshAgent()->nnodes() << " " << GetMeshAgent()->ncells() << endl;
      GetMultiLevelSolver()->GetSolver()->OutputSettings();
      
      ///////////
      /////////// Init Zeitgitter
      ///////////
      DataFormatHandler DFH;
      DFH.insert("Tstart",&__Tstart);
      DFH.insert("Tstop",&__Tstop);
      FileScanner FS(DFH,_paramfile,"Loop");

      // aquidistant
      __DT.clear(); __THETA.clear(); __TIME.clear();
      __TIME.push_back(__Tstart);
      double dt = (__Tstop-__Tstart)/static_cast<double>(_niter);

      for (int i=0;i<_niter;++i)
	{
	  __DT.push_back(dt);
	  __THETA.push_back(__theta0);
	  __TIME.push_back(__TIME[i]+__DT[i]);
	}
      assert(fabs(__TIME[_niter]-__Tstop)<1.e-12);
      SetTimeData(__TIME[0], __DT[0], __THETA[0]);


      ///////////// Vorwaerts
      PrimalLoop(u,f);

      ///////////// Funktionale
      SetTimeData(__Tstop, __DT[0], __THETA[0]);
      assert(__primal_solution.size()==_niter+1);
      func_time.push_back(GetTimeFunctional(f,u));

      cerr.precision(10);
      

      cerr << "TIME: " << func_time[func_time.size()-1] << "\t" ;
      Extra(func_time);
      cerr << endl;
      
      
    }
}
