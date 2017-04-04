#include "loop.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "functionals.h"


using namespace std;
using namespace Gascoigne;

double __DT,__TIME, __THETA;
double __DT_OLD, __THETA_OLD;

bool __FIRSTSTEP,__MIDDLESTEP,__LASTSTEP;

bool __average_functional=false;


nvector<double> exakt;

double PQ0(double t,double t0,double t1,double t2)  { return (t-t1)*(t-t2)/(t0-t1)/(t0-t2); }
double PQ1(double t,double t0,double t1,double t2)  { return (t-t0)*(t-t2)/(t1-t0)/(t1-t2); }
double PQ2(double t,double t0,double t1,double t2)  { return (t-t0)*(t-t1)/(t2-t0)/(t2-t1); }


double dtPQ0(double t,double t0,double t1,double t2) { return ((t-t1)+(t-t2))/(t0-t1)/(t0-t2); } 
double dtPQ1(double t,double t0,double t1,double t2) { return ((t-t0)+(t-t2))/(t1-t0)/(t1-t2); } 
double dtPQ2(double t,double t0,double t1,double t2) { return ((t-t0)+(t-t1))/(t2-t0)/(t2-t1); } 

double QUAD(double t, double t0, double t1, double t2, double x0, double x1, double x2)
{
  return x0*PQ0(t,t0,t1,t2) + x1*PQ1(t,t0,t1,t2) + x2*PQ2(t,t0,t1,t2);
}
double LIN(double t, double t0, double t1, double x0, double x1)
{
  return (t-t1)/(t0-t1) * x0 + (t-t0)/(t1-t0) * x1;
}


double dtQUAD(double t, double t0, double t1, double t2, double x0, double x1, double x2)
{
  return x0*dtPQ0(t,t0,t1,t2) + x1*dtPQ1(t,t0,t1,t2) + x2*dtPQ2(t,t0,t1,t2);
}
double dtLIN(double t, double t0, double t1, double x0, double x1)
{
  return 1.0/(t0-t1) * x0 + 1.0/(t1-t0) * x1;
}


void Loop::InitPrimal(GlobalVector& U, int m, int q)
{
  assert(m<__dt.size());
  assert((q>=0)&&(q<2));
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};
  
  for (int i=0;i<U.n();++i)
    {
      U(i,0) = __U[m+1](i,0);
      for (int c=1;c<3;++c)
	U(i,c) = (1.0-chi[q]) * __U[m](i,c) + chi[q] * __U[m+1](i,c);
    }
  
  GetMultiLevelSolver()->GetSolver()->GetDiscretization()->HNAverage(U);
}

void Loop::InitPrimalTime(GlobalVector& dtU, int m)
{
  assert(m<__dt.size());

  double dt = __dt[m];
  
  for (int i=0;i<dtU.n();++i)
    {
      dtU(i,0) = 0.0;
      for (int c=1;c<3;++c)
	dtU(i,c) = (__U[m+1](i,c)-__U[m](i,c))/dt;
    }
  
  GetMultiLevelSolver()->GetSolver()->GetDiscretization()->HNAverage(dtU);
}


void Loop::InitDual(GlobalVector& Z, int m, int q)
{
  assert(m<__dt.size());
  assert((q>=0)&&(q<2));
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};
  double theta = __theta[m];
  double gzw[2] = {1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * chi[0] - 1.0) ,
		   1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * chi[1] - 1.0) };
  
  for (int i=0;i<Z.n();++i)
    {
      for (int c=0;c<3;++c)
	Z(i,c) = gzw[q] * __Z[m+1](i,c);
    }
  GetMultiLevelSolver()->GetSolver()->GetDiscretization()->HNAverage(Z);
}

void Loop::InitDualHigher(GlobalVector& Z, int m, int q)
{
  assert(m<__dt.size());
  assert((q>=0)&&(q<2));
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};
  double theta = __theta[m];
  double time  = __time[m] + __dt[m] * chi[q];
  double gzw[2] = {1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * chi[0] - 1.0) ,
		   1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * chi[1] - 1.0) };

  
  // quadratic??
  int m0 = m-m%2;
  for (int i=0;i<Z.n();++i)
    for (int c=0;c<3;++c)
      Z(i,c) = QUAD(time, __time[m0],__time[m0+1],__time[m0+2],
		    __Z[m0](i,c),__Z[m0+1](i,c),__Z[m0+2](i,c))
	- __Z[m+1](i,c)*gzw[q];
  GetMultiLevelSolver()->GetSolver()->GetDiscretization()->HNAverage(Z);
}

void Loop::InitDualHigherFST(GlobalVector& Z, int m, int q)
{
  assert(m<__dt.size());
  assert((q>=0)&&(q<2));
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};
  double theta = __theta[m];
  double time  = __time[m] + __dt[m] * chi[q];
  double gzw[2] = {1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * chi[0] - 1.0) ,
		   1.0 + 3.0 * (2.0 * theta-1.0) * (2.0 * chi[1] - 1.0) };

  
  // quadratic??
  int m0 = m-m%6;
  for (int i=0;i<Z.n();++i)
    for (int c=0;c<3;++c)
      Z(i,c) = QUAD(time, __time[m0],__time[m0+3],__time[m0+6],
		    __Z[m0](i,c),__Z[m0+3](i,c),__Z[m0+6](i,c))
	- __Z[m+1](i,c)*gzw[q];
  GetMultiLevelSolver()->GetSolver()->GetDiscretization()->HNAverage(Z);
}

void Loop::InitPrimalHigher(GlobalVector& U, int m, int q)
{
  assert(m<__dt.size());
  assert((q>=0)&&(q<2));
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};
  double dt    = __dt[m];
  double time  = __time[m] + dt*chi[q];
  
  int m0 = m-m%2; // start of macro-time step, interpolate quadratic in m0->m0+2
  assert(m0+2 < __U.size());


  for (int i=0;i<U.n();++i)
    {
      // quadratic interpolation on macro-step for velocities
      U(i,0) = QUAD(time,        __time[m0],    __time[m0+1], __time[m0+2],
		    __U[m0](i,0),__U[m0+1](i,0),__U[m0+2](i,0))
	- __U[m+1](i,0);
      for (int c=1;c<3;++c)
	U(i,c) = QUAD(time,        __time[m0],    __time[m0+1], __time[m0+2],
		      __U[m0](i,c),__U[m0+1](i,c),__U[m0+2](i,c))
	  - LIN(time, __time[m],__time[m+1], __U[m](i,c), __U[m+1](i,c));
    }
}


void Loop::InitPrimalHigherFST(GlobalVector& U, int m, int q)
{
  int N = __dt.size();
  assert(N%3==0);
  
  assert(m<__dt.size());
  assert((q>=0)&&(q<2));
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};
  double dt    = __dt[m];
  double time  = __time[m] + dt*chi[q];
  
  int m0 = m-m%6; // start of macro-time step, interpolate quadratic in m0, m0+3, m0+6
  assert(m0+6 < __U.size());


  for (int i=0;i<U.n();++i)
    {
      // quadratic interpolation on macro-step for velocities
      U(i,0) = QUAD(time,        __time[m0],    __time[m0+3], __time[m0+6],
		    __U[m0](i,0),__U[m0+3](i,0),__U[m0+6](i,0))
	- __U[m+1](i,0);
      for (int c=1;c<3;++c)
	U(i,c) = QUAD(time,        __time[m0],    __time[m0+3], __time[m0+6],
		      __U[m0](i,c),__U[m0+3](i,c),__U[m0+6](i,c))
	  - LIN(time, __time[m],__time[m+1], __U[m](i,c), __U[m+1](i,c));
    }
}


//  am Zeitpunkt t=T fuer die rechte Seite
void Loop::InitPrimalHigherEnd(GlobalVector& U)
{
  int M = __dt.size();
  assert(__time.size()==M+1);
  
  assert(__U.size()==M+1);
  
  for (int i=0;i<U.n();++i)
    {
      // first, linear interpolation of pressure unknown at right end of interval

      // Interpolation durch Mitte.
      U(i,0) = 0;// .5*(__U[M](i,0)-__U[M-1](i,0));
      U(i,1)=0;
      U(i,2)=0;
    }

}


void Loop::InitPrimalTimeHigher(GlobalVector& U, int m, int q)
{
  assert(m<__dt.size());
  assert((q>=0)&&(q<2));
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};
  double dt    = __dt[m];
  double time  = __time[m] + dt*chi[q];
  
  int m0 = m-m%2; // start of macro-time step, interpolate quadratic in m0->m0+2
  assert(m0+2 < __U.size());
  

  for (int i=0;i<U.n();++i)
    {
      // quadratic interpolation on macro-step for velocities
      U(i,0) = 0.0;
      for (int c=1;c<3;++c)
	U(i,c) = dtQUAD(time, __time[m0], __time[m0+1], __time[m0+2], __U[m0](i,c),__U[m0+1](i,c),__U[m0+2](i,c)) 
	  - dtLIN(time, __time[m],__time[m+1], __U[m](i,c), __U[m+1](i,c));
    }

}
void Loop::InitPrimalTimeHigherFST(GlobalVector& U, int m, int q)
{
  assert(m<__dt.size());
  assert((q>=0)&&(q<2));
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};
  double dt    = __dt[m];
  double time  = __time[m] + dt*chi[q];
  
  int m0 = m-m%6; // start of macro-time step, interpolate quadratic in m0->m0+2
  assert(m0+6 < __U.size());
  

  for (int i=0;i<U.n();++i)
    {
      // quadratic interpolation on macro-step for velocities
      U(i,0) = 0.0;
      for (int c=1;c<3;++c)
	U(i,c) = dtQUAD(time, __time[m0], __time[m0+3], __time[m0+6], __U[m0](i,c),__U[m0+3](i,c),__U[m0+6](i,c)) 
	  - dtLIN(time, __time[m],__time[m+1], __U[m](i,c), __U[m+1](i,c));
    }

}


// computed with two point Gauss
nvector<double> Loop::GetMeanFunctional()
{
  VectorInterface tmp1("tmp1"), tmp2("tmp2");
  GetMultiLevelSolver()->ReInitVector(tmp1);
  GetMultiLevelSolver()->ReInitVector(tmp2);
  
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};

  nvector<double> sum_j;
  assert(__U.size()==_niter+1);
  GlobalVector U = __U[0];
  for (int i=0;i<_niter;++i) // time interval [i,i+1]
    {
      double dt   = __dt[i];
      for (int q=0;q<2;++q) //two gauss points
	{
	  __TIME = (1.0-chi[q]) * __time[i] + chi[q] * __time[i+1];
	  InitPrimal(U,i,q); // discrete solution in Gauss Point
	  GetMultiLevelSolver()->GetSolver()->GetGV(tmp1) = U;
	  nvector<double> j = ComputeFunctionals(tmp2,tmp1);
	  if (sum_j.size()!=j.size()) sum_j.resize(j.size());
	  sum_j.add(0.5*dt,j);
	}
    }
  GetMultiLevelSolver()->DeleteVector(tmp1);
  GetMultiLevelSolver()->DeleteVector(tmp2);
  return sum_j;
}

// computed with two point Gauss
void Loop::WriteFunctionals()
{
  ostringstream strstr;
  double dt = __dt[0];
  strstr << "functionals_" << dt << ".txt";
  ofstream OUT(strstr.str().c_str());
  
  VectorInterface tmp1("tmp1"), tmp2("tmp2");
  GetMultiLevelSolver()->ReInitVector(tmp1);
  GetMultiLevelSolver()->ReInitVector(tmp2);
  
  nvector<double> sum_j;
  assert(__U.size()==_niter+1);
  GlobalVector U = __U[0];
  for (int i=0;i<=_niter;++i) // time interval [i,i+1]
    {
      double time=__time[i];
      GetMultiLevelSolver()->GetSolver()->GetGV(tmp1) = __U[i];
      nvector<double> j = ComputeFunctionals(tmp2,tmp1);
      OUT << time << "\t" << j << endl;
    }
  GetMultiLevelSolver()->DeleteVector(tmp1);
  GetMultiLevelSolver()->DeleteVector(tmp2);
  OUT.close();
}

nvector<double> Loop::GetEndFunctional(const GlobalVector& U)
{
  VectorInterface tmp1("tmp1"), tmp2("tmp2");
  GetMultiLevelSolver()->ReInitVector(tmp1);
  GetMultiLevelSolver()->ReInitVector(tmp2);
  
  // set time to last step
  __DT = __dt[__dt.size()-1];
  __THETA = __theta[__theta.size()-1];
  
  
  GetMultiLevelSolver()->GetSolver()->GetGV(tmp1) = U;
  nvector<double> j = ComputeFunctionals(tmp2,tmp1);
  GetMultiLevelSolver()->DeleteVector(tmp1);
  GetMultiLevelSolver()->DeleteVector(tmp2);

  return j;
}


void Loop::InitUniformTimeMesh(double start, double stop, double dt,double theta)
{
  //  cerr << "Init uniform time mesh from " << start << " to " << stop 
  //       << "\t" << dt << "/" << theta << endl;
  
  assert(stop>start);
  assert(dt>0);
  assert(theta>0);
  
  _niter = static_cast<int> ( (stop-start+1.e-12)/dt );

  __theta.resize(_niter);
  __dt.resize   (_niter);
  __start.resize(_niter);
  __stop.resize (_niter);
  __time.resize(_niter+1);
  
  for (int i=0;i<_niter;++i)
    {
      __theta[i] = theta;
      __dt[i]    = dt;
      __start[i]   = static_cast<double> (i)*dt;
      __time[i]   = static_cast<double> (i)*dt;
      __stop[i]    = static_cast<double> (i+1)*dt;
    }
  __time[_niter]   = static_cast<double> (_niter)*dt;
  
  assert(__stop[_niter-1]==stop);
  assert(__start[0]==start);
  __U.resize(_niter+1);
  __Z.resize(_niter+1);
}


void Loop::AdaptFSTMesh(nvector<double>& fst)
{
  double gamma = 1.0;
  do
    {
      int M = fst.size();
      assert(6*M==__start.size());
      
      double avg = fst.norm_l1()/static_cast<double> (M);
      nvector<double> odt = __dt;
      nvector<double> otime = __time;
      
      cerr << "refine at: " << 0.1*avg << " < " << avg << " < " << gamma*avg << endl;
      
      double alpha = 1.0 - 1.0/sqrt(2.0);
      double theta = 0.55;
      double oldlast = __time[__time.size()-1];
      
      __time.clear();
      __dt.clear();
      __theta.clear();
      __start.clear();
      __stop.clear();
      __time.push_back(otime[0]);
      for (int i=0;i<M;++i)
	{
	  assert(6*i+6<otime.size());
	  double olddt=otime[6*i+3]-otime[6*i];

	  bool coarse=false;
	  // coarse????
	  if (0)
	  if (i+1<M)
	    if (fabs(fst[i])+fabs(fst[i+1])<0.1*avg)
	      {
		double otc = otime[6*i+12]-otime[6*i];
		coarse=true;
		__dt.push_back(0.5*otc*alpha);            __theta.push_back(theta);
		__dt.push_back(0.5*otc*(1.0-2.0*alpha));  __theta.push_back(1.0-theta);
		__dt.push_back(0.5*otc*alpha);            __theta.push_back(theta);
		__dt.push_back(0.5*otc*alpha);            __theta.push_back(theta);
		__dt.push_back(0.5*otc*(1.0-2.0*alpha));  __theta.push_back(1.0-theta);
		__dt.push_back(0.5*otc*alpha);            __theta.push_back(theta);
		++i;
	      }
	  if (coarse) continue;
	  if (fabs(fst[i])<gamma*avg) // don't refine
	    {
	      __dt.push_back(olddt*alpha);            __theta.push_back(theta);
	      __dt.push_back(olddt*(1.0-2.0*alpha));  __theta.push_back(1.0-theta);
	      __dt.push_back(olddt*alpha);            __theta.push_back(theta);
	      __dt.push_back(olddt*alpha);            __theta.push_back(theta);
	      __dt.push_back(olddt*(1.0-2.0*alpha));  __theta.push_back(1.0-theta);
	      __dt.push_back(olddt*alpha);            __theta.push_back(theta);
	    }
	  else
	    {
	      __dt.push_back(0.5*olddt*alpha);            __theta.push_back(theta);
	      __dt.push_back(0.5*olddt*(1.0-2.0*alpha));  __theta.push_back(1.0-theta);
	      __dt.push_back(0.5*olddt*alpha);            __theta.push_back(theta);
	      
	      __dt.push_back(0.5*olddt*alpha);            __theta.push_back(theta);
	      __dt.push_back(0.5*olddt*(1.0-2.0*alpha));  __theta.push_back(1.0-theta);
	      __dt.push_back(0.5*olddt*alpha);            __theta.push_back(theta);
	      
	      __dt.push_back(0.5*olddt*alpha);            __theta.push_back(theta);
	      __dt.push_back(0.5*olddt*(1.0-2.0*alpha));  __theta.push_back(1.0-theta);
	      __dt.push_back(0.5*olddt*alpha);            __theta.push_back(theta);
	      
	      __dt.push_back(0.5*olddt*alpha);            __theta.push_back(theta);
	      __dt.push_back(0.5*olddt*(1.0-2.0*alpha));  __theta.push_back(1.0-theta);
	      __dt.push_back(0.5*olddt*alpha);            __theta.push_back(theta);
	    }
	}
      double t = __time[0];
      for (int i=0;i<__dt.size();++i)
	{
	  __start.push_back(t);
	  t+=__dt[i];
	  __time.push_back(t);
	  __stop.push_back(t);
	}
      
      _niter = __dt.size();
      assert(__time.size()==_niter+1);
      assert(fabs(__time[_niter]-oldlast)<1.e-10);
      __time[_niter] = oldlast;
      __stop[_niter-1] = oldlast;
      
      __U.resize(_niter+1);
      __Z.resize(_niter+1);
      gamma*=0.8;
    }
  while (_niter==fst.size()*6);
    
}
  
  
void Loop::InitUniformFSTMesh(double start, double stop, double dt)
{
  //  cerr << "Init uniform time mesh from " << start << " to " << stop 
  //       << "\t" << dt << "/" << theta << endl;
  
  assert(stop>start);
  assert(dt>0);
  
  _niter = 3*static_cast<int> ( (stop-start+1.e-12)/dt ); // dt is step length of macro time step
  double alpha = 1.0 - 1.0/sqrt(2.0);
  double theta = 0.55;

  __theta.resize(_niter);
  __dt.resize   (_niter);
  __start.resize(_niter);
  __stop.resize (_niter);
  __time.resize(_niter+1);
  
  __start[0]=start;
  for (int i=0;i<_niter;++i)
    {
      __time[i]=__start[i];
      if ((i%3==0)||(i%3==2))
	{
	  __theta[i] = theta;
	  __dt[i]    = alpha*dt;
	  __stop[i]=__start[i]+alpha*dt;
	}
      else
	{
	  __theta[i] = 1.0-theta;
	  __dt[i]    = (1.0-2.0*alpha)*dt;
	  __stop[i]=__start[i]+(1.0-2.0*alpha)*dt;
	}

      if ((i+1)<_niter) __start[i+1] = __stop[i];
    }
  assert(fabs(stop-__stop[_niter-1])<1.e-10);
  __stop[_niter-1]=stop;
  __time[_niter] = stop;
  assert(__stop[_niter-1]==stop);
  assert(__start[0]==start);
  __U.resize(_niter+1);
  __Z.resize(_niter+1);
}

string Loop::Solve(VectorInterface& u, VectorInterface& f, string name)
{
  GetMultiLevelSolver()->GetSolver()->Zero(f);
  GetMultiLevelSolver()->GetSolver()->Rhs(f,__DT);

  GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
  string status = GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());
  Output(u,name);
  return status;
}



void Loop::PrimalLoop(VectorInterface& u, VectorInterface& f, VectorInterface& old)
{
  std::cerr << "Primal: ";
  
  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  ofstream OUT("out");

  for (int i=0;i<_niter;++i)
    {
      _iter = i+1;
      
      __DT    = __dt[i];
      __THETA = __theta[i];

//       if (i%(_niter/10)==0)
// 	{
// 	  std::cerr << "\t" << i/(_niter/10)*10 << "%";
// 	}
       std::cerr << i << "/" << _niter << "\t" << __TIME+__DT
       		<< "\t" << __dt[i] << "\t theta: " << __theta[i] << endl;
           
      // old
      __TIME  = __start[i];      
      GetMultiLevelSolver()->GetSolver()->GetGV(old) = __U[i];
      GetMultiLevelSolver()->SolutionTransfer(old);
      GetMultiLevelSolver()->GetSolver()->Zero(f);
      GetMultiLevelSolver()->GetSolver()->Rhs(f,(1.0-__THETA) * __DT);
      
      // new
      __TIME  = __stop[i];      
      GetMultiLevelSolver()->AddNodeVector("old",old);
      GetMultiLevelSolver()->GetSolver()->Rhs(f,__THETA * __DT);

      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
      
      // Solve
      GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());
      //assert(GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo())=="converged");
      Output(u,"Results/u");
  
      // Store solution
      __U[i+1] = GetMultiLevelSolver()->GetSolver()->GetGV(u);
      GetMultiLevelSolver()->DeleteNodeVector("old");
    }
  std::cerr << endl;
  
  OUT.close();
}




void Loop::AdjointLoop(VectorInterface& z, VectorInterface& uu, VectorInterface& f, VectorInterface& zold)
{
  std::cerr << "Adjoint: " ;
  
  ofstream OUT("out");
  for (_iter=_niter;_iter>=0;--_iter)
    {
      // INIT
      __FIRSTSTEP=__LASTSTEP=__MIDDLESTEP=false;
      if (_iter==_niter) { __FIRSTSTEP=true;  }
      else if (_iter==0) { __LASTSTEP=true; }
      else __MIDDLESTEP=true;
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      if (_iter%(_niter/10)==0) std::cerr << "\t" <<  _iter/(_niter/10)*10 << "%";

      
      if (_iter>0) {__DT = __dt[_iter-1]; __THETA = __theta[_iter-1]; __DT_OLD=-1; __THETA_OLD = -1;}
      if (_iter<_niter) { __DT_OLD = __dt[_iter]; __THETA_OLD = __theta[_iter]; }

      // old adjoint solution Z[i+1]
      if (_iter<_niter)
	{
	  GetMultiLevelSolver()->GetSolver()->GetGV(zold) = __Z[_iter+1];
	  GetMultiLevelSolver()->SolutionTransfer(zold);
	}
      GetMultiLevelSolver()->AddNodeVector("zold",zold);
      
      // init primal solution U[i]
      __TIME = __time[_iter];
      GetMultiLevelSolver()->GetSolver()->GetGV(uu) = __U[_iter];
      GetMultiLevelSolver()->SolutionTransfer(uu);
      GetMultiLevelSolver()->AddNodeVector("uu",uu);

      std::cerr << _iter << "/" << _niter << "\t" << __TIME+__DT
       		<< "\t" << __DT << "\t theta: " << __THETA << endl;
      /////////// RHS
      GetMultiLevelSolver()->GetSolver()->Zero(f);

      if (__average_functional)
	{
	  if (__FIRSTSTEP||__MIDDLESTEP)
	    {
	      assert(_iter<__U.size());
	      GetMultiLevelSolver()->GetSolver()->GetDiscretization()->AddNodeVector("U",&__U[_iter]);
	      GetMultiLevelSolver()->GetSolver()->Rhs(f,0.5*__dt[_iter-1]);
	      GetMultiLevelSolver()->GetSolver()->GetDiscretization()->DeleteNodeVector("U");
	    }
	  if (__LASTSTEP||__MIDDLESTEP)
	    {
	      assert(_iter<__dt.size());
	      GetMultiLevelSolver()->GetSolver()->GetDiscretization()->AddNodeVector("U",&__U[_iter]);
	      GetMultiLevelSolver()->GetSolver()->Rhs(f,0.5*__dt[_iter]);
	      GetMultiLevelSolver()->GetSolver()->GetDiscretization()->DeleteNodeVector("U");
	    }
	}
      else // Functional at t=T
	{
	  if (__FIRSTSTEP)
	    {
	      GetMultiLevelSolver()->GetSolver()->GetDiscretization()->AddNodeVector("U",&__U[__U.size()-1]);
	      GetMultiLevelSolver()->GetSolver()->Rhs(f);
	      GetMultiLevelSolver()->GetSolver()->GetDiscretization()->DeleteNodeVector("U");
	    }
	  
	}
      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
      GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(z);

      
      /////////////// SOLVE (own function due to rhs in first step only
      //assert(GetMultiLevelSolver()->Solve(z,f,GetSolverInfos()->GetNLInfo())=="converged");
      GetMultiLevelSolver()->Solve(z,f,GetSolverInfos()->GetNLInfo());
      Output(z,"Results/z");
      ///////////////

      // save Z[i]
      __Z[_iter] = GetMultiLevelSolver()->GetSolver()->GetGV(z);
      
      GetMultiLevelSolver()->DeleteNodeVector("uu");
      GetMultiLevelSolver()->DeleteNodeVector("zold");
    }
  
  std::cerr << endl;
  
  OUT.close();
}











////////////////////////////////////////////////// ESTIMATE



void Loop::EstimateQuadrature(nvector<double>& estA, VectorInterface& u, VectorInterface& f)
{
  int M = __dt.size();
  assert(M==__theta.size());assert(M+1==__time.size());
  assert(M+1==__U.size());assert(M+1==__Z.size());
  estA.resize(M);
  estA.zero();

  // calculate rho(u_h)(z_h) of Galerkin formulation in two Gauss-Points
  
  GlobalVector U=__U[0],dtU=__U[0], Z=__U[0];
  
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};

  GalerkinResidual GR(_paramfile);
  
  const StdSolver* S = dynamic_cast<const StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(S);
  
  for (int m=0;m<M;++m)
    {
      // set time interval
      __DT = __dt[m]; __THETA = __theta[m];
      
      // two quadrature points
      for (int q=0;q<2;++q)
	{
	  __TIME = (1.0-chi[q]) * __time[m] + chi[q] * __time[m+1];
	  InitPrimal(U,m,q);     // linear interpolation of primal solution
	  InitPrimalTime(dtU,m);
	  InitDual(Z,m,q);
	  
	  S->GetDiscretization()->AddNodeVector("dtU",&dtU);
	  S->GetDiscretization()->AddNodeVector("Z",&Z);
	  estA[m] += __DT  * 0.5 * S->GetDiscretization()->ComputeDomainFunctional(U,GR);
	  S->GetDiscretization()->DeleteNodeVector("dtU");
	  S->GetDiscretization()->DeleteNodeVector("Z");
	}
    }
  // Startpunkt + (u^0_k,z^0_k) fliegt raus, da u^0_k=0.
  InitialResidual IR;
  S->GetDiscretization()->AddNodeVector("Z",&__Z[0]);
  estA[0] += S->GetDiscretization()->ComputeDomainFunctional(__U[0],IR); 
  S->GetDiscretization()->DeleteNodeVector("Z"); 
}


void Loop::EstimatePrimal(nvector<double>& estP, VectorInterface& u, VectorInterface& f)
{
  int M = __dt.size();
  assert(M==__theta.size());assert(M+1==__time.size());
  assert(M+1==__U.size());assert(M+1==__Z.size());
  estP.resize(M);
  estP.zero();

  // calculate rho(u_h)(z_h) of Galerkin formulation in two Gauss-Points
  
  GlobalVector U=__U[0],dtU=__U[0], ZH=__U[0];
  
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};

  GalerkinResidual GR(_paramfile);
  
  const StdSolver* S = dynamic_cast<const StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(S);
  
  for (int m=0;m<M;++m)
    {
      // set time interval
      __DT = __dt[m]; __THETA = __theta[m];
      
      // two quadrature points
      for (int q=0;q<2;++q)
	{
	  __TIME = (1.0-chi[q]) * __time[m] + chi[q] * __time[m+1];
	  InitPrimal(U,m,q);     // linear interpolation of primal solution
	  InitPrimalTime(dtU,m);
	  InitDualHigher(ZH,m,q);
	  
	  S->GetDiscretization()->AddNodeVector("dtU",&dtU);
	  S->GetDiscretization()->AddNodeVector("Z",&ZH);
	  estP[m] += __DT  * 0.5 * S->GetDiscretization()->ComputeDomainFunctional(U,GR);
	  S->GetDiscretization()->DeleteNodeVector("dtU");
	  S->GetDiscretization()->DeleteNodeVector("Z");
	}
    }
  // 0 bei meiner Interpolation von Zhigher

}


void Loop::EstimateDual(nvector<double>& estD, VectorInterface& u, VectorInterface& f)
{
  int M = __dt.size();
  assert(M==__theta.size());assert(M+1==__time.size());
  assert(M+1==__U.size());assert(M+1==__Z.size());
  estD.resize(M);
  estD.zero();

  GlobalVector U=__U[0],UH = __U[0], dtUH=__U[0], Z=__U[0];
  
  double chi[2] = {0.5-sqrt(1.0/12.0), 0.5+sqrt(1.0/12.0)};

  GalerkinAdjointResidual GAR(_paramfile);

  
  
  const StdSolver* S = dynamic_cast<const StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(S);
  
  for (int m=0;m<M;++m)
    {
      // set time interval
      __DT = __dt[m]; __THETA = __theta[m];
      
      // two quadrature points
      for (int q=0;q<2;++q)
	{
	  __TIME = (1.0-chi[q]) * __time[m] + chi[q] * __time[m+1];
	  InitPrimal(U,m,q);     // linear interpolation of primal solution
	  InitPrimalHigher(UH,m,q);
	  InitPrimalTimeHigher(dtUH,m,q);     
	  InitDual(Z,m,q);
	  
	  S->GetDiscretization()->AddNodeVector("UH",  &UH);
	  S->GetDiscretization()->AddNodeVector("dtUH",&dtUH);
	  S->GetDiscretization()->AddNodeVector("Z",   &Z);
	  estD[m] += __DT  * 0.5 * S->GetDiscretization()->ComputeDomainFunctional(U,GAR);
	  S->GetDiscretization()->DeleteNodeVector("UH");
	  S->GetDiscretization()->DeleteNodeVector("dtUH");
	  S->GetDiscretization()->DeleteNodeVector("Z");
	}
    }

  if (__average_functional) /////// Mean-Functional
    {
      Blift lift(_paramfile);
      set<int> cols;cols.insert(80);
      
      for (int m=0;m<M;++m)
	{
	  __DT = __dt[m]; __THETA = __theta[m];
	  // two quadrature points
	  for (int q=0;q<2;++q)
	    {
	      __TIME = (1.0-chi[q]) * __time[m] + chi[q] * __time[m+1];
	      InitPrimalHigher(UH,m,q);
	      estD[m] += __DT  * 0.5 * S->GetDiscretization()->ComputeBoundaryFunctional(UH,cols,lift);
	    }
	}
      
    }
  

  // Funktional zum Endzeitpunkt taucht nicht auf, da PrimalHigher=0 in t=T
}


////////////////////////////////////////////////// ESTIMATE






void Loop::run(const std::string& problemlabel)
{
  double start,stop,theta,dt;
  
  DataFormatHandler DFH;
  DFH.insert("start_time" ,    &start , 0.0);
  DFH.insert("stop_time" ,    &stop , 0.0);
  DFH.insert("theta" ,    &theta , 0.0);
  DFH.insert("dt" ,    &dt,0.0);
  FileScanner FS(DFH, _paramfile, "Equation");

  




  for (int IT=0;IT<10;++IT)
    {
      cerr << endl;
      
      cerr << "=============== " << IT << "\t" << dt << std::endl;
      
      //      theta = 0.5+2.0*dt;
      
      //      InitUniformTimeMesh(start,stop,dt,theta);
      InitUniformFSTMesh(start,stop,dt);
      
      VectorInterface z("z"), u("u"), f("f"), old("old");
      GetMultiLevelSolver()->ReInit(problemlabel);
      GetMultiLevelSolver()->ReInitVector(z);
      GetMultiLevelSolver()->ReInitVector(u);
      GetMultiLevelSolver()->ReInitVector(old);
      GetMultiLevelSolver()->ReInitVector(f);
      InitSolution(u);
      __U[0] = GetMultiLevelSolver()->GetSolver()->GetGV(u);
      
      //      cerr << endl << "****************** Primal ****************** " << endl;
      PrimalLoop(u,f,old);
      
      ///////// FUNC
      nvector<double> j = GetEndFunctional(__U[__U.size()-1]);
      //nvector<double> j = GetMeanFunctional();
      
      //      cerr << endl << "****************** Adjoint ****************** "  << endl;
      GetMultiLevelSolver()->ReInit("adjoint");
      AdjointLoop(z,u,f,old);
      
      cerr << endl << "****************** Estimate ****************** " << endl;
      nvector<double> estA, estP, estD;
      double error = exakt[0] - j[0];
      //      cerr << "Error: " << error << endl;
      

      EstimateQuadrature(estA, u,f);
      EstimatePrimal(estP, u,f);
      EstimateDual(estD, u,f);
      double ea = estA.sum();
      double ep = 0.5*estP.sum();
      double ed = 0.5*estD.sum();
      double est = ea+ep+ed;
      
      
      assert(estA.size()==estP.size());
      assert(estA.size()==estD.size());
      assert(__time.size()==estD.size()+1);


      ostringstream stream;
      stream << "est_" << dt << ".txt";
      ofstream EST(stream.str().c_str());
      EST.precision(20);
      for (int i=0;i<estA.size();++i)
      	EST << __time[i+1] << "\t"
      	    << estA[i] << " " << estP[i] << " " << estD[i] << " " <<  estA[i]+0.5*estP[i]+0.5*estD[i] << endl;
      EST.close();
      
      
      cerr << "Estimate: " << est << "\t";
      cerr << ea << "\t";
      cerr << ep << "\t";
      cerr << ed << endl;
      
      
      cerr << "Functionals: " <<  j << endl;
      cerr << "Errors: ";
      for (int i=0;i<j.size();++i)
	cerr << error << "\t";
      cerr << endl;

      cerr << "\033[31m" << "\033[1m" 
	   << "Efficiency: " << est / error
	   << "\033[0m" << endl;
      
      dt *=0.5;
    }
  
  
}



void Loop::run_fst_adaptive(const std::string& problemlabel)
{
  double start,stop,theta,dt;
  
  DataFormatHandler DFH;
  DFH.insert("start_time" ,    &start , 0.0);
  DFH.insert("stop_time" ,    &stop , 0.0);
  DFH.insert("theta" ,    &theta , 0.0);
  DFH.insert("dt" ,    &dt,0.0);
  FileScanner FS(DFH, _paramfile, "Equation");

  InitUniformFSTMesh(start,stop,dt);
  
  for (int IT=0;IT<20;++IT)
    {
      cerr << endl;
      
      cerr << "=============== " << IT << "\t" << dt << "\t M=" << __dt.size()/3 << std::endl;
      
      VectorInterface z("z"), u("u"), f("f"), old("old");
      GetMultiLevelSolver()->ReInit(problemlabel);
      GetMultiLevelSolver()->ReInitVector(z);
      GetMultiLevelSolver()->ReInitVector(u);
      GetMultiLevelSolver()->ReInitVector(old);
      GetMultiLevelSolver()->ReInitVector(f);
      InitSolution(u);
      __U[0] = GetMultiLevelSolver()->GetSolver()->GetGV(u);
      
      ///////// SOLVE
      PrimalLoop(u,f,old);
      
      ///////// FUNC
      nvector<double> j;
      if (__average_functional) j =  GetMeanFunctional();
      else j = GetEndFunctional(__U[__U.size()-1]);

      ///////// ADJOINT
      GetMultiLevelSolver()->ReInit("adjoint");
      AdjointLoop(z,u,f,old);
      
      ///////// ESTIMATE
      cerr << endl << "****************** Estimate ****************** " << endl;
      nvector<double> estA, estP, estD;
      double error = exakt[0] - j[0];
      EstimateQuadrature(estA, u,f);
      EstimatePrimal(estP, u,f);
      EstimateDual(estD, u,f);
      double ea = estA.sum();
      double ep = 0.5*estP.sum();
      double ed = 0.5*estD.sum();
      double est = ea+ep+ed;
      assert(estA.size()==estP.size());
      assert(estA.size()==estD.size());
      assert(__time.size()==estD.size()+1);


      ///////// combine for fst - patch!
      int N = estA.size();assert(N%6==0);
      int M = N/6;
      nvector<double> fstA(M),fstP(M),fstD(M),fst(M);
      for (int i=0;i<M;++i)
	{
	  fstA[i]=0; for (int c=0;c<6;++c) fstA[i]+=estA[6*i+c];
	  fstP[i]=0; for (int c=0;c<6;++c) fstP[i]+=estP[6*i+c];
	  fstD[i]=0; for (int c=0;c<6;++c) fstD[i]+=estD[6*i+c];
	  fst[i]=fstA[i]+0.5*fstP[i]+0.5*fstD[i];
	}
      /////////// OUTPUT
      if (1)
	{
	  ostringstream stream;
	  stream << "est_" << IT  << ".txt";
	  ofstream EST(stream.str().c_str());
	  EST.precision(20);
	  for (int i=0;i<M;++i)
	    {
	      assert(6*i+6<N+1);
	      EST << __time[6*i] << "\t" << __time[6*i+6]-__time[6*i] << "\t" 
		  << fstA[i] << " " << fstP[i] << " " << fstD[i] << " " <<  fst[i] << endl; 
	      EST << __time[6*i+6] << "\t"<< __time[6*i+6]-__time[6*i] << "\t" 
		  << fstA[i] << " " << fstP[i] << " " << fstD[i] << " " <<  fst[i] << endl; 
	      EST << endl << endl;
	    }
	  EST.close();
	}      
      


      cerr << "Estimate: " << est << "\t";
      cerr << ea << "\t";
      cerr << ep << "\t";
      cerr << ed << endl;
      
      
      cerr << "Functionals: " <<  j << endl;
      cerr << "Errors: ";
      for (int i=0;i<j.size();++i)
	cerr << error << "\t";
      cerr << endl;

      cerr << "\033[31m" << "\033[1m" 
	   << "Efficiency: " << est / error
	   << "\033[0m" << endl;


      // ADAPT
      

      ///
      for (int i=0;i<fst.size();++i)
       	fst[i]=2.0*fabs(fstA[i])+fabs(fstP[i])+fabs(fstD[i]);
      AdaptFSTMesh(fst);
      
    }
  
  
}

void Loop::extrapolate(const vector<nvector<double> >& j)
{
  cerr.precision(20);
  vector<double> exact;
  if (j.size()<3) 
    {
      int nc = j[0].size();
      for (int c=0;c<nc;++c)
	cerr << j[j.size()-1][c] << "\t";
      cerr << endl;
    }
  else
    {
      int nc = j[0].size();
      for (int c=0;c<nc;++c)
	{
	  double a0 = j[j.size()-3][c];
	  double a1 = j[j.size()-2][c];
	  double a2 = j[j.size()-1][c];
	  double q = -log((a1-a2)/(a0-a1))/log(2.0);
	  double a = (a0*a2-a1*a1)/(a0+a2-2.0*a1);
	  cerr << "functional " << c << endl;
	  for (int ii=0;ii<j.size();++ii)
	    cerr << j[ii][c] << endl;
	  cerr << "extrapolate: " << a << " / " << a-a2 << " (" << q << ")\t"; 
	  cerr << endl;
	}
      cerr << endl;
    }
}


void Loop::run_exact(const std::string& problemlabel)
{
  double start,stop;
   
  double dt = 0.1;  
  DataFormatHandler DFH;
  DFH.insert("start_time" ,    &start , 0.0);
  DFH.insert("stop_time" ,    &stop , 0.0);
  DFH.insert("dt" ,    &dt , 0);
  FileScanner FS(DFH, _paramfile, "Equation");

  double theta = 0.5;
  
  
  
  
  

  
  
  vector<nvector<double> > jmean,jend;
  
  do
    {
      theta = 0.5+dt;
      //      theta = 1.0;
      
  
      InitUniformTimeMesh(start,stop,dt,theta);
      
      VectorInterface z("z"), u("u"), f("f"), old("old");
      GetMultiLevelSolver()->ReInit(problemlabel);
      GetMultiLevelSolver()->ReInitVector(z);
      GetMultiLevelSolver()->ReInitVector(u);
      GetMultiLevelSolver()->ReInitVector(old);
      GetMultiLevelSolver()->ReInitVector(f);
      InitSolution(u);
      __U[0] = GetMultiLevelSolver()->GetSolver()->GetGV(u);
      
      PrimalLoop(u,f,old);
      WriteFunctionals();
      jmean.push_back(GetMeanFunctional());
      jend.push_back(GetEndFunctional(__U[__U.size()-1]));

      cerr << "MEAN: " << dt <<  "\t" << endl;
      
      extrapolate(jmean);
      cerr << "END:  " << dt <<  "\t" << endl;
      extrapolate(jend);
      dt *= 0.5;
      
      cerr << endl;
      
      if (dt<1.e-6) break;
    }
  while (1);
}


void Loop::run_exact_fst(const std::string& problemlabel)
{
  double start,stop;
   
  double dt = 0.1;  
  DataFormatHandler DFH;
  DFH.insert("start_time" ,    &start , 0.0);
  DFH.insert("stop_time" ,    &stop , 0.0);
  DFH.insert("dt" ,    &dt , 0);
  FileScanner FS(DFH, _paramfile, "Equation");

  vector<nvector<double> > jmean,jend;
  
  do
    {
  
      InitUniformFSTMesh(start,stop,dt);
      
      VectorInterface z("z"), u("u"), f("f"), old("old");
      GetMultiLevelSolver()->ReInit(problemlabel);
      GetMultiLevelSolver()->ReInitVector(z);
      GetMultiLevelSolver()->ReInitVector(u);
      GetMultiLevelSolver()->ReInitVector(old);
      GetMultiLevelSolver()->ReInitVector(f);
      InitSolution(u);
      __U[0] = GetMultiLevelSolver()->GetSolver()->GetGV(u);
      
      PrimalLoop(u,f,old);
      WriteFunctionals();
      jmean.push_back(GetMeanFunctional());
      jend.push_back(GetEndFunctional(__U[__U.size()-1]));

      cerr << "MEAN: " << dt <<  "\t" << endl;
      
      extrapolate(jmean);
      cerr << "END:  " << dt <<  "\t" << endl;
      extrapolate(jend);
      dt *= 0.5;
      
      cerr << endl;
      
      if (dt<1.e-6) break;
    }
  while (1);
}
