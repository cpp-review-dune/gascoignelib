#include "aleloop.h"
#include "stdsolver.h"
#include "stdmultilevelsolver.h"
#include "dwrq1q2.h"
#include "dwr.h"
#include  "malteadaptor.h"
#include <sys/time.h>
#include <sys/resource.h> 
#include "dwr_boundary.h"
#include <map>

using namespace std;


extern bool __ADJOINT;
extern bool __ADJOINT1;
extern double __EXACT;
extern bool __ESTIMATE;

namespace Gascoigne
{




  void AleLoop::SolveDualProblem(const std::string& duallabel,
				 VectorInterface& z, VectorInterface& u, VectorInterface& f)
  {
    VectorInterface w("w");
    GetMultiLevelSolver()->ReInitVector(w);
    GetMultiLevelSolver()->AddNodeVector("U",u);    

    __ADJOINT = true;
    
    
    GetMultiLevelSolver()->Zero(z);
    GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(z);

    GetSolverInfos()->GetLInfo().statistics().reset();
    GetSolverInfos()->GetLInfo().control().reset();

    for (int ii=0;ii<2;++ii)
      {
	// residual
	GetMultiLevelSolver()->Zero(f);
	GetMultiLevelSolver()->GetSolver()->Rhs(f);

	GetMultiLevelSolver()->GetSolver()->Form(f,z,-1.0);
	
	GetMultiLevelSolver()->GetSolver()->SetBoundaryVectorZero(f);
	GetMultiLevelSolver()->Zero(w);

	GetSolverInfos()->GetLInfo().reset();	
	dynamic_cast<StdMultiLevelSolver*>
	  (GetMultiLevelSolver())->NewtonLinearSolve(w,f,GetSolverInfos()->GetLInfo());
	
	// update
	GetMultiLevelSolver()->GetSolver()->Add(z,1.0,w);
      }
    


    // // Delete dual 'velocity' in solid domain
    // const AleSolver* ALES = dynamic_cast<const AleSolver*> (GetMultiLevelSolver()->GetSolver());
    // const HASHSET<int>& IN = dynamic_cast<const AleSolver*> (GetMultiLevelSolver()->GetSolver())->GetInterfaceNodes();
    // const std::vector<int>&  SL2G = ALES->GetSolidL2G();
    // for (int i=0;i<SL2G.size();++i)
    //   if (IN.find(SL2G[i])==IN.end())
    // 	for (int c=0;c<2;++c)
    // 	  GetMultiLevelSolver()->GetSolver()->GetGV(z)(SL2G[i],c+1+2)=0;
    
    

    cout.width(4);
    cout.setf(ios::scientific);
    cout.precision(4);
    // cout << "Linear solve: [" 
    // 	 << GetSolverInfos()->GetLInfo().control().firstresidual() << " -> "
    // 	 << GetSolverInfos()->GetLInfo().control().residual() << "]  {"
    // 	 << GetSolverInfos()->GetLInfo().statistics().totaliter() << "}" << endl;
    
   
    
    string fname = "Results/" + duallabel;
    GetMultiLevelSolver()->GetSolver()->Visu(fname,z,_iter);    

    __ADJOINT = false;

    GetMultiLevelSolver()->DeleteNodeVector("U");
    GetMultiLevelSolver()->DeleteVector(w);
  }
  


	
  double AleLoop::Adjoint(const std::string& problemlabel,
			  const std::string& label, DoubleVector& eta,
			  const FunctionalContainer* FC,
			  VectorInterface& z, VectorInterface& u, VectorInterface& f)
  {
    assert(FC);
    assert(FC->find(label)!=FC->end());
    
    
    GetMultiLevelSolver()->SetProblem(label);
    SolveDualProblem(label, z,u,f);
    
    Dwr dwr(*GetMultiLevelSolver()->GetSolver(), 
	    GetMultiLevelSolver()->GetProblemContainer()->GetProblem(problemlabel),
	    GetMultiLevelSolver()->GetProblemContainer()->GetProblem(label));
    double estimate = dwr.Estimator(eta,f,u,z);
    GetMultiLevelSolver()->SetProblem(problemlabel);
    return estimate;
  }
	





  
  void process_mem_usage(double& vm_usage, double& resident_set)
  {
    using std::ios_base;
    using std::ifstream;
    using std::string;
    
    vm_usage     = 0.0;
    resident_set = 0.0;
    
    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat",ios_base::in);
    
    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;
    
    // the two fields we want
    //
    unsigned long vsize;
    long rss;
    
    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
		>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
		>> utime >> stime >> cutime >> cstime >> priority >> nice
		>> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
    
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage     = vsize / 1024.0;
    resident_set = rss * page_size_kb;
  }
  
  extern StopWatch S1;
  extern StopWatch S2;
  extern StopWatch S25;
  extern StopWatch S3;
  extern StopWatch S4;
  extern StopWatch S5;
  extern StopWatch S6;
  extern StopWatch S_LIN;
  StopWatch SDUAL;
  
  
  void AleLoop::AdaptMesh(const DoubleVector& eta)
  {
    if(_refiner=="etainterface") 
      {
	IntVector refnodes, coarsenodes,dummynodes;
	
	MalteAdaptor A(_paramfile,eta);
	A.refine(refnodes,coarsenodes);

	const HASHSET<int>& IN = dynamic_cast<const AleSolver*> (GetMultiLevelSolver()->GetSolver())->GetInterfaceNodes();
	  
	for (HASHSET<int>::const_iterator it = IN.begin();it!=IN.end();++it)
	  refnodes.push_back(*it);
	GetMeshAgent()->refine_nodes(refnodes);
      }
    else if(_refiner=="mean") 
      {
	IntVector refnodes;
	
	double mean = 0.0;
	mean = eta.norm_l1()/eta.size();
	
	for (int i=0;i<eta.size();++i)
	  if (eta[i]>12.0*mean) refnodes.push_back(i);
	if (refnodes.size()==0)
	  for (int i=0;i<eta.size();++i)
	    if (eta[i]>6.0*mean) refnodes.push_back(i);
	if (refnodes.size()==0)
	  for (int i=0;i<eta.size();++i)
	    if (eta[i]>2.0*mean) refnodes.push_back(i);
	if (refnodes.size()==0)
	  for (int i=0;i<eta.size();++i)
	    if (eta[i]>mean) refnodes.push_back(i);

	GetMeshAgent()->refine_nodes(refnodes);
      }
    else if(_refiner=="interface") 
      {
	IntVector refnodes, coarsenodes,dummynodes;
	
	const HASHSET<int>& IN = dynamic_cast<const AleSolver*> (GetMultiLevelSolver()->GetSolver())->GetInterfaceNodes();
	  
	for (HASHSET<int>::const_iterator it = IN.begin();it!=IN.end();++it)
	  refnodes.push_back(*it);

	const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*>
	  (GetMultiLevelSolver()->GetSolver()->GetMesh());
	assert(M);
	for (int i=0;i<M->nnodes();++i)
	  if (fabs(M->vertex2d(i).x())==6.0) refnodes.push_back(i);
		

	GetMeshAgent()->refine_nodes(refnodes);
      }
    else StdLoop::AdaptMesh(eta);
  }

  
  /*-------------------------------------------------*/

  void AleLoop::run(const std::string& problemlabel)
  {
    VectorInterface u("u"), f("f"), z("z"), help("help");
    GlobalVector    ualt;

    string fun_name = "log_" + _estimator + ".txt";
    string est_name = "est_" + _estimator + ".txt";
    

    map<string,vector<double> > func_history;
    

    for (_iter=1; _iter<=_niter; _iter++)
      {
	__ADJOINT  = false;
	__ADJOINT1 = false;
	__ESTIMATE = false;
	
	
	
	cout << "\n================== " << _iter << " ================";
	PrintMeshInformation();

	GetMultiLevelSolver()->ReInit(problemlabel);
	GetMultiLevelSolver()->ReInitVector(u);
	GetMultiLevelSolver()->ReInitVector(f);
	GetMultiLevelSolver()->ReInitVector(z);
	GetMultiLevelSolver()->ReInitVector(help);
	GetMultiLevelSolver()->InterpolateSolution(u,ualt);
	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
	
	// statistics... nr. nodes
	int nnodes = GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes();
	
	if (_iter==1) 
	  {
	    GetMultiLevelSolver()->GetSolver()->OutputSettings();
	    InitSolution(u);
	  }

	S1.reset();S2.reset();S3.reset();S4.reset();S5.reset();S6.reset();
	S25.reset();
	
	S_LIN.reset();
	
	
	S1.start();
	string solve_name = "Results/u_" + _estimator;
	Solve(u,f,solve_name);
	S1.stop();

	DoubleVector eta;
	////////////////////////////////////////////////// FUNCTIONALS
	// compute adjoint matrix
	for (int l=0;l<GetMultiLevelSolver()->nlevels();++l)
	  {
	    GetMultiLevelSolver()->GetSolver(l)->MatrixZero();
	    dynamic_cast<AleSolver*>
	      (GetMultiLevelSolver()->GetSolver(l))->AssembleDualMatrix(u,1.0);
	  }
	
	const FunctionalContainer* FC = GetMultiLevelSolver()->GetFunctionalContainer();
	assert(FC);

	for (map<string,const Functional*>::const_iterator it = FC->begin();it!=FC->end();++it)
	  {
	    const string name = it->first;

	    std::cout << endl << "Functional: " << name << endl;
	    
	    double value = GetMultiLevelSolver()->ComputeFunctional(f,u,name);
	    double error = it->second->ExactValue() - value;

	    double extra = 0.0;
	    double order = 0.0;
	    // Extrapolate?
	    func_history[name].push_back(value);	    
	    if (func_history[name].size()>2)
	      {
		double a2 = func_history[name][func_history[name].size()-1];
		double a1 = func_history[name][func_history[name].size()-2];
		double a0 = func_history[name][func_history[name].size()-3];
		extra =  (-a1*a1+a0*a2)/(a0+a2-2*a1);
		order =  -log((a1-a2)/(a0-a1))/log(2);
	      }

	    
	    cout << "\t value:       " << value << endl
		 << "\t error:       " << error << endl
		 << "\t extrapolate: " << extra << "\t [" << order << "]" << endl
		 << "\t err (extra): " << extra - value << endl;
	    
	    DoubleVector e1;
	    double estimate = Adjoint(problemlabel,name,e1, FC,z,u,f);
	    if (name==_estimator) 
	      {
		eta = e1;
		string etaname = "Results/eta_"+name;
		EtaVisu(etaname,_iter,eta);
		
	      }
	    
	    cout << "\t estimate:    " << estimate << endl
		 << "\t eff:         " << estimate / error << endl
		 << "\t eff (extra): " << estimate / (extra-value) << endl
		 << endl;

	    string fname = name + ".txt";
	    ofstream OUT(fname.c_str(), ios::app);
	    OUT << GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes()
		<< "\t" << value << " " << error << " " << extra << " " 
		<< estimate << endl;
	    OUT.close();
	  }

	// extrapolate
    
    
	if (_iter<_niter) 
	  {
	    CopyVector(ualt,u);
	    AdaptMesh(eta);
	  }
      }
    
  }



}
