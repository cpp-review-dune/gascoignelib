#include "loop.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "dwrq1q2.h"


using namespace std;


double __DT,__TIME, __THETA;





namespace Gascoigne
{

  
  
  template<int DIM>
  void Loop<DIM>::run(const std::string& problemlabel)
  {
    double STOP_TIME;
    int _initial_refine=0;
    

    if (1)
      {
	DataFormatHandler DFH;
	DFH.insert("start_time" ,    &__TIME , 0.0);
	DFH.insert("stop_time" ,    &STOP_TIME , 0.0);
	DFH.insert("theta" ,    &__THETA , 0.0);
	DFH.insert("dt" ,    &__DT,0.0);
	FileScanner FS(DFH, _paramfile, "Equation");
	assert(STOP_TIME>__TIME);
	assert(__DT>0);
	assert(__THETA>0);
      }
    if (1)
      {
	DataFormatHandler DFH;
	DFH.insert("initialrefine", &_initial_refine,0);
	FileScanner FS(DFH, _paramfile, "Loop");
      }
    


    _niter = static_cast<int> ( (STOP_TIME-__TIME+1.e-12)/__DT );

  
  

    VectorInterface u("u"), f("f"), old("old");
  
    GetMultiLevelSolver()->ReInit(problemlabel);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(old);
    GetMultiLevelSolver()->ReInitVector(f);
    InitSolution(u);

    for (int i=0;i<_initial_refine;++i)
      {	
	cout << "Mesh with " << GetMultiLevelSolver()->GetSolver()->GetGV(u).n() << " nodes" << endl;

	GlobalVector ualt;
	CopyVector(ualt,u);

	nvector<int> ref;

	// refine at chi
	// const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*> (GetMultiLevelSolver()->GetSolver()->GetMesh());
	// assert(M);
	// Chi chi;
	// for (int j=0;j<M->nnodes();++j)
	//   if (chi(M->vertex3d(j))==0) ref.push_back(j);

	// refine with primal estimator?
	DoubleVector eta;
	DwrQ1Q2 dwr(*(GetMultiLevelSolver()->GetSolver()));
	GetMultiLevelSolver()->GetSolver()->HNAverage(u);
	GetMultiLevelSolver()->AddNodeVector("old",old);
	GetMultiLevelSolver()->GetSolver()->HNAverage(old);
	dwr.EstimatorEnergy(eta,f,u);
	GetMultiLevelSolver()->DeleteNodeVector("old");
	assert(eta.size()==GetMultiLevelSolver()->GetSolver()->GetGV(u).n());
	double avg = 1.0/eta.size() * eta.norm_l1();
	for (int j=0;j<eta.size();++j)
	  if (fabs(eta[j])>avg*2.0) ref.push_back(j);
	
	
	
	GetMeshAgent()->refine_nodes(ref);

	//	GetMeshAgent()->global_refine(1);
	
	GetMultiLevelSolver()->ReInit(problemlabel);
	GetMultiLevelSolver()->ReInitVector(u);
	GetMultiLevelSolver()->ReInitVector(old);
	GetMultiLevelSolver()->ReInitVector(f);
	GetMultiLevelSolver()->InterpolateSolution(u,ualt);
      }
    cout << "Mesh with " << GetMultiLevelSolver()->GetSolver()->GetGV(u).n() << " nodes" << endl;
	

    

    ofstream func_log("functional.txt");
    func_log.precision(12);
  

    // damp at 1, 2, a.s.o.
    int DI = static_cast<int> (1.0 / __DT);
    cout << DI << endl;
    //_iteration, _solve, _matrix, _form, _ilu, _smooth;
    

    StopWatch _all;

    GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,0);    

    for (_iter=1; _iter<=_niter; _iter++)
      {

	
      
	// // Rannacher-Time-Stepping

	// if ( (_iter%DI==1) || (_iter%DI==2) )
	// 	{
	// 	  __THETA=1.0;
	// 	}
	// else __THETA = thetasave;
      
	cout << "========== " << _iter << ": " << __TIME << " -> " 
	     << __TIME+__DT << "  (" << __THETA << ")" << endl;
	__TIME += __DT;

	GetMultiLevelSolver()->Equ(old,1.0,u);

	GetMultiLevelSolver()->AddNodeVector("old",old);

	//	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

	assert(Solve(u,f)=="converged");

	if (_iter%100==0)
	  {
	    GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,_iter);
	    WriteMeshAndSolution("Results/u",u);
	  }

      
	DoubleVector juh = Functionals(u,f);
	GetMultiLevelSolver()->DeleteNodeVector("old");




	double teff = __TIME;
	int Ti=0;
	while (teff>5.0){ ++Ti; teff -=5.0; }
	
	// transition by 0.1
	double vnew = 1.0+0.1*Ti;
	double vold = 1.0+0.1*(Ti-1.0);
	if (vold<1.0) vold = 1.0;
	double sc = 1.0;
	if (teff<0.1)
	  sc = 0.5*(1.0-cos(M_PI*teff/0.1));
	//	vold = vnew = 1.4;
	
	double veff = vold + sc*(vnew-vold);

	func_log <<__TIME << "\t" << veff << "\t" << juh << endl;

	if (_iter%1==1)
	  {	
	    cout << "Mesh with " << GetMultiLevelSolver()->GetSolver()->GetGV(u).n() << " nodes" << endl;
	    
	    GlobalVector ualt;
	    CopyVector(ualt,u);
	    
	    nvector<int> ref;
	    
	    DoubleVector eta;
	    DwrQ1Q2 dwr(*(GetMultiLevelSolver()->GetSolver()));
	    GetMultiLevelSolver()->GetSolver()->HNAverage(u);
	    GetMultiLevelSolver()->AddNodeVector("old",old);
	    GetMultiLevelSolver()->GetSolver()->HNAverage(old);
	    dwr.EstimatorEnergy(eta,f,u);
	    GetMultiLevelSolver()->DeleteNodeVector("old");
	    assert(eta.size()==GetMultiLevelSolver()->GetSolver()->GetGV(u).n());
	    double avg = 1.0/eta.size() * eta.norm_l1();
	    for (int j=0;j<eta.size();++j)
	      if (fabs(eta[j])>avg*2.0) ref.push_back(j);
	    
	    
	    GetMeshAgent()->refine_nodes(ref);
	    
	    GetMultiLevelSolver()->ReInit(problemlabel);
	    GetMultiLevelSolver()->ReInitVector(u);
	    GetMultiLevelSolver()->ReInitVector(old);
	    GetMultiLevelSolver()->ReInitVector(f);
	    GetMultiLevelSolver()->InterpolateSolution(u,ualt);
	  }
	cout << "Mesh with " << GetMultiLevelSolver()->GetSolver()->GetGV(u).n() << " nodes" << endl;
      }
    
    func_log.close();

  }
  

  template class Loop<2>;
  template class Loop<3>;


}





