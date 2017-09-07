#include "loop.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include  "compose_name.h"
using namespace std;


double __DT,__TIME;


//extern bool BOUNDARY;



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
	DFH.insert("dt" ,    &__DT,0.0);
	FileScanner FS(DFH, _paramfile, "Equation");
	assert(STOP_TIME>__TIME);
	assert(__DT>0);
      }
    if (1)
      {
	DataFormatHandler DFH;
	DFH.insert("initialrefine", &_initial_refine,0);
	FileScanner FS(DFH, _paramfile, "Loop");
      }



    _niter = static_cast<int> ( (STOP_TIME-__TIME+1.e-12)/__DT );




    VectorInterface u("u"), f("f"), u0("u0"), sigma0("sigma0"), old("old"), sigma("sigma"), sigmaold("sigmaold");



    GetMultiLevelSolver()->ReInit(problemlabel);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(old);
    GetMultiLevelSolver()->ReInitVector(f);
    GetMultiLevelSolver()->ReInitVector(sigma);
    GetMultiLevelSolver()->ReInitVector(sigmaold);
    GetMultiLevelSolver()->ReInitVector(sigma0);
    GetMultiLevelSolver()->ReInitVector(u0);


    GetMultiLevelSolver()->SetProblem("vel");
    InitSolution(u);
    GetMultiLevelSolver()->SetProblem("stress");
    GetMultiLevelSolver()->GetSolver()->GetGV(sigma).CompEq(0,1.0);
    GetMultiLevelSolver()->GetSolver()->GetGV(sigma).CompEq(2,1.0);

    //    GetMultiLevelSolver()->GetSolver()->Read(u,"startu.bup");
    //    GetMultiLevelSolver()->GetSolver()->Read(sigma,"starts.bup");

    ofstream func_log("functional.txt");
    func_log.precision(12);



    GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,0);
    GetMultiLevelSolver()->GetSolver()->Visu("Results/s",sigma,0);
    GetMultiLevelSolver()->GetSolver()->Visu("Results/u0",u0,0);
    GetMultiLevelSolver()->GetSolver()->Visu("Results/s0",sigma0,0);


    double writeevery=0.1;
    double writeiter =0.0;

    for (_iter=1; _iter<=_niter; _iter++)
      {


	cout << "================================ " << _iter << ": " << __TIME << " -> "
	     << __TIME+__DT <<  endl;
	__TIME += __DT;
	writeiter+=__DT;


	GetMultiLevelSolver()->Equ(old,1.0,u);
	GetMultiLevelSolver()->Equ(u0,1.0,u);
	GetMultiLevelSolver()->Equ(sigmaold,1.0,sigma);
	GetMultiLevelSolver()->Equ(sigma0,1.0,sigma);

  /*
	// Half-Step
	__DT *= 0.5;
	///////////////////////////////////////////// 1a
	cout << endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1A" << endl;
	GetMultiLevelSolver()->SetProblem("vel"); 	BOUNDARY = true;
	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

	GetMultiLevelSolver()->AddNodeVector("old",old);
	GetMultiLevelSolver()->AddNodeVector("sigma",sigma);
	assert(Solve(u0,f)=="converged");
	GetMultiLevelSolver()->DeleteNodeVector("sigma");
	GetMultiLevelSolver()->DeleteNodeVector("old");

	///////////////////////////////////////////// 1b
	cout << endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1B" << endl;
	GetMultiLevelSolver()->SetProblem("stress");	BOUNDARY = false;
	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

	GetMultiLevelSolver()->AddNodeVector("sigmaold",sigmaold);
	GetMultiLevelSolver()->AddNodeVector("V",old);
	assert(Solve(sigma0,f)=="converged");
	GetMultiLevelSolver()->DeleteNodeVector("sigmaold");
	GetMultiLevelSolver()->DeleteNodeVector("V");
	*/

	// Full-Step
	//__DT *= 2.0;
	///////////////////////////////////////////// 2a
	cout << endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2A" << endl;
	GetMultiLevelSolver()->SetProblem("vel");	BOUNDARY = true;
	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

	GetMultiLevelSolver()->AddNodeVector("old",old);
	GetMultiLevelSolver()->AddNodeVector("sigma",sigma0);
	assert(Solve(u,f)=="converged");
	GetMultiLevelSolver()->DeleteNodeVector("sigma");
	GetMultiLevelSolver()->DeleteNodeVector("old");
	nvector<double> jkin = 	Functionals(u,f);

	///////////////////////////////////////////// 2b
	cout << endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2B" << endl;
	GetMultiLevelSolver()->SetProblem("stress");	BOUNDARY = false;
	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

	GetMultiLevelSolver()->AddNodeVector("V",u0);
	GetMultiLevelSolver()->AddNodeVector("sigmaold",sigmaold);
	assert(Solve(sigma,f)=="converged");
	GetMultiLevelSolver()->DeleteNodeVector("sigmaold");
	GetMultiLevelSolver()->DeleteNodeVector("V");


	if (writeiter+1.e-8>writeevery)
	  {
	    GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,_iter);
	    GetMultiLevelSolver()->GetSolver()->Visu("Results/s",sigma,_iter);
	    GetMultiLevelSolver()->GetSolver()->Visu("Results/u0",u0,_iter);
	    GetMultiLevelSolver()->GetSolver()->Visu("Results/s0",sigma0,_iter);

	    string name;
	    name = "Results/u";compose_name(name,_iter);
	    GetMultiLevelSolver()->GetSolver()->Write(u,name);
	    name = "Results/sigma";compose_name(name,_iter);
	    GetMultiLevelSolver()->GetSolver()->Write(sigma,name);
	    writeiter=0;
	  }

	nvector<double> jel = 	Functionals(sigma,f);

	func_log << __TIME << "\t" << jkin[0] << "\t" << jel[1] << endl;
      }

    func_log.close();
  }


  template class Loop<2>;
  template class Loop<3>;


}

