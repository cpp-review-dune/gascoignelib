#include "loop.h"

using namespace std;

double TIME,DT,THETA;

namespace Gascoigne
{

  
  void Loop::run(const std::string& problemlabel)
  {
    VectorInterface u("u"), f("f"), old("old");
    

    TIME = 0.0;

    DataFormatHandler DFH;
    DFH.insert("dt" ,    &DT , 0.0);
    DFH.insert("theta" , &THETA , 0.0);
    FileScanner FS(DFH, _paramfile, "Equation");
    assert(DT>0);
    assert(THETA>0);


    GetMultiLevelSolver()->ReInit(problemlabel);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(old);
    GetMultiLevelSolver()->ReInitVector(f);
    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

    InitSolution(u);
	
    
    for (_iter=1; _iter<=_niter; _iter++)
      {
	
	TIME += DT;
	cout << "\n Time step " << TIME-DT << " -> " << TIME << "\t dt = " << DT << endl;

	GetMultiLevelSolver()->Equ(old,1.0,u);
	
	
	GetMultiLevelSolver()->AddNodeVector("OLD",old);
	Solve(u,f);
	GetMultiLevelSolver()->DeleteNodeVector("OLD");
	
	DoubleVector juh = Functionals(u,f);
      }
  }

  
  
}





