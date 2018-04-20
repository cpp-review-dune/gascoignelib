#include "loop.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"

using namespace std;


double __DT,__TIME, __THETA;
bool InterfaceResidual=false;




namespace Gascoigne
{
  template<int DIM>
  void Loop<DIM>::run(const std::string& problemlabel)
  {
    double STOP_TIME;

	DataFormatHandler DFH;
	DFH.insert("start_time" ,    &__TIME , 0.0);
	DFH.insert("stop_time" ,    &STOP_TIME , 0.0);
	DFH.insert("theta" ,    &__THETA , 0.0);
	DFH.insert("dt" ,    &__DT,0.0);
	FileScanner FS(DFH, _paramfile, "Equation");
	assert(STOP_TIME>__TIME);
	assert(__DT>0);
	assert(__THETA>0);
      
     
    _niter = static_cast<int> ( (STOP_TIME-__TIME+1.e-12)/__DT );
	ofstream func_log("functional.txt");
    func_log.precision(12);
  

    // damp at 1, 2, a.s.o.
    int DI = static_cast<int> (1.0 / __DT);
    cout << DI << endl;
    StopWatch _all;

  
/*
     for (int i=0;i<2;++i)
      {
				GlobalVector ualt;
				CopyVector(ualt,u);   
				IntVector empty(0);  

				const IntVector* refnodes1 = GetMultiLevelSolver()->GetSolver()->GetMesh()->VertexOnBoundary(80);
				const IntVector* refnodes2 = GetMultiLevelSolver()->GetSolver()->GetMesh()->VertexOnBoundary(81);
				const IntVector* refnodes3 = GetMultiLevelSolver()->GetSolver()->GetMesh()->VertexOnBoundary(90);

				IntVector refnodescorner;
				refnodescorner.insert(refnodescorner.begin(),refnodes1->begin(),refnodes1->end());
				refnodescorner.insert(refnodescorner.begin(),refnodes2->begin(),refnodes2->end());
				refnodescorner.insert(refnodescorner.begin(),refnodes3->begin(),refnodes3->end());
		
				GetMeshAgent()->refine_nodes(refnodescorner,empty);	
	
				GetMultiLevelSolver()->ReInit(problemlabel);
				GetMultiLevelSolver()->ReInitVector(u);
				GetMultiLevelSolver()->ReInitVector(old);
				GetMultiLevelSolver()->ReInitVector(f);
				GetMultiLevelSolver()->ReInitVector(def);
				GetMultiLevelSolver()->ReInitVector(defold);
				GetMultiLevelSolver()->InterpolateSolution(u,ualt);
		}
		////////////////////////////////////////////////////////////////////////////////////
  */
  
    //////////////////////////////////////////////////////////////////////////////
    //Initialiseren der Vektoren
    VectorInterface u("u"), f("f"), old("old"), vel("vel"), def("def"), defold("defold"),def_pres("def_pres");
 

	
    GetMultiLevelSolver()->ReInit("fsi");
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(old);
    GetMultiLevelSolver()->ReInitVector(vel);
    GetMultiLevelSolver()->ReInitVector(f);
    GetMultiLevelSolver()->ReInitVector(def);
    GetMultiLevelSolver()->ReInitVector(defold);
    GetMultiLevelSolver()->ReInitVector(def_pres);
	//////////////////////////////////////////////////////////////////////////////	
 	//////////////////////////////////////////////////////////////////////////////
    cout<<"================================================="<<endl;
    cout<<"Only Fluid Stat to compute initial condition"<<endl;

    GetMultiLevelSolver()->SetProblem("fluid_stat");
	//// SOLVE
	GetMultiLevelSolver()->GetSolver()->GetGV(u).zero();
	GetMultiLevelSolver()->GetSolver()->GetGV(f).zero();
	GetMultiLevelSolver()->GetSolver()->Rhs(f);
	GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
	// set offset first, so nodes that are both periodic and dirichlet will become dirichlet
	GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(u);
	GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);

    
	string solved=GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());
	assert(solved=="converged");
		
	Output(u,"Results/init_v");
	
    //Abspeichern des Geschwindigkeit als Anfangsbed
    GetMultiLevelSolver()->Equ(vel,1.0,u);
   
   	//////////////////////////////////////////////////////////////////////////////	
 	//////////////////////////////////////////////////////////////////////////////
 	cout<<"================================================="<<endl;
    cout<<"Only Solid Euler to compute prestress"<<endl;
    

    GetMultiLevelSolver()->SetProblem("solid_euler");
	GetMultiLevelSolver()->GetSolver()->GetGV(u).zero();
	GetMultiLevelSolver()->AddNodeVector("VEL",vel);
	//// SOLVE

	GetMultiLevelSolver()->GetSolver()->GetGV(f).zero();
	GetMultiLevelSolver()->GetSolver()->Rhs(f);
	GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
	// set offset first, so nodes that are both periodic and dirichlet will become dirichlet
	GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(u);
	GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
	solved=GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());
	Output(u,"Results/def_pres");
	assert(solved=="converged");
	
    GetMultiLevelSolver()->Equ(def_pres,-1.0,u);
	GetMultiLevelSolver()->DeleteNodeVector("VEL");
	
   	//////////////////////////////////////////////////////////////////////////////	
 	//////////////////////////////////////////////////////////////////////////////	
  	cout<<"================================================="<<endl;
    cout<<"Prestressed Fluid-Structure Interaction Problem"<<endl;		
    
    GetMultiLevelSolver()->SetProblem("fsi");
    GetMultiLevelSolver()->GetSolver()->Zero(def);
    
    GetMultiLevelSolver()->Equ(u,1.0,vel);
	Output(u,"Results/init");
    //InitSolution(u);
    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    for (_iter=1; _iter<=_niter; _iter++)
      {
	cout << "========== " << _iter << ": " << __TIME << " -> " 
			 << __TIME+__DT << "  (" << __THETA << ")" << endl;
	__TIME += __DT;

	GetMultiLevelSolver()->Equ(old,1.0,u);
	GetMultiLevelSolver()->Equ(defold,1.0,def);

	GetMultiLevelSolver()->AddNodeVector("OLD",old);
	GetMultiLevelSolver()->AddNodeVector("DEFOLD",defold);
	GetMultiLevelSolver()->AddNodeVector("DEF",def);
    GetMultiLevelSolver()->AddNodeVector("DEF_PRES",def_pres);

	//	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
	//assert(Solve(u,f)=="converged");
	//// SOLVE
	GetMultiLevelSolver()->GetSolver()->Zero(f);
	GetMultiLevelSolver()->GetSolver()->Rhs(f);
	GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
	// set offset first, so nodes that are both periodic and dirichlet will become dirichlet
	GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(u);
	GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(u);
	dynamic_cast<FSIMultiLevelSolver<DIM>*>(GetMultiLevelSolver())->UpdateDeformation(u);

	solved=GetMultiLevelSolver()->Solve(u,f,GetSolverInfos()->GetNLInfo());
	assert(solved=="converged");
	Output(u,"Results/uuuuuu");
  
	if (_iter%5==0)
	{
	GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,_iter);
	WriteMeshAndSolution("Results/u",u);
	}

	/*		
	DoubleVector juh = Functionals(u,f);
	DoubleVector juh2 = Functionals(def,f);
	
	InterfaceResidual=true;
	DoubleVector juh3 = Functionals(u,f);
	InterfaceResidual=false;
	*/
	GetMultiLevelSolver()->DeleteNodeVector("OLD");
	GetMultiLevelSolver()->DeleteNodeVector("DEFOLD");
	GetMultiLevelSolver()->DeleteNodeVector("DEF");
	GetMultiLevelSolver()->DeleteNodeVector("DEF_PRES");
	//func_log <<__TIME << "\t" << 0.0 << "\t" << juh<<"\t"<<juh3[5]<<"  "<<juh3[7]<<"  "<<juh2[8]<<"  "<<juh2[9]<< endl;

    }

    func_log.close();

  }
  

  template class Loop<2>;
  template class Loop<3>;


}





