#include "loop.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "Solver/multilevelsolvers.h"
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

    
  
    FSIMultiLevelSolver<DIM>* FSI_MLS= dynamic_cast<FSIMultiLevelSolver<DIM>*>(GetMultiLevelSolver());
	
	FSI_MLS->ReInit("blubbb");
	//Anlegen von vier Loesern  
	// fsi_main 		Zum Abspeichern von Vektoren. Keine Matrix!	  	Anlegen mit Problem: "fsi_main"
	// fsi_reduced 		Matrix mit (DIM+1)*(DIM+1)						Anlegen mit Problem: "fsi_reduced"
	// meshmotion 		Matrix mit DIM*DIM								Anlegen mit Problem: "meshmotion"
	// def_solid		Matrix mit DIM*DIM								Anlegen mit Problem: "def_solid"
	
	int _initial_refine=3;
	GlobalVector vel_zw, def_pres_zw;
	string solved; 
	//////////////////////////////////////////////////////////////////////////////	
 	//////////////////////////////////////////////////////////////////////////////
    cout<<"================================================="<<endl;
    cout<<"Only Fluid Stat to compute initial condition on fluid domain(!!!)"<<endl;
        
    VectorInterface  f("f"), vel("VEL");    
    VectorInterface   def_pres("def_pres"); 
   for (int i=0;i<_initial_refine;++i)
      {
      	GetMultiLevelSolver()->ReInit("blubbb");	
        FSI_MLS->SetSolverLabel("fsi_reduced");
    	GetMultiLevelSolver()->SetProblem("fluid_stat");
    
	    GetMultiLevelSolver()->ReInitVector(vel);
        GetMultiLevelSolver()->ReInitVector(f);
		if(i==0)     GetMultiLevelSolver()->GetSolver()->GetGV(vel).zero();					
		else GetMultiLevelSolver()->InterpolateSolution(vel,vel_zw);
		
		//// SOLVE
		GetMultiLevelSolver()->GetSolver()->GetGV(f).zero();
		GetMultiLevelSolver()->GetSolver()->Rhs(f);
		GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
		// set offset first, so nodes that are both periodic and dirichlet will become dirichlet
		GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(vel);
		GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(vel);
		GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
	
		solved=GetMultiLevelSolver()->Solve(vel,f,GetSolverInfos()->GetNLInfo());
		assert(solved=="converged");
		Output(vel,"Results/init_v");
		   
	   	//////////////////////////////////////////////////////////////////////////////	
	 	//////////////////////////////////////////////////////////////////////////////
	 	cout<<"================================================="<<endl;
		cout<<"Only Solid Euler to compute prestress on solid domain!!"<<endl;
		
		FSI_MLS->SetSolverLabel("fsi_reduced");
		GetMultiLevelSolver()->SetProblem("solid_euler");
		 
   
		 
		GetMultiLevelSolver()->ReInitVector(def_pres);
		
		if(i==0)     GetMultiLevelSolver()->GetSolver()->GetGV(def_pres).zero();					
		else GetMultiLevelSolver()->InterpolateSolution(def_pres,def_pres_zw);
			
		GetMultiLevelSolver()->AddNodeVector("VEL",vel);
	
		//// SOLVE
		GetMultiLevelSolver()->GetSolver()->GetGV(f).zero();
		GetMultiLevelSolver()->GetSolver()->Rhs(f);
		GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(f);
		// set offset first, so nodes that are both periodic and dirichlet will become dirichlet
		GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(def_pres);
		GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(def_pres);
		GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
	
		solved=GetMultiLevelSolver()->Solve(def_pres,f,GetSolverInfos()->GetNLInfo());
		Output(def_pres,"Results/def_pres");
		assert(solved=="converged");
	

		GetMultiLevelSolver()->DeleteNodeVector("VEL");
	

		
		if(i<_initial_refine-1)
		{
			cout<<"============== refine mesh ================="<<endl;
			CopyVector(vel_zw,vel);
			CopyVector(def_pres_zw,def_pres);
			
			GetMeshAgent()->global_refine(1);
		}
	}	
	
	GetMultiLevelSolver()->Equ(def_pres,-1.0,def_pres);
   	//////////////////////////////////////////////////////////////////////////////	
 	//////////////////////////////////////////////////////////////////////////////	
  	cout<<"================================================="<<endl;
    cout<<"Prestressed Fluid-Structure Interaction Problem"<<endl;		
    

     SolverInfos*   _SI_Solid_Disp= new SolverInfos;
  	_SI_Solid_Disp->BasicInit(_paramfile);
    SolverInfos*   _SI_Solid_MeshMotion= new SolverInfos;
  	_SI_Solid_MeshMotion->BasicInit(_paramfile);
  	
  	_SI_Solid_Disp->GetNLInfo().control().matrixmustbebuild() = 1;
  	_SI_Solid_MeshMotion->GetNLInfo().control().matrixmustbebuild() = 1;
    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    
    FSI_MLS->SetSolverLabel("fsi_main");
    GetMultiLevelSolver()->SetProblem("fsi_main");
    
    FSI_MLS->SetSolverLabel("meshmotion");
    GetMultiLevelSolver()->SetProblem("meshmotion");
    
    FSI_MLS->SetSolverLabel("def_solid");
    GetMultiLevelSolver()->SetProblem("def_solid");
        
    FSI_MLS->SetSolverLabel("fsi_reduced");
    GetMultiLevelSolver()->SetProblem("fsi_reduced");
    GetMultiLevelSolver()->AddNodeVector("DEF_PRES",def_pres);
    
    
    FSI_MLS->SetSolverLabel("fsi_main")	;
    VectorInterface U_Vec("U_Vec"),UOLD_Vec("UOLD_Vec"), F("F");
    
    FSI_MLS->ReInitVector(U_Vec);
    FSI_MLS->GetSolver()->GetGV(U_Vec).zero();
    FSI_MLS->ReInitVector(UOLD_Vec);
    FSI_MLS->GetSolver()->GetGV(UOLD_Vec).zero();
    FSI_MLS->ReInitVector(F);
    FSI_MLS->GetSolver()->GetGV(F).zero();

    
    FSI_MLS->AddNodeVectorinAllSolvers("U_Vec",U_Vec);
    FSI_MLS->AddNodeVectorinAllSolvers("UOLD_Vec",UOLD_Vec);

	//Init Condition aus Vel Vektor ubertragen
	GlobalVector &U_GV =GetMultiLevelSolver()->GetSolver()->GetGV(U_Vec); 
	for (int node=0;node<GetMultiLevelSolver()->GetSolver()->GetMesh()->nnodes();++node)
	  {
		for(int i=0;i<DIM+1;i++)
	  		U_GV(node,i) = FSI_MLS->GetSolver("fsi_reduced")->GetGV(vel)(node,i); 	
	  	for(int i=0;i<DIM;i++)
	  		U_GV(node,DIM+1+i)=0.0	;	
	  }
	  
	Output(U_Vec,"Results/init");
    //InitSolution(u);

    for (_iter=1; _iter<=_niter; _iter++)
      {
		cout << "========== " << _iter << ": " << __TIME << " -> " 
				 << __TIME+__DT << "  (" << __THETA << ")" << endl;
		__TIME += __DT;

		GetMultiLevelSolver()->Equ(UOLD_Vec,1.0,U_Vec);

		//// SOLVE
		GetMultiLevelSolver()->GetSolver()->Zero(F);
		GetMultiLevelSolver()->GetSolver()->Rhs(F);
		GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(F);
		// set offset first, so nodes that are both periodic and dirichlet will become dirichlet
		GetMultiLevelSolver()->GetSolver()->SetPeriodicVector(U_Vec);
		GetMultiLevelSolver()->GetSolver()->SetBoundaryVector(U_Vec);
		
		solved=FSI_MLS->Solve(U_Vec,F,GetSolverInfos()->GetNLInfo(),_SI_Solid_Disp->GetNLInfo(),_SI_Solid_MeshMotion->GetNLInfo());
		assert(solved=="converged");
		Output(U_Vec,"Results/uuuuuu");
	  
		if (_iter%5==0)
		{
		GetMultiLevelSolver()->GetSolver()->Visu("Results/U",U_Vec,_iter);
		WriteMeshAndSolution("Results/U",U_Vec);
		}

		/*		
		DoubleVector juh = Functionals(u,f);
		DoubleVector juh2 = Functionals(def,f);
	
		InterfaceResidual=true;
		DoubleVector juh3 = Functionals(u,f);
		InterfaceResidual=false;
		*/

		//func_log <<__TIME << "\t" << 0.0 << "\t" << juh<<"\t"<<juh3[5]<<"  "<<juh3[7]<<"  "<<juh2[8]<<"  "<<juh2[9]<< endl;
    }
    FSI_MLS->SetSolverLabel("fsi_main"); 
    	FSI_MLS->DeleteNodeVectorinAllSolvers("U_Vec");
    	FSI_MLS->DeleteNodeVectorinAllSolvers("UOLD_Vec");
    
    FSI_MLS->SetSolverLabel("fsi_reduced")	;
		GetMultiLevelSolver()->DeleteNodeVector("DEF_PRES");

    func_log.close();

  }
  

  template class Loop<2>;
  template class Loop<3>;


}





