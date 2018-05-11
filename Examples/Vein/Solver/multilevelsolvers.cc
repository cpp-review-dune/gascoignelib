#include "multilevelsolvers.h"
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include <algorithm>
#include  "solvers.h"


using namespace std;

extern double __DT,__THETA;



namespace Gascoigne
{
	template<int DIM>
    FSIMultiLevelSolver<DIM>::FSIMultiLevelSolver(): StdMultiLevelSolver() , _f("_f"), _u("_u")
    {}
    

	template<int DIM>
    FSIMultiLevelSolver<DIM>::~FSIMultiLevelSolver()
	{
	  //ViewProtocoll();
		
	  for(auto itt : _MSP)
	  {
	  	  SetSolverLabel(itt.first);
		  if(DataP) delete DataP;
		  DataP=NULL;

		  for(int i=0;i<GetSolverPointers().size();i++) 
			{ 
			  if (GetSolverPointers()[i]) 
				{
				  delete GetSolverPointers()[i]; 
				  GetSolverPointers()[i]=NULL; 
				}
			}
	  }
	   for(int i=0; i<_Interpolator.size(); i++)  
		{
		  if (_Interpolator[i]) 
		    {
		      delete _Interpolator[i]; 
		      _Interpolator[i]=NULL;
		    }
		}
	}
	/*-----------------------------------------------------------------------*/
    template<int DIM>
    void FSIMultiLevelSolver<DIM>::ReInit(const std::string& problemlabel)
	{
	  DataP->CountResidual() = 0;
			// DataP->GetNLInfo().control().matrixmustbebuild() = 1;
	  //MgInterpolator		
	  NewMgInterpolator();

	  //Solver1
	  SetSolverLabel("fsi_main");
	  NewSolvers("fsi_main");
	  SolverNewMesh();
	  MgInterpolatorToSolver();
	  SetProblem("fsi_main");
	  //RegisterMatrix();
	  //ReInitMatrix();	  	  
	  
	  //In jedem Solver gibt es jetzt diese Vektoren bis au fsi_main!!

	  SetSolverLabel("fsi_reduced");
	  NewSolvers("fsi_reduced");
	  SolverNewMesh();
	  MgInterpolatorToSolver();
	  SetProblem("fsi_reduced");
	  RegisterMatrix();
	  ReInitMatrix();	  

        
	  //In jedem Solver gibt es jetzt diese Vektoren!!
	  RegisterVectors();	
	  ReInitVector(_cor);
	  ReInitVector(_res);
	  ReInitVector(_mg0);
	  ReInitVector(_mg1);
	  ReInitVector(_f);
	  ReInitVector(_u);
	  	  
	  SetSolverLabel("meshmotion");
	  NewSolvers("meshmotion");
	  SolverNewMesh();
	  MgInterpolatorToSolver();
	  SetProblem("meshmotion");
	  RegisterMatrix();
	  ReInitMatrix();	  	  

	  //In jedem Solver gibt es jetzt diese Vektoren!!
	  RegisterVectors();	
	  ReInitVector(_cor);
	  ReInitVector(_res);
	  ReInitVector(_mg0);
	  ReInitVector(_mg1);
	  ReInitVector(_f);
	  ReInitVector(_u);
	  	  
	  SetSolverLabel("def_solid");
	  NewSolvers("def_solid");
	  SolverNewMesh();
	  MgInterpolatorToSolver();
	  SetProblem("def_solid");
	  RegisterMatrix();
	  ReInitMatrix();	  	  

	  //In jedem Solver gibt es jetzt diese Vektoren!!
	  RegisterVectors();	
	  ReInitVector(_cor);
	  ReInitVector(_res);
	  ReInitVector(_mg0);
	  ReInitVector(_mg1);
	  ReInitVector(_f);
	  ReInitVector(_u);
	  
	  SetSolverLabel("fsi_reduced");
	}
	/*-----------------------------------------------------------------------*/ 
	template<int DIM>
    void FSIMultiLevelSolver<DIM>::RegisterVectors() 
	{
  		assert(nlevels()==GetSolverPointers().size());
	  for (int level=0; level<nlevels(); ++level)  
		{
		  GetSolver(level)->RegisterVector(_cor);
		  GetSolver(level)->RegisterVector(_res);
		  GetSolver(level)->RegisterVector(_mg0);
		  GetSolver(level)->RegisterVector(_mg1);
		  GetSolver(level)->RegisterVector(_f);
		  GetSolver(level)->RegisterVector(_u);
		}
	}
	/*-----------------------------------------------------------------------*/ 
	
	   
   template<int DIM>
    void FSIMultiLevelSolver<DIM>::NewMgInterpolator()
	{
	  for (int i=0;i<_Interpolator.size();++i)
		{
		  assert(_Interpolator[i]!=NULL);
		  delete _Interpolator[i];
		  _Interpolator[i]=NULL;
		}
	  _Interpolator.resize(nlevels()-1,NULL);

	  for(int l=0; l<nlevels()-1; ++l)  
		{
		  //_Interpolator[l] = new MgInterpolatorMatrix;
		  _Interpolator[l] = new MgInterpolatorNested;
		}
	}
	/*-----------------------------------------------------------------------*/
        template<int DIM>
    void FSIMultiLevelSolver<DIM>::MgInterpolatorToSolver()
	{
	  for (int level=0;level<nlevels()-1;++level)
		{
		  int sl = nlevels()-level-2;

		  const MeshTransferInterface* MT = GetMeshAgent()->GetTransfer(sl);
		  assert(MT);
		  assert(_Interpolator[level]);
		  GetSolver(level)->ConstructInterpolator(_Interpolator[level],MT);
		}
	}
	/*-----------------------------------------------------------------------*/
	
    template<int DIM>
    void FSIMultiLevelSolver<DIM>::NewSolvers(std::string solverlabel)
	{
	  assert(_SolverLabel==solverlabel);
	  if(_MSP.find(solverlabel)==_MSP.end())
	  {
	  	std::vector<SolverInterface*>	zw_solver_int;
	  	_MSP.insert(pair<std::string,std::vector<SolverInterface*>> (solverlabel,zw_solver_int));
	  }
	  
	  oldnlevels = GetSolverPointers().size();

	  if (oldnlevels>nlevels())
		{
		  for (int l=oldnlevels-1; l>=nlevels(); l--)
		    {
		      delete GetSolverPointers()[l];
		      GetSolverPointers()[l] = NULL;
		    }
		}
	  GetSolverPointers().resize(nlevels(),NULL);
	  ComputeLevel = GetSolverPointers().size()-1;

	  for(int level=0; level<nlevels(); ++level)  
		{
		  int solverlevel = nlevels()-1-level;

		  // new Solvers
		  if(GetSolver(solverlevel)==NULL) 
		    {
		      GetSolverPointer(solverlevel) = NewSolver(solverlevel);
		      cout<<"ATTENTION Solver"<<solverlabel<<" gets paramfile from problem"<< solverlabel <<"!!!!!!"<<endl;
		      GetSolver(solverlevel)->BasicInit(GetProblemContainer()->GetProblem(solverlabel)->GetParamFile(),GetMeshAgent()->GetDimension());
		      dynamic_cast<FSISolver<DIM>*>(GetSolver(solverlevel))->SetSolverLabel(solverlabel);
		    }
		}
			  
	} 
  	/*-----------------------------------------------------------------------*/

  
    template<int DIM>
  void FSIMultiLevelSolver<DIM>::NewtonUpdateDisplacement(GlobalVector& U_Vec_GV,NLInfo& info_solid_disp,NLInfo& info_meshmotion)
  {		
		SetSolverLabel("def_solid");
		NewtonVectorZero(_cor);
		NewtonVectorZero(_f);
	  	double rrr = NewtonResidual(_res,_cor,_f);
	  	

	 	NewtonMatrixControl(_cor,info_solid_disp);
		NewtonVectorZero(_cor);
		NewtonLinearSolve(_cor,_res,info_solid_disp.GetLinearInfo());
		
		for (int node=0;node<GetSolver(ComputeLevel)->GetMesh()->nnodes();++node)
		  {
			for(int i=0;i<DIM;i++)
		  		U_Vec_GV(node,1+DIM+i) += GetSolver(ComputeLevel)->GetGV(_cor)(node,i); 		
		  }
		rrr = NewtonResidual(_res,_cor,_f);	  

		SetSolverLabel("meshmotion");
		NewtonVectorZero(_cor);
		NewtonVectorZero(_f);
	  	double rrrr = NewtonResidual(_res,_cor,_f);

	 	NewtonMatrixControl(_cor,info_meshmotion);
		NewtonVectorZero(_cor);
		NewtonLinearSolve(_cor,_res,info_meshmotion.GetLinearInfo());
	 	
		for (int node=0;node<GetSolver(ComputeLevel)->GetMesh()->nnodes();++node)
		  {
			for(int i=0;i<DIM;i++)
		  		U_Vec_GV(node,1+DIM+i)+=GetSolver(ComputeLevel)->GetGV(_cor)(node,i); 		
		  }
		rrrr = NewtonResidual(_res,_cor,_f);

	  	
		SetSolverLabel("fsi_reduced");
	
  }
  /*-----------------------------------------------------------------------*/
  
  template<int DIM>
  void FSIMultiLevelSolver<DIM>::AddNodeVectorinAllSolvers(const string& name, VectorInterface& gq)
{
  Transfer(ComputeLevel,1,gq);
  for(int l=0; l<nlevels(); l++)
    {
    	//Geht davon aus, dass gq in aktuellem Solver gespeichert
      for(auto it: _MSP)	
      	dynamic_cast<FSISolver<DIM>*>(it.second[l])->AddNodeVector(name,GetSolver(l)->GetGV(gq));
      	
    }
}
  /*-----------------------------------------------------------------------*/
  
  template<int DIM>
  void FSIMultiLevelSolver<DIM>::DeleteNodeVectorinAllSolvers(const string& name)
{
  for(int l=0; l<nlevels(); l++)
    {
    	//Geht davon aus, dass gq in aktuellem Solver gespeichert
      for(auto it: _MSP)	
      	dynamic_cast<FSISolver<DIM>*>(it.second[l])->DeleteNodeVector(name);      	
    }

}
/*-----------------------------------------------------------------------*/
  template<int DIM>
  void FSIMultiLevelSolver<DIM>::UpdateVelocityandPressureinU(GlobalVector& U_Vec_GV, double factor, VectorInterface& dx)
	{

		SetSolverLabel("fsi_reduced")	;
		for (int node=0;node<GetSolver(ComputeLevel)->GetMesh()->nnodes();++node)
		  {
			for(int i=0;i<DIM+1;i++)
		  		U_Vec_GV(node,i) +=factor*GetSolver(ComputeLevel)->GetGV(dx)(node,i); 		
		  }  

	}
	
 /*-----------------------------------------------------------------------*/
  template<int DIM>
  void FSIMultiLevelSolver<DIM>::VelocityandPressureUto_u(GlobalVector& U_Vec_GV, VectorInterface& _u)
	{
			SetSolverLabel("fsi_reduced")	;
			for (int node=0;node<GetSolver(ComputeLevel)->GetMesh()->nnodes();++node)
			  {
				for(int i=0;i<DIM+1;i++)
			  		GetSolver(ComputeLevel)->GetGV(_u)(node,i) = U_Vec_GV(node,i); 		
			  }  
		
	}	
  /*-----------------------------------------------------------------------*/

  template<int DIM>
  double FSIMultiLevelSolver<DIM>::NewtonUpdate(GlobalVector& U_Vec_GV, double& rr, VectorInterface& x, VectorInterface& dx, VectorInterface& r, const VectorInterface& f, NLInfo& nlinfo,NLInfo& info_solid_disp,NLInfo& info_meshmotion)
  {

    const CGInfo& linfo = nlinfo.GetLinearInfo();
    bool lex  = linfo.control().status()=="exploded";
    
    double nn = NewtonNorm(dx);

    double nr = GetSolver(ComputeLevel)->Norm(r);
    
    if (nn>1.e30)  lex =1;
    if (!(nn>=0.)) lex =1;
    if (nr>1.e30)  lex =1;
    if (!(nr>=0.)) lex =1;
    
    if(lex)
      {
		nlinfo.control().status()="diverged";
		cerr << "linear : " << linfo.control().status() << endl;
		cerr << "nonlinear : " << nn << endl;
		return NewtonNorm(dx);
      }
    
    double omega = 0.7;
    double relax = 1.;
    
    GetSolver(ComputeLevel)->SetPeriodicVectorZero(dx);
    
    //VectorInterface XXX("old");
    //GlobalVector OLD = GetSolver()->GetGV(XXX);
    //For LPS Stabilization
    GetSolver(ComputeLevel)->Add(x,relax,dx);
	UpdateVelocityandPressureinU(U_Vec_GV,relax, dx);
	
    NewtonUpdateDisplacement(U_Vec_GV,info_solid_disp,info_meshmotion);
    
    NewtonResidual(r,x,f);
    rr = NewtonNorm(r);
    
 
    string message = "";
    for(int iter=0;iter<nlinfo.user().maxrelax();iter++)
      { 
		message = nlinfo.check_damping(iter,rr);
	
		if (message=="ok")       break;
		if (message=="continue") 
		  {
			//For LPS Stabilization
			GetSolver(ComputeLevel)->Add(x,relax*(omega-1.),dx);
			UpdateVelocityandPressureinU(U_Vec_GV,relax*(omega-1.), dx);
		
			NewtonUpdateDisplacement(U_Vec_GV,info_solid_disp,info_meshmotion);	    
			
			NewtonResidual(r,x,f);
			rr = NewtonNorm(r);

			relax *= omega;
			continue;
		  }
		if (message=="exploded")
		  {
			//For LPS Stabilization
			GetSolver(ComputeLevel)->Add(x,-relax,dx);
			UpdateVelocityandPressureinU(U_Vec_GV,-relax, dx);
		
			NewtonUpdateDisplacement(U_Vec_GV,info_solid_disp,info_meshmotion);
			relax = 0.;
			cout << "Damping exploded !!!!!" << endl;
			nlinfo.control().status() = "diverged";
			break;
		  }
      }
    return NewtonNorm(dx);
  }
  
 
  /*-----------------------------------------------------------------------*/
  
  template<int DIM>
  void FSIMultiLevelSolver<DIM>::newton(VectorInterface& U_Vec, const VectorInterface& f, VectorInterface& r, VectorInterface& w, NLInfo& info_fsi_reduced ,NLInfo& info_solid_disp,NLInfo&  info_meshmotion)
	{
	  info_fsi_reduced.reset();
	  info_solid_disp.reset();
	  info_meshmotion.reset();

	  SetSolverLabel("fsi_main");
	  GlobalVector &U_Vec_GV =GetSolver(ComputeLevel)->GetGV(U_Vec); 
	  
	  NewtonUpdateDisplacement(U_Vec_GV,info_solid_disp,info_meshmotion);	
  	  info_solid_disp.control().iteration()=0;
	  info_meshmotion.control().iteration()=0;
	  
	  cout<<"Res solid_disp: "; NewtonOutput(info_solid_disp);
	  cout<<"Res meshmotion: ";	NewtonOutput(info_meshmotion);
	  
	  SetSolverLabel("fsi_reduced");
	  //For LPS Stabilization
	  VelocityandPressureUto_u(U_Vec_GV,_u);
	  double rr = NewtonResidual(_res,_u,_f);
	  //GetSolver()->Visu("newton",_res,GetSolver()->GetGV(_res).n());
	  bool reached = info_fsi_reduced.check(0,rr,0.);
	  cout<<"fsi reduced: "; NewtonOutput(info_fsi_reduced);

	  for(int it=1; !reached; it++)
		{
		  NewtonMatrixControl(_u,info_fsi_reduced);
		  NewtonVectorZero(_cor);
		  NewtonLinearSolve(_cor,_res,info_fsi_reduced.GetLinearInfo());
		  double rw = NewtonUpdate(U_Vec_GV, rr,_u,_cor,_res,_f,info_fsi_reduced,info_solid_disp,info_meshmotion);
		  reached = info_fsi_reduced.check(it,rr,rw);
		  info_solid_disp.control().iteration()=it;
		  info_meshmotion.control().iteration()=it;

		  cout<<"fsi reduced: "; NewtonOutput(info_fsi_reduced);
		  cout<<"solid_disp:  "; NewtonOutput(info_solid_disp);
	  	  cout<<"meshmotion:  "; NewtonOutput(info_meshmotion);
		}
	SetSolverLabel("fsi_main");	

	}
	
	/*-----------------------------------------------------------------------*/
	  template<int DIM>
  string FSIMultiLevelSolver<DIM>::Solve( VectorInterface& U_Vec, const VectorInterface& b, NLInfo& info_fsi_reduced ,NLInfo& info_solid_disp,NLInfo& info_meshmotion)
	{
	  ComputeLevel = nlevels()-1;

	  string status;

	  assert(DataP->NonLinearSolve() == "newton");

	  GetSolver(ComputeLevel)->HNAverage(U_Vec);
	  newton(U_Vec,b,_res,_cor,info_fsi_reduced, info_solid_disp, info_meshmotion);
	  GetSolver(ComputeLevel)->HNZero(U_Vec);
	  return info_fsi_reduced.CheckMatrix();

	}
	
	
		/*-----------------------------------------------------------------------*/
	  template<int DIM>
  void FSIMultiLevelSolver<DIM>::AssembleMatrix(VectorInterface& u)
	{
	
	  //SolutionTransfer(u);
	    Transfer(ComputeLevel,1,u);

 		//for(int l=ComputeLevel;l>=1;l--)
    	  //GetSolver(l-1)->SetBoundaryVector(u);

	  for(int l=0;l<=ComputeLevel;l++)
		{
		  GetSolver(l)->MatrixZero();
		  GetSolver(l)->AssembleMatrix(u,1.);
		}
	}
	/*-----------------------------------------------------------------------*/
	 template<int DIM>
	 void FSIMultiLevelSolver<DIM>::ComputeIlu(VectorInterface& u)
		{
			  //SolutionTransfer(u);
				Transfer(ComputeLevel,1,u);

		 		//for(int l=ComputeLevel;l>=1;l--)
				//GetSolver(l-1)->SetBoundaryVector(u);

		  for(int l=0;l<=ComputeLevel;l++)
			{
			  GetSolver(l)->ComputeIlu(u);
			}
		}
  template class FSIMultiLevelSolver<2>;
  template class FSIMultiLevelSolver<3>;
}
