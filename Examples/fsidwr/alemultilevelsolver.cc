#include "alemultilevelsolver.h"
#include "aleinterpolator.h"


using namespace std;

namespace Gascoigne
{
  extern StopWatch S4;
  extern StopWatch S_LIN;
  // //////////////////////////////////////////////////

  void AleMultiLevelSolver::smooth(int l, int niter, VectorInterface& x, VectorInterface& b, VectorInterface& h)
  {
    S4.start();
    
    assert(dynamic_cast<const StdSolver*> (GetSolver(l)));
    bool direct = dynamic_cast<const StdSolver*> (GetSolver(l))->DirectSolver();
    bool splitting_smoother  = GetAleSolver(l)->SplittingSmoother();
    bool splitting_fluid_exact  = GetAleSolver(l)->SplittingFluidExact();
    bool splitting_solid_exact  = GetAleSolver(l)->SplittingSolidExact();
    
    
    if (direct)  
      GetSolver(l)->smooth_exact(x,b,h);       // Exact Coupled Solving
    else         
      {
	if (splitting_smoother)                // Splitting Smoother
	  {
	    for (int iter=0;iter<niter;++iter)
	      {

		if (0)
		  {
		    GetSolver(l)->MatrixResidual(h,x,b);
		    std::cout << l << "\t" << iter << "\t 3: "
			      << GetSolver(l)->GetGV(h).norm() << std::endl;
		  }
		


		GetSolver(l)->MatrixResidual(h,x,b);	
		if(0)
		  std::cout << l << "\t" << iter << "\t 1: "
			    << GetSolver(l)->GetGV(h).norm() << std::endl;
		






		if (splitting_solid_exact) GetAleSolver(l)->SolveSolidExact(x,h);
		else                       GetAleSolver(l)->SolveSolidIterative(20,0.0001, x,h);

		GetSolver(l)->MatrixResidual(h,x,b);

		if (splitting_fluid_exact) GetAleSolver(l)->SolveFluidExact(x,h);
		else                       GetAleSolver(l)->SolveFluidIterative(2000,0.0001, x,h);
	      }
	  }
	else                                   // Coupled Std Smoother
	  {
	    GetSolver(l)->smooth(niter,x,b,h);  
	  }
      }
    
    S4.stop();
    
  }
  




  
  /*-------------------------------------------------------------*/
 
  void AleMultiLevelSolver::mgstep(vector<double>& res, vector<double>& rw, 
				   int l, int finelevel, int coarselevel, string& p0, string p,
				   VectorInterface& u, VectorInterface& b, VectorInterface& v)
  {
    assert(dynamic_cast<const StdSolver*> (GetSolver(l)));
    int pre   = dynamic_cast<const StdSolver*> (GetSolver(l))->GetSolverData().GetIterPre();
    int post  = dynamic_cast<const StdSolver*> (GetSolver(l))->GetSolverData().GetIterPost();
    int exact = dynamic_cast<const StdSolver*> (GetSolver(l))->GetSolverData().GetIterExact();
    
    if(l==coarselevel)
      {
	if(p=="F") {p0="V";}
	smooth(l,exact,u,b,v);
	
	if(l==finelevel)
	  {
	    GetSolver(l)->MatrixResidual(v, u, b);
	    res[l] = GetSolver(l)->Norm(v);
	  }
      }
    else
      {
	smooth(l,pre,u,b,v);

	GetSolver(l)->MatrixResidual(v,u,b);
      
	_Interpolator[l-1]-> restrict_zero(GetSolver(l-1)->GetGV(b),GetSolver(l)->GetGV(v));
	GetSolver(l-1)->HNDistribute(b);
	GetSolver(l-1)->SetBoundaryVectorZero(b);
	GetSolver(l-1)->Zero(u);
      
	int j = 0;
	if (p0=="V") j = 1;
	if (p0=="W") j = 2;
	if (p0=="F") j = 3;
	for (int i = 0; i<j; i++)
	  {
	    mgstep(res,rw,l-1,finelevel,coarselevel,p0,p,u,b,v);
	  }
	if ((l==0)&&(p=="F")) { p0="W";}
	rw[l] = GetSolver(l-1)->Norm(u);

	GetSolver(l)->Zero(v);
	GetSolver(l-1)    -> HNAverage(u);
	_Interpolator[l-1]-> prolongate_add(GetSolver(l)->GetGV(v),GetSolver(l-1)->GetGV(u));
	GetSolver(l-1) -> HNZero(u);
	GetSolver(l)   -> HNZero(v);
     
	GetSolver(l)   -> SetBoundaryVectorZero(v);

	GetSolver(l)->Add(u,DataP->MgOmega(),v);
    
	smooth(l,post,u,b,v);
	
	GetSolver(l)->MatrixResidual(v,u,b);
	res[l] = GetSolver(l)->Norm(v);
      }
  }

  /*-------------------------------------------------------------*/
  
  void AleMultiLevelSolver::NewtonLinearSolve(VectorInterface& x, const VectorInterface& b, CGInfo& info)
  {
    S_LIN.start();
    StdMultiLevelSolver::NewtonLinearSolve(x,b,info);
    S_LIN.stop();
  }
  
  void AleMultiLevelSolver::NewMgInterpolator()
  {
    for (int i=0;i<_Interpolator.size();++i)
      {
	assert(_Interpolator[i]!=NULL);
	delete _Interpolator[i];
	_Interpolator[i]=NULL;
      }
    _Interpolator.resize(nlevels()-1,NULL);
    
    for(int l=0; l<nlevels()-1; ++l)  
      _Interpolator[l] = new AleInterpolator;
    
    for (int level=0;level<nlevels()-1;++level)
      {
	int sl = nlevels()-level-2;
	
	const MeshTransferInterface* MT = GetMeshAgent()->GetTransfer(sl);
	assert(MT);
	assert(_Interpolator[level]);
	GetSolver(level)->ConstructInterpolator(_Interpolator[level],MT);
      }
  }
}
