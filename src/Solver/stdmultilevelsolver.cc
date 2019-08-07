/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2011 by the Gascoigne 3D authors
 *
 * This file is part of Gascoigne 3D
 *
 * Gascoigne 3D is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version.
 *
 * Gascoigne 3D is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * Please refer to the file LICENSE.TXT for further information
 * on this license.
 *
 **/


#include  "stdmultilevelsolver.h"
#include  "stdtimesolver.h"
#include  "compose_name.h"
#include  "gascoignemultigridmesh.h"
#include  "cg.h"
#include  <iomanip>
#include  "mginterpolatormatrix.h"
#include  "mginterpolatornested.h"
#include  "gmres.h"

using namespace std;

/*-------------------------------------------------------------*/

namespace Gascoigne
{

  extern Timer GlobalTimer;
  
  StdMultiLevelSolver::~StdMultiLevelSolver()
  {
    //ViewProtocoll();

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
    for(int i=0; i<_Interpolator.size(); i++)  
      {
	if (_Interpolator[i]) 
	  {
	    delete _Interpolator[i]; 
	    _Interpolator[i]=NULL;
	  }
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::ViewProtocoll() const
  {
    // old...
    GlobalTimer.print100();
    assert(0);
    // cout << "\n************************************************************************\n\n";
    // cout << "StdMultiLevelSolver\tTIME\n";
    // cout << "  Residual\t\t" << _clock_residual.read() << endl;
    // cout << "  Solve\t\t\t" << _clock_solve.read() << endl;
    // cout << "MATRIX\t\t\tTIME\n";

    // double vm=0.;
    // double il=0.;
    // double so=0.;
    // double ca=0., ci=0., cs=0.;
    // double re=0;
    // for(int level=0;level<GetSolverPointers().size();level++)
    //   {
    // 	const StdSolver* S = dynamic_cast<const StdSolver*>(GetSolver(level));
    // 	assert(S);
    // 	vm += S->clock_vmult();
    // 	il += S->clock_ilu();
    // 	so += S->clock_solve();
    // 	ca += S->clock_computematrix();
    // 	ci += S->clock_computeilu();
    // 	cs += S->clock_computesolver();
    // 	re += S->clock_residual();
    //   }

    // cout << "  vmult\t\t\t" << vm << endl;
    // cout << "  smooth\t\t\t" << il << endl;
    // cout << "  solve\t\t\t" << so << endl;
    // cout << "  compute matrix\t" << ca << endl;
    // cout << "  compute ilu\t\t" << ci << endl;
    // cout << "  compute solver\t" << cs << endl;
    // cout << "VECTOR\t\t\tTIME\n";
    // cout << "  residual\t\t" << re << endl;
    // cout << "\n************************************************************************\n";
  }

  /*-------------------------------------------------------------*/

  StdMultiLevelSolver::StdMultiLevelSolver() : 
    _MAP(NULL), _cor("cor"), _res("res"), _mg0("mg0"), _mg1("mg1"),_trustdirection0("trustdirection0"),_trustdirection1("trustdirection1"),_trustdirection2("trustdirection2"),
    oldnlevels(-1), _paramfile(NULL), MON(NULL), DataP(NULL), _PD(0), _PC(0)
  {
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::BasicInit(const MeshAgentInterface* MAP, const ParamFile* paramfile,
				      const ProblemContainer* PC,
				      const FunctionalContainer* FC)
  {
    _MAP = MAP;

    _paramfile = paramfile;

    SetProblemContainer(PC);
    SetFunctionalContainer(FC);
  
    if(!DataP)
      {
	DataP = new StdMultiLevelSolverData;
      }
    DataP->BasicInit(_paramfile);
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::SetProblem(const std::string& label)
  {
    _PD = GetProblemContainer()->GetProblem(label);
    for(int level=0; level<nlevels(); ++level)  
      {
	int solverlevel = nlevels()-1-level;
	assert(GetSolver(solverlevel)) ;
	GetSolver(solverlevel)->SetProblem(*_PD);
      }
  }

  /*-------------------------------------------------------------*/

  StdSolver* StdMultiLevelSolver::NewSolver(int solverlevel) 
  { 
    if(DataP->Solver()=="instat")
      {
	return new StdTimeSolver;
      }
    else
      {
	return new StdSolver;
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::RegisterVectors() 
  {
    assert(nlevels()==GetSolverPointers().size());
    for (int level=0; level<nlevels(); ++level)  
      {
	GetSolver(level)->RegisterVector(_cor);
	GetSolver(level)->RegisterVector(_res);
	GetSolver(level)->RegisterVector(_mg0);
	GetSolver(level)->RegisterVector(_mg1);
	GetSolver(level)->RegisterVector(_trustdirection0);
	GetSolver(level)->RegisterVector(_trustdirection1);
	GetSolver(level)->RegisterVector(_trustdirection2);
      }
  }

  /*-------------------------------------------------------------*/


  void StdMultiLevelSolver::ReInitVector(VectorInterface& v)
  {
    for (int level=0; level<nlevels(); ++level)  
      {
	GetSolver(level)->ReInitVector(v);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::ReInitVector(VectorInterface& v, int comp)
  {
    for (int level=0; level<nlevels(); ++level)  
      {
	GetSolver(level)->ReInitVector(v,comp);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::NewSolvers()
  {
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
	    GetSolver(solverlevel)->BasicInit(_paramfile,GetMeshAgent()->GetDimension());
	  }
      }
  }

  /*-------------------------------------------------------------*/
  
  void StdMultiLevelSolver::RegisterMatrix()
  {
    for(int level=0; level<nlevels(); ++level)  
      {
	GetSolver(level)->RegisterMatrix();
      }
  }

  void StdMultiLevelSolver::ReInitMatrix()
  {
    for(int level=0; level<nlevels(); ++level)  
      {
	GetSolver(level)->ReInitMatrix();
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::SolverNewMesh()
  {
    for(int level=0; level<nlevels(); ++level)  
      {
	const GascoigneMesh* MIP = GetMeshAgent()->GetMesh(level);
	assert(MIP);

	int solverlevel = nlevels()-1-level;
	GetSolver(solverlevel)->NewMesh(MIP);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::NewMgInterpolator()
  {
    if (DataP->LinearSolve()=="direct") return; // No interpolator for direct solver
    
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
    //
    // Interpolator [l] :   interpoliert   (l+1)->l  (fein->grob)
    //
    for (int level=0;level<nlevels()-1;++level)
      {
	int sl = nlevels()-level-2;

	const MeshTransferInterface* MT = GetMeshAgent()->GetTransfer(sl);
	assert(MT);
	assert(_Interpolator[level]);
	GetSolver(level)->ConstructInterpolator(_Interpolator[level],MT);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::ReInit(const std::string& problemlabel)
  {
    DataP->CountResidual() = 0;
    // DataP->GetNLInfo().control().matrixmustbebuild() = 1;
  
    NewSolvers();
    SolverNewMesh();
    NewMgInterpolator();
    SetProblem(problemlabel);
    RegisterMatrix();
    RegisterVectors();
    ReInitMatrix();

    ReInitVector(_cor);
    ReInitVector(_res);
    ReInitVector(_mg0);
    ReInitVector(_mg1);
    ReInitVector(_trustdirection0);
    ReInitVector(_trustdirection1);
    ReInitVector(_trustdirection2);
  }

  /*-------------------------------------------------------------*/

  const DoubleVector StdMultiLevelSolver::GetExactValues() const
  {
    if (!GetFunctionalContainer()) return DoubleVector(0);
  
    int n = GetFunctionalContainer()->size();
    DoubleVector j(n,0.);
    int i = 0;
    for (FunctionalContainer::const_iterator it = GetFunctionalContainer()->begin();
	 it!=GetFunctionalContainer()->end();++it,++i)
      j[i] = it->second->ExactValue();
    return j;
  }

  /*-------------------------------------------------------------*/

  const DoubleVector StdMultiLevelSolver::ComputeFunctionals(VectorInterface& f, const VectorInterface& u,
							     FunctionalContainer* FC)
  {
    if (!FC) return DoubleVector(0);
    int n = FC->size();
    DoubleVector j(n,0.);
    int i = 0;
    for (FunctionalContainer::const_iterator it = FC->begin(); it!=FC->end();++it,++i)
      {
	j[i] = GetSolver(ComputeLevel)->ComputeFunctional(f,u,it->second);
      }
    return j;
  }

  /*-------------------------------------------------------------*/

  const DoubleVector StdMultiLevelSolver::ComputeFunctionals(VectorInterface& f, const VectorInterface& u)
  {
    if (!GetFunctionalContainer()) return DoubleVector(0);
    int n = GetFunctionalContainer()->size();
    DoubleVector j(n,0.);
    int i = 0;
    for (FunctionalContainer::const_iterator it = GetFunctionalContainer()->begin();
	 it!=GetFunctionalContainer()->end();++it,++i)
      {
	cout << it->first << " ";
	j[i] = GetSolver(ComputeLevel)->ComputeFunctional(f,u,it->second);
      }
    cout << endl;
    return j;
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::InterpolateSolution(VectorInterface& u, const GlobalVector& uold) const
  {
    if(oldnlevels<=0) return;

    GetSolver(FinestLevel())->InterpolateSolution(u,uold);

    for(int l=0; l<FinestLevel(); l++)
      {
	GetSolver(l)->Zero(u);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::vmulteq(VectorInterface& y, const VectorInterface& x) const
  {
    GetSolver(ComputeLevel)->vmulteq(y,x,1.);
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::LinearMg(int finelevel, int coarselevel, VectorInterface& u, const VectorInterface& f, CGInfo& info)
  {
    int clevel=coarselevel;
    for(int level=coarselevel;level<nlevels();level++)
      {
	if(GetSolver(level)->DirectSolver()) clevel=level;
      }
    // wir haben auf einem hoeheren level einen direkten loeser...
    if(clevel>finelevel)  {clevel=finelevel;}

    assert(finelevel>=clevel);

    int nl = nlevels();
    DoubleVector res(nl,0.), rw(nl,0.);

    GetSolver(finelevel)->Equ(_mg0,1.,f);

    GetSolver(finelevel)->MatrixResidual(_mg1,u,_mg0);
    res[finelevel] = GetSolver(finelevel)->Norm(_mg1);
    rw[finelevel] = 0.;
    info.check(res[finelevel],rw[finelevel]);

    bool reached = false; // mindestens einen schritt machen

    for(int it=0; !reached; it++)
      {
	string p = DataP->MgType();
	string p0 = p;
	if(p=="F") p0="W";
	mgstep(res,rw,finelevel,finelevel,clevel,p0,p,u,_mg0,_mg1);
	reached = info.check(res[finelevel],rw[finelevel]);
      }
  }

  /*-------------------------------------------------------------*/
 
  void StdMultiLevelSolver::mgstep(vector<double>& res, vector<double>& rw, 
				   int l, int finelevel, int coarselevel, string& p0, string p,
				   VectorInterface& u, VectorInterface& b, VectorInterface& v)
  {
    if(l==coarselevel)
      {
	if(p=="F") {p0="V";}
	GetSolver(l)->smooth_exact(u,b,v);
	if(l==finelevel)
	  {
	    GetSolver(l)->MatrixResidual(v, u, b);
	    res[l] = GetSolver(l)->Norm(v);
	  }
      }
    else
      {
	GetSolver(l)->smooth_pre(u,b,v);
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
    
	GetSolver(l)->smooth_post(u,b,v);
	GetSolver(l)->MatrixResidual(v,u,b);
	res[l] = GetSolver(l)->Norm(v);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::Cg(VectorInterface& x, const VectorInterface& f, CGInfo& info)
  {
    std::cerr << "\"StdMultiLevelSolver::Cg\" not written!" << std::endl;
    abort();
    //   CG<SolverInterface,StdMultiLevelSolver,GhostVector> cg(*(GetSolver(ComputeLevel)),*this);
  
    //   cg.solve(x,f,info);
  }


  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::newton(VectorInterface& u, const VectorInterface& f, VectorInterface& r, VectorInterface& w, NLInfo& info)
  {
  	double sigma;
    info.reset();
    double rr = NewtonResidual(r,u,f);
    bool reached = info.check(0,rr,0.);
    NewtonOutput(info);
    NewtonPreProcess(u,f,info);
    for(int it=1; !reached; it++)
      {
	NewtonMatrixControl(u,info);
	NewtonVectorZero(w);
	NewtonLinearSolve(w,r,info.GetLinearInfo());
	//if(it==1) sigma=min(GetSolver(ComputeLevel)->Norm(w),1.0e10);
	double rw = NewtonUpdate(rr,u,w,r,f,info);
	//double rw = NewtonUpdateTrustRegion(rr,u,w,r,f,info,sigma);
	reached = info.check(it,rr,rw);
	NewtonOutput(info);
      }
    NewtonPostProcess(u,f,info);
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::NewtonOutput(NLInfo& nlinfo) const
  {
    assert(MON!=0);
    MON->nonlinear_step(nlinfo.GetLinearInfo(),nlinfo);
  }
  /*-------------------------------------------------------------*/
 
  void StdMultiLevelSolver::NewtonPreProcess(VectorInterface& u, const VectorInterface& f,NLInfo& info) const  
  { ;
  }

  /*-------------------------------------------------------------*/
 
  void StdMultiLevelSolver::NewtonPostProcess(VectorInterface& u, const VectorInterface& f,NLInfo& info) const  
  { ;
  }

  /*-------------------------------------------------------------*/
 
  void StdMultiLevelSolver::NewtonVectorZero(VectorInterface& w) const
  {
    GetSolver(FinestLevel())->Zero(w);
  }

  /*-------------------------------------------------------------*/
 
  double StdMultiLevelSolver::NewtonResidual(VectorInterface& y, const VectorInterface& x,const VectorInterface& b) const
  {
    _clock_residual.start();
    DataP->CountResidual()++;
    GetSolver(ComputeLevel)->Equ(y,1.,b);
    GetSolver(ComputeLevel)->Form(y,x,-1.);
    GetSolver(ComputeLevel)->SetPeriodicVectorZero(y);
    GetSolver(ComputeLevel)->SetBoundaryVectorZero(y);
    GetSolver(ComputeLevel)->SubtractMeanAlgebraic(y);
    _clock_residual.stop();
    return NewtonNorm(y);
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::NewtonMatrixControl(VectorInterface& u, NLInfo& nlinfo)
  {
    MON->new_matrix() = 0;

    int nm1 = nlinfo.control().newmatrix();
    int nm2 = nlinfo.control().matrixmustbebuild();

    if (nm1+nm2==0) return;
  
    MON->new_matrix() = 1;
  
    AssembleMatrix(u,nlinfo);
    ComputeIlu(u);
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::NewtonLinearSolve(VectorInterface& x, const VectorInterface& b, CGInfo& info)
  {
    _clock_solve.start();
    info.reset();
    GetSolver(FinestLevel())->Zero(x);

    assert(DataP->LinearSolve()=="mg"  || DataP->LinearSolve()=="gmres" || DataP->LinearSolve()=="direct");

    //GlobalVector H = GetSolver()->GetGV(b);
    //dynamic_cast<const StdSolver*>(GetSolver())->GetMatrix()->vmult(H,GetSolver()->GetGV(x),-1.0);
    //double n1 = H.norm();

    if (DataP->LinearSolve()=="mg")
      {
	int clevel=std::max(DataP->CoarseLevel() ,0);
	if(DataP->CoarseLevel() == -1) clevel = FinestLevel(); 
	LinearMg(ComputeLevel,clevel,x,b, info);
      }
    else if (DataP->LinearSolve()=="direct")
      {
	int clevel = FinestLevel(); 
	LinearMg(ComputeLevel,clevel,x,b, info);
      }
    else if (DataP->LinearSolve()=="gmres")
      {
	Gmres(x,b,info);
      }

    //H = GetSolver()->GetGV(b);
    //dynamic_cast<const StdSolver*>(GetSolver())->GetMatrix()->vmult(H,GetSolver()->GetGV(x),-1.0);
    //double n2 = H.norm();
    //cout << "NORM: " << n1 << " " << n2 << "\t" << n2/n1 << endl;


    GetSolver(ComputeLevel)->SubtractMean(x);
    _clock_solve.stop();
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::Gmres(VectorInterface& x, const VectorInterface& f, CGInfo& info)
  {
    int n = DataP->GmresMemSize();

    StdSolver* S = dynamic_cast<StdSolver*>(GetSolver(ComputeLevel));
    GMRES<StdSolver,StdMultiLevelSolver,VectorInterface> gmres(*S,*this,n);

    gmres.solve(x,f,info);
  }

  /*-------------------------------------------------------------*/

  double StdMultiLevelSolver::NewtonUpdate(double& rr, VectorInterface& x, VectorInterface& dx, VectorInterface& r, const VectorInterface& f, NLInfo& nlinfo)
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

    GetSolver(ComputeLevel)->Add(x,relax,dx);
    NewtonResidual(r,x,f);
    rr = NewtonNorm(r);

    string message = "";
    for(int iter=0;iter<nlinfo.user().maxrelax();iter++)
      {
	message = nlinfo.check_damping(iter,rr);

	if (message=="ok")       break;
	if (message=="continue") 
	  {
	    GetSolver(ComputeLevel)->Add(x,relax*(omega-1.),dx);  

	    NewtonResidual(r,x,f);
	    rr = NewtonNorm(r);
	    relax *= omega;
	    continue;
	  }
	if (message=="exploded")
	  {
	    //GetSolver(ComputeLevel)->Add(x,-relax,dx);
	    //relax = 0.;
	    //cout << "Damping exploded !!!!!" << endl;
	    //nlinfo.control().status() = "diverged";
	    //break;
	    cout << "Damping exploded !!!!!" << endl;
	    //nlinfo.control().status() = "diverged";
	    GetSolver(ComputeLevel)->Add(x,relax*(omega-1.),dx);  

	    NewtonResidual(r,x,f);
	    rr = NewtonNorm(r);
	    relax *= omega;
	    continue;
	  }
      }

    // NewtonUpdateShowCompResiduals(nlinfo.control().iteration(), x, r, f,dx);

    return NewtonNorm(dx);
  }

  double StdMultiLevelSolver::NewtonUpdateTrustRegion(double& rr, VectorInterface& x, VectorInterface& dx, VectorInterface& r, const VectorInterface& f, NLInfo& nlinfo,double& sigma)
  {
    cout<<"===================================================="<<endl;
    //const CGInfo& linfo = nlinfo.GetLinearInfo();
    double nn = NewtonNorm(dx);
    double nr = GetSolver(ComputeLevel)->Norm(r);
    
    bool trust_failed=true;
    
    
    GetSolver(ComputeLevel)->vmulteq(_trustdirection0,dx,1.0);
    GetSolver(ComputeLevel)->Add(_trustdirection0,-1.0,r);


    while(trust_failed)
    {
		double tau=0;
		double tau1=0;
		double tau2=0;
		GetSolver(ComputeLevel)->Equ(r,-1.0,r);
		//cout<<"Norm(dx)"<<GetSolver(ComputeLevel)->Norm(dx)<<endl;
		/*if(GetSolver(ComputeLevel)->Norm(dx)<sigma  )
		{
		  cout<<"Newton ";
		  GetSolver(ComputeLevel)->Equ(_trustdirection0,1.0,dx);
		}
		else*/
		{
			//Cauchy Point: _trustdirection0
			//Newton direction: dx
		
			GetSolver(ComputeLevel)->vmult_transposedeq(_trustdirection0,r,1.0);
			double nJTxr = GetSolver(ComputeLevel)->Norm(_trustdirection0);

			GetSolver(ComputeLevel)->vmulteq(_trustdirection1,_trustdirection0,1.0);
		
			//GetSolver(ComputeLevel)->vmult_transposedeq(_trustdirection0,_trustdirection1,1.0);
			//GetSolver(ComputeLevel)->vmulteq(_trustdirection1,_trustdirection0,1.0);
			//double tau= std::min(1.0,pow(nJTxr,3)/(sigma*GetSolver(ComputeLevel)->ScalarProduct(r,_trustdirection1)));

			tau= std::min(1.0,pow(nJTxr,3)/(sigma*GetSolver(ComputeLevel)->ScalarProduct(_trustdirection1,_trustdirection1)));	
	
			//GetSolver(ComputeLevel)->vmult_transposedeq(_trustdirection0,r,1.0);
			GetSolver(ComputeLevel)->Equ(_trustdirection0,-tau*sigma/nJTxr,_trustdirection0);
			//cout<<"tau: "<<tau<<endl;
			if(tau!=1 )
			{ //Cauchy Punkt nicht auf Rand der TrustRegion

				double n1=GetSolver(ComputeLevel)->ScalarProduct(_trustdirection0,_trustdirection0);
				double n2=GetSolver(ComputeLevel)->ScalarProduct(_trustdirection0,dx);
				double n3=GetSolver(ComputeLevel)->ScalarProduct(dx,dx);
				
				tau1=((2.0*n1-2.0*n2)+sqrt((2.0*n1-2.0*n2)*(2.0*n1-2.0*n2)-4.0*(n1-2.0*n2+n3)*(n1-sigma*sigma)))/(2.0*(n1-2.0*n2+n3));
				tau2=((2.0*n1-2.0*n2)-sqrt((2.0*n1-2.0*n2)*(2.0*n1-2.0*n2)-4.0*(n1-2.0*n2+n3)*(n1-sigma*sigma)))/(2.0*(n1-2.0*n2+n3));
				
				
				// Cauchy Point _trustdirection0
				// dx Newton Point
				
				//Residuums vector r(x)+J*p_cauchy
			   	GetSolver(ComputeLevel)->vmulteq(_trustdirection1,_trustdirection0,1.0);
				GetSolver(ComputeLevel)->Add(_trustdirection1,1.0,r);
				//double r_model_cauchy=GetSolver(ComputeLevel)->Norm(_trustdirection1);
				//cout<<"Model Norm Cauchy Residuum: "<<r_model_cauchy<<endl;
								
				//Residuums vector r(x)+J*p_Newt
			   	GetSolver(ComputeLevel)->vmulteq(_trustdirection2,dx,1.0);
				GetSolver(ComputeLevel)->Add(_trustdirection2,1.0,r);
				//double r_model_newt=GetSolver(ComputeLevel)->Norm(_trustdirection2);
				//cout<<"Model Norm Newton Residuum: "<<r_model_newt<<endl;
				

				GetSolver(ComputeLevel)->Add(_trustdirection2,-1.0,_trustdirection1);
				
				n1= -GetSolver(ComputeLevel)->ScalarProduct(_trustdirection1,_trustdirection2);
				n2= GetSolver(ComputeLevel)->ScalarProduct(_trustdirection2,_trustdirection2);
				
				double gamma = n1/n2;
				//if(gamma<0.1) gamma=0;
				double gamma_min=min(tau1,tau2);
				double gamma_plus=max(tau1,tau2);
				//if(max(tau1,tau2)<0 && min(tau1,tau2)>1.0)
				//{
				//	cout<< "error in Trust Region"<<endl; abort();
				//}
				//else 
				//{
					GetSolver(ComputeLevel)->Equ(_trustdirection1,1.0,dx);
					GetSolver(ComputeLevel)->Add(_trustdirection1,-1.0,_trustdirection0);

					//GetSolver(ComputeLevel)->Add(_trustdirection0,max(tau1,tau2),_trustdirection1);
					GetSolver(ComputeLevel)->Add(_trustdirection0,max(gamma_min,min(gamma,gamma_plus)),_trustdirection1);
					cout<<"gamma_min: "<<gamma_min<<" gamma: "<<gamma<<" gamma_plus: "<<gamma_plus<<"used gamma: "<<max(gamma_min,min(gamma,gamma_plus))<<endl;
				//}
				//cout<<"n1: "<<n1<<" n2: "<<n2<<" n3: "<<n3<<" ||tr_dir||: "<<GetSolver(ComputeLevel)->Norm(_trustdirection0)<<endl;
			}
		}
		//Update TrustRegion
		//||r(x)||
		double r_old=GetSolver(ComputeLevel)->Norm(r);
		double r_old_infty=NewtonNorm(r);
		//||r(x)+J*p||
	   	GetSolver(ComputeLevel)->vmulteq(_trustdirection1,_trustdirection0,1.0);
		GetSolver(ComputeLevel)->Add(r,1.0,_trustdirection1);
		
		double r_model=GetSolver(ComputeLevel)->Norm(r);
		double r_model_infty=NewtonNorm(r);
		
		//Update 		
		GetSolver(ComputeLevel)->SetPeriodicVectorZero(_trustdirection0);
		GetSolver(ComputeLevel)->Add(x,1.0,_trustdirection0);
		//||r(x+p)||
		NewtonResidual(r,x,f);
		double r_new = GetSolver(ComputeLevel)->Norm(r);  
		double r_new_infty = NewtonNorm(r);  
		//NewtonNorm = l_infty norm!!!!
		rr= NewtonNorm(r);
		
		//Check auf Abstieg:
		string message = nlinfo.check_damping(0,rr);
		
		

		//cout<<"r_old*r_old-r_model*r_model: "<<r_old*r_old-r_model*r_model<<endl;
		//cout<<"r_old*r_old-r_new*r_new: "<<r_old*r_old-r_new*r_new<<endl;
		
		double rho=(r_old-r_new)/(r_old-r_model);
		//double rho_infty=(r_old_infty*r_old_infty-r_new_infty*r_new_infty)/(r_old_infty*r_old_infty-r_model_infty*r_model_infty);
		
		  double sigma_max=1.0e10;
		  double sigma_min=1.0e-6;
		  double eps = 0.25;
		  double k1  = 2.0;
		  double k2  = 0.25;

		  double c=fabs(rho-1);
		  //double c_infty=fabs(rho_infty-1);
		  if((r_old*r_old-r_model*r_model<=0) && (r_old*r_old-r_new*r_new<=0))
			cout<<"no reduction in model!!! Singularity in Matrix??????"<<endl;
		  
		  if(rho<=0.0 && sigma>sigma_min)
		  {
		  	
		  	sigma*=k2*k2;
		  	sigma=max(sigma,sigma_min);
			//cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cout<<" !Reject ! ";
			//cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			GetSolver(ComputeLevel)->Add(x,-1.0,_trustdirection0);
			NewtonResidual(r,x,f);
			rr = NewtonNorm(r);
			message = nlinfo.check_damping(0,rr);

			
		  }
		  else
		  {
		  	trust_failed=false;	
			  if(c<eps/8.)
			  {
				if(k1*k1*sigma>sigma_max)
				{
				  sigma=sigma_max;
				}
				else
				{
				  sigma*=k1*k1;
				}
			  }
			  else if (rho>1.0-eps )
			  {
				if(k1*sigma>sigma_max)
				{
				  sigma=sigma_max;
				}
				else
				{
				  sigma*=k1;
				}
			  }
			  else if (rho>0.5 )
			  {
				if(1.2*sigma>sigma_max)
				{
				  sigma=sigma_max;
				}
				else
				{
				  sigma*=1.2;
				}
			  }
			  else if (rho<0.1 )
			  {
				sigma*=k2;
				sigma=max(sigma,sigma_min);
			  }
			
		}
		cout<<" rho: "<<rho<<" sigma_new: "<<sigma<<"||r||_2: "<<r_new<<endl;
	}

    return NewtonNorm(_trustdirection0);
  }
  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::AssembleMatrix(VectorInterface& u)
  {
    SolutionTransfer(u);
    for(int l=0;l<=ComputeLevel;l++)
      {
	GetSolver(l)->MatrixZero();
	GetSolver(l)->AssembleMatrix(u,1.);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::AssembleMatrix(VectorInterface& u, NLInfo& nlinfo)
  {
    AssembleMatrix(u);
    nlinfo.control().matrixmustbebuild() = 0;
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::ComputeIlu()
  {
    for(int l=0;l<=ComputeLevel;l++)
      {
	GetSolver(l)->ComputeIlu();
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::ComputeIlu(VectorInterface& u)
  {
    SolutionTransfer(u);
    for(int l=0;l<=ComputeLevel;l++)
      {
	GetSolver(l)->ComputeIlu(u);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::BoundaryInit(VectorInterface& u) const
  {
    for(int l=0;l<GetSolverPointers().size();l++)
      GetSolver(l)->BoundaryInit(u);
  }

  /*-------------------------------------------------------------*/

  string StdMultiLevelSolver::LinearSolve(int level, VectorInterface& u, const VectorInterface& b, CGInfo& info)
  {
    ComputeLevel = level;

    GetSolver(ComputeLevel)->HNAverage(u);
  
    info.reset();
  
    int clevel=std::max(DataP->CoarseLevel() ,0);
    if(DataP->CoarseLevel() == -1) clevel = FinestLevel(); 

    LinearMg(ComputeLevel,clevel,u,b,info);

    GetSolver(ComputeLevel)->SubtractMean(u);

    GetSolver(ComputeLevel)->HNZero(u);
    string status = info.control().status();

    return status;
  }

  /*-------------------------------------------------------------*/

  string StdMultiLevelSolver::Solve(int level, VectorInterface& u, const VectorInterface& b, NLInfo& nlinfo)
  {
    ComputeLevel = level;

    string status;

    assert(DataP->NonLinearSolve() == "newton");
    //   if(DataP->NonLinearSolve() == "newton")
    //     {
    GetSolver(ComputeLevel)->HNAverage(u);
    newton(u,b,_res,_cor,nlinfo);
    GetSolver(ComputeLevel)->HNZero(u);
    return nlinfo.CheckMatrix();
    //     }
  }

  /*-------------------------------------------------------------*/

  double StdMultiLevelSolver::ComputeFunctional(VectorInterface& f, const VectorInterface& u, const std::string& label)
  {
    const Functional* FP = GetFunctionalContainer()->GetFunctional(label);
    return GetSolver(ComputeLevel)->ComputeFunctional(f,u,FP);
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::Transfer(int high, int low, VectorInterface& u) const
  {
    if (DataP->LinearSolve()=="direct") return; // No transfer for direct solver
    
    for(int l=high;l>=low;l--)
      {
	GetSolver(l)->HNAverage(u);

	assert(_Interpolator[l-1]);
	_Interpolator[l-1]->SolutionTransfer(GetSolver(l-1)->GetGV(u),
					     GetSolver(l)  ->GetGV(u));
      
	GetSolver(l)->HNZero(u);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::Transfer(VectorInterface& u) const
  {
    Transfer(ComputeLevel,1,u);
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::SolutionTransfer(int high, int low, VectorInterface& u) const
  {
    if (DataP->LinearSolve()=="direct") return; // No transfer for direct solver

    Transfer(high,low,u);

    for(int l=high;l>=low;l--)
      GetSolver(l-1)->SetBoundaryVector(u);
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::SolutionTransfer(VectorInterface& u) const
  {
    SolutionTransfer(ComputeLevel,1,u);
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::AssembleDualMatrix(VectorInterface& u)
  {
    for(int l=0; l<nlevels(); l++)
      {
	GetSolver(l)->AssembleDualMatrix(u,1.);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::DeleteVector(VectorInterface& v)
  {
    for(int l=0; l<nlevels(); ++l)  
      {
	GetSolver(l)->DeleteVector(v);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::precondition(VectorInterface& x, VectorInterface& y)
  {
    CGInfo& precinfo = DataP->GetPrecInfo();
    precinfo.reset();
    precinfo.check(0.,0.);

    int clevel=std::max(DataP->CoarseLevel(),0);
    if(DataP->CoarseLevel() == -1) clevel = FinestLevel(); 

    LinearMg(ComputeLevel,clevel,x,y,precinfo);
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::Equ(VectorInterface& dst, double s, const VectorInterface& src) const
  {
    for(int l=0; l<nlevels(); l++)
      {
	GetSolver(l)->Equ(dst,s,src);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::Zero(VectorInterface& dst) const
  {
    for(int l=0; l<nlevels(); l++)
      {
	GetSolver(l)->Zero(dst);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::AddNodeVector(const string& name, VectorInterface& gq)
  {
    Transfer(ComputeLevel,1,gq);
    for(int l=0; l<nlevels(); l++)
      {
	GetSolver(l)->AddNodeVector(name,gq);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::DeleteNodeVector(const string& name)
  {
    for(int l=0; l<nlevels(); l++)
      {
	GetSolver(l)->DeleteNodeVector(name);
      }
  }

  /*-------------------------------------------------------------*/

  void StdMultiLevelSolver::InterpolateCellSolution(VectorInterface& u, const GlobalVector& uold) const
  {
    assert(u.GetType()=="cell");
    GlobalVector &uu = GetSolver()->GetGV(u);

    if(!GetMeshAgent()->Goc2nc())
      {
	cerr << "No cell interpolation activated"<< endl;
	abort();
      }
    int cells = GetMeshAgent()->ncells();
    uu.ReInit(uold.ncomp(),cells);

    uu.zero();
    //Jetzt Interpolieren wir die Loesung auf das neue Gitter
    for(int i = 0; i < uold.n(); i++)
      {
	set<int> kinder = GetMeshAgent()->Cello2n(i);
	if(!kinder.empty())
	  {
	    //es wurde nicht vergroebert
	    for(set<int>::iterator p = kinder.begin(); p != kinder.end(); p++)
	      {
		for(int c=0; c<uold.ncomp();++c)
		  {
		    uu(*p,c) = uold(i,c);
		  }
	      }
	  }
	else
	  {
	    //Es wurde vergroebert
	    int patchsize = 1<<GetMeshAgent()->GetMesh(FinestLevel())->dimension();
	    uu[GetMeshAgent()->Cello2nFather(i)] += uold[i]/static_cast<double>(patchsize);
	  }

      }
  }

  /*-------------------------------------------------------------*/
}
