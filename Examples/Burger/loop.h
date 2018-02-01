/*----------------------------   loop.h     ---------------------------*/
#ifndef __loop_H
#define __loop_H
/*----------------------------   loop.h     ---------------------------*/


#include "stdloop.h"
#include  "periodicmapping.h"
#include "meshagent.h"
#include "gascoignemesh2d.h"
#include  "stdmultilevelsolver.h"
#include  "stdsolver.h"
#include  "simplematrix.h"
#include  "sparseblockmatrix.h"
#include  "fmatrixblock.h"

#include "solvers.h"


namespace Gascoigne
{

  
  

  class Loop : public StdLoop
  {
    
    
 

    int _M;
    
  public:

    string SolvePrimalSingle(VectorInterface& u, VectorInterface& f, string name);
    void SolvePrimalProblem(vector<GlobalVector> &Utotal, nvector<double>& Jtotal, VectorInterface& u, VectorInterface& oldu, VectorInterface& f, int ADAITER);

    void Reconstruction(GlobalVector& U2, const vector<GlobalVector>& U,int start);

    void TrapezInt(GlobalVector& avg, const vector<GlobalVector>& U, int start, int stopp,int DTM );
    void MittelInt(GlobalVector& avg_old,GlobalVector& avg, const vector<GlobalVector>& U, double start, double stopp,double DTM);
    void Gauss_Q2(GlobalVector& avg_old,GlobalVector& avg, const vector<GlobalVector>& U, int start, int stopp,int DTM) ;
 
  
   
    
    void SolveDualProblem(vector<GlobalVector>& Ztotal, VectorInterface& f, VectorInterface& u,  VectorInterface& oldu, VectorInterface& newu, VectorInterface& z,  VectorInterface& oldz,const vector<GlobalVector>& Pu_k,int ADAITER, vector<double>& DT_M);


 void EstimateDualError(DoubleVector& eta,
			    DoubleVector& eta0,
			     DoubleVector& eta1,
			     DoubleVector& eta11,
			     DoubleVector& eta2,
			     DoubleVector& eta22,
			     DoubleVector& eta23,
			     DoubleVector& eta3,
			   DoubleVector& eta4,
			      DoubleVector& eta5,
			   vector<GlobalVector>& Utotal,
			vector<GlobalVector>& U_2,
			   vector<GlobalVector>& Ztotal,
			   vector<GlobalVector>& Pu_kM,
			   vector<GlobalVector>& Pu_M,
			VectorInterface& u,
		   	   VectorInterface& oldu,
			VectorInterface& newu,
			   VectorInterface& z,
			  VectorInterface& oldz,
			   VectorInterface& f, vector<double>& DT_M,vector<double>& T);

    void EstimateDWRprim(DoubleVector& eta, int m, const GlobalVector& Pu_kM,   vector<GlobalVector>& U,
		     GlobalVector& Z,VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);

     void EstimateDWRdual(DoubleVector& eta, int m, vector<GlobalVector>& Pu_kM, GlobalVector& Pu_M,
			const GlobalVector& OLDZ, GlobalVector& Z,VectorInterface& u, VectorInterface& oldu,VectorInterface& newu,VectorInterface& z,VectorInterface& oldz,VectorInterface& f, vector<double>& DT_M);

   void EstimateAvg(DoubleVector& eta,  GlobalVector& Pu, const GlobalVector &Puold,
		     const GlobalVector& U, const GlobalVector &Uold, GlobalVector& Z,
		     VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);


   void EstimateRest(DoubleVector& eta,int m,
		     const GlobalVector& Pu, const GlobalVector &Puold,
		     const GlobalVector& Puk, const GlobalVector &Pukold,
		     const GlobalVector& U, const GlobalVector &Uold, GlobalVector& Z,
		     VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);

    void EstimateNonU(DoubleVector& eta,
		       vector<GlobalVector>& Utotal, GlobalVector& Z,
		 VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f, int start, int stopp);
    void EstimateNonPu(DoubleVector& eta, 
		      vector<GlobalVector>& Pu_k, GlobalVector& Z,
		       VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f,int m, double DTM_PU);

    void EstimateNonMeanU(DoubleVector& eta, int m,
		    GlobalVector& Pu, GlobalVector& Pu_k,vector<GlobalVector>& U,vector<GlobalVector>& U_2, GlobalVector& Z,
		 VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f,int start, int stopp);
    
    void EstimateNonMeanPu(DoubleVector& eta,int m,
		 GlobalVector& Pu,vector<GlobalVector>& Pu_k, vector<GlobalVector>& U, GlobalVector& Z,
		 VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f,double DTM_U);


 MySolver* GetMySolver()
    {
      assert(dynamic_cast<MySolver*> (GetMultiLevelSolver()->GetSolver()));
      return dynamic_cast<MySolver*> (GetMultiLevelSolver()->GetSolver());
    }
    const MySolver* GetMySolver() const
    {
      assert(dynamic_cast<const MySolver*> (GetMultiLevelSolver()->GetSolver()));
      return dynamic_cast<const MySolver*> (GetMultiLevelSolver()->GetSolver());
    }
      
    
    void BasicInit(const ParamFile* paramfile, const ProblemContainer* PC, const FunctionalContainer* FC)
    {
      GetMeshAgentPointer() = new MeshAgent();
      GetMultiLevelSolverPointer() = new MyMLS();
      
      StdLoop::BasicInit(paramfile, PC, FC);
    }
    
   void run(const std::string& problemlabel); 
  
  };
  
}



/*----------------------------   loop.h     ---------------------------*/
#endif
/*----------------------------   loop.h     ---------------------------*/
