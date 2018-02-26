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

    double CompFunctional(vector<GlobalVector> &Utota,VectorInterface& u,VectorInterface& f);

    string SolveTransportSingle(VectorInterface& h, VectorInterface& f, string name);
    string SolvePrimalSingle(VectorInterface& u, VectorInterface& f, string name);
    void SolvePrimalProblem(vector<GlobalVector> &Utotal, VectorInterface& u, VectorInterface& oldu,VectorInterface& newu,VectorInterface& dtu3, VectorInterface& f,vector<GlobalVector> &Htotal, VectorInterface& h, VectorInterface& oldh, int ADAITER);

    
    void MittelInt(GlobalVector& avg_old,GlobalVector& avg, const vector<GlobalVector>& U, int  start, int stopp,double DTM);
    
  
   void SolveDualProblem(vector<GlobalVector>& Ztotal, vector<GlobalVector>& Wtotal,vector<GlobalVector>& Htotal,VectorInterface& f, VectorInterface& u,  VectorInterface& oldu, VectorInterface& newu,VectorInterface& dtu3, VectorInterface& z,  VectorInterface& oldz,VectorInterface& oldoldz,VectorInterface& w,VectorInterface& oldw,VectorInterface& h,VectorInterface& oldh,VectorInterface& newh,const vector<GlobalVector>& Pu_k,int ADAITER, vector<double>& DT_M,vector<double>& T);
    
string SolveDualSingle(vector<GlobalVector>& Ztotal,vector<GlobalVector>& Wtotal,vector<GlobalVector>& Htotal,VectorInterface& f, VectorInterface& u,  VectorInterface& oldu, VectorInterface& newu,VectorInterface& dtu3, VectorInterface& z,  VectorInterface& oldz,VectorInterface& oldoldz,VectorInterface& w,VectorInterface& oldw,VectorInterface& h,VectorInterface& oldh,VectorInterface& newh,const vector<GlobalVector>& Pu_k,int m, vector<double>& DT_M,vector<double>& T,string name);


 string SolveDualTransportSingle(vector<GlobalVector>& Ztotal,vector<GlobalVector>& Wtotal,vector<GlobalVector>& Htotal,VectorInterface& f, VectorInterface& u,  VectorInterface& oldu, VectorInterface& newu,VectorInterface& dtu3, VectorInterface& z,  VectorInterface& oldz,VectorInterface& w,VectorInterface& oldw,VectorInterface& h,VectorInterface& oldh,const vector<GlobalVector>& Pu_k,int m, vector<double>& DT_M,vector<double>& T,string name);


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
			DoubleVector& eta6,
			DoubleVector& eta7,
			   vector<GlobalVector>& Utotal,
               vector<GlobalVector>& Htotal,
               vector<GlobalVector>& Wtotal,
			vector<GlobalVector>& U_2,
			   vector<GlobalVector>& Ztotal,
			   vector<GlobalVector>& Pu_kM,
			   vector<GlobalVector>& Pu_M,
               vector<GlobalVector>& Ph_M,
			VectorInterface& u,
		   	   VectorInterface& oldu,
			VectorInterface& newu,
			VectorInterface& dtu3,
			   VectorInterface& z,
			  VectorInterface& oldz,
             VectorInterface& h ,
              VectorInterface& oldh ,
               VectorInterface& w ,
                VectorInterface& oldw ,
			VectorInterface& f, vector<double>& DT_M,vector<double>& T,vector<double>& eta_time);
 
 
  
    void EstimateDWRprim(DoubleVector& eta, int m, const GlobalVector& Pu_kM,  vector<GlobalVector>& U,vector<GlobalVector>& H,
			 GlobalVector& Z,VectorInterface& u, VectorInterface& oldu,VectorInterface& h,VectorInterface& oldh,VectorInterface& z,VectorInterface& f,vector<double>& T);

     void EstimateDWRdual(DoubleVector& eta, int m, vector<GlobalVector>& Pu_kM,vector<GlobalVector>&Htotal,vector<GlobalVector>&Wtotal, GlobalVector& Pu_M,
			  const GlobalVector& OLDZ, GlobalVector& Z,VectorInterface& u, VectorInterface& oldu,VectorInterface& newu,VectorInterface& dtu3,VectorInterface& z,VectorInterface& oldz,VectorInterface& h,VectorInterface& oldh,VectorInterface& w,VectorInterface& oldw,VectorInterface& f, vector<double>& DT_M,vector<double>& T);


    
void EstimateAvg(DoubleVector& eta, GlobalVector& Pu, const GlobalVector &Puold,GlobalVector& Ph, const GlobalVector &Phold,
		       const GlobalVector& U, const GlobalVector &Uold, GlobalVector& Z,GlobalVector& W, GlobalVector& H,GlobalVector& OLDH,
		 VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& h,VectorInterface& f);
     
  


 void EstimateRest(DoubleVector& eta, int m,
			  const GlobalVector& Pu, const GlobalVector &Puold,
			  const GlobalVector& Puk, const GlobalVector &Pukold,
			  const GlobalVector& U, const GlobalVector &Uold, GlobalVector& Z,GlobalVector& W,const GlobalVector &H,const GlobalVector &OLD,const GlobalVector &Phk,const GlobalVector &Phkold,
		   VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& h,VectorInterface& oldh,VectorInterface& f);


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

 void  EstimateKonsistenz(DoubleVector& eta, int m,
			      vector<GlobalVector>& U,vector<GlobalVector>& H,  GlobalVector& W,
		     VectorInterface& u, VectorInterface& oldu,VectorInterface& h,VectorInterface& z,VectorInterface& f, int start , int stopp);
    


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
