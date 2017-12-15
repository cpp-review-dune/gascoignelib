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

  class MyPeriodicMapping : public virtual PeriodicMapping
  {
    public:
    std::string GetName() const { return "Periodic Data"; }
    
    void transformCoords(Vertex2d& w, const Vertex2d& v) const
    {
      w.x()=v.x()-40;
      w.y()=v.y();
    }
    void transformCoords(Vertex3d& w, const Vertex3d& v) const
    {
      abort();
    }

  };
  

  
  

  class Loop : public StdLoop
  {
    SimpleMatrix    MM;
    nvector<double> LMM;
    
    MyPeriodicMapping mymap;

    int _M;
    
  public:

    string SolvePrimalSingle(VectorInterface& u, VectorInterface& f, string name);
    void SolvePrimalProblem(vector<GlobalVector> &Utotal, nvector<double>& Jtotal, VectorInterface& u, VectorInterface& oldu, VectorInterface& f);
    void BoxInt(GlobalVector& avg, const vector<GlobalVector>& U, int start, int stopp);
    void TrapezInt(GlobalVector& avg, const vector<GlobalVector>& U, int start, int stopp);
    void MittelInt(GlobalVector& avg_old,GlobalVector& avg, const vector<GlobalVector>& U, int start, int stopp);
    void Gauss_Q2(GlobalVector& avg_old,GlobalVector& avg, const vector<GlobalVector>& U, int start, int stopp);
    
    
    
    // Fehlerschaetzer
    void SolveDualProblem(VectorInterface& z, VectorInterface& f,double DTM);
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
			   vector<GlobalVector>& Ztotal,
			   vector<GlobalVector>& Pu_kM,
			   vector<GlobalVector>& Pu_M,
			   VectorInterface& u,
			   VectorInterface& oldu,
			   VectorInterface& z,
			   VectorInterface& f);

    void EstimateDWRprim(DoubleVector& eta, int m, const GlobalVector& Pu_kM,   vector<GlobalVector>& U,
		     GlobalVector& Z,VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);
   void EstimateDWRdual(DoubleVector& eta, int m, const GlobalVector& Pu_kM, GlobalVector& Pu_M,
			const GlobalVector& OLDZ, GlobalVector& Z,VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);
    void EstimateInitial(DoubleVector& eta, const GlobalVector& U, GlobalVector &Z,
			 VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);
    void EstimateAvg(DoubleVector& eta,  GlobalVector& Pu, const GlobalVector &Puold,
		     const GlobalVector& U, const GlobalVector &Uold, GlobalVector& Z,
		     VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);
    void EstimateRest(DoubleVector& eta,int m,
		     const GlobalVector& Pu, const GlobalVector &Puold,
		     const GlobalVector& Puk, const GlobalVector &Pukold,
		     const GlobalVector& U, const GlobalVector &Uold, GlobalVector& Z,
		     VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);

    void EstimateNonU(DoubleVector& eta,int m,
		       vector<GlobalVector>& Utotal, GlobalVector& Z,
		 VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);

 void EstimateNonPu(DoubleVector& eta,
		       const GlobalVector& U, GlobalVector& Z,
		 VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);
 
void EstimateNonMeanPu(DoubleVector& eta,int m,
		 GlobalVector& Pu, vector<GlobalVector>& U, GlobalVector& Z,
		 VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);


 void EstimateNonMeanU(DoubleVector& eta, int m,
		    GlobalVector& Pu,  vector<GlobalVector>& U, GlobalVector& Z,
		 VectorInterface& u, VectorInterface& oldu,VectorInterface& z,VectorInterface& f);
 
 


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
      GetMultiLevelSolverPointer() = new MyMLS;
      
      StdLoop::BasicInit(paramfile, PC, FC);
    }
    
   string PrecondCGMass(GlobalVector& u, GlobalVector& f, const TimePattern& TP, double s);
   void run(const std::string& problemlabel);
   void runwater(const std::string& problemlabel);
   void AssembleMassMatrix();
   void SolveDIV(VectorInterface& div,VectorInterface& vel,VectorInterface& f);
   void SolveTransport(VectorInterface& h, VectorInterface& hl, VectorInterface& hh, 
		       VectorInterface& div, VectorInterface& vel, VectorInterface& f);
  
  };
  
}



/*----------------------------   loop.h     ---------------------------*/
#endif
/*----------------------------   loop.h     ---------------------------*/
