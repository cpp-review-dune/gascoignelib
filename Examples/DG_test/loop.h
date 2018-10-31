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
#include "umfilu.h"
#include "pointmatrix.h"




namespace Gascoigne
{

  
  class MySO : public StdSolver
  {
    MatrixInterface  *A1, *A2, *A3, *A4;
    IluInterface     *I1, *I2, *I3, *I4;

  public:

  MySO() : StdSolver(), A1(NULL),A2(NULL),A3(NULL),A4(NULL),
      I1(NULL),I2(NULL),I3(NULL),I4(NULL)
      {}
    MatrixInterface* GetMatrix() const
    {
      if (GetProblemDescriptor()->GetName() == "LaplaceTProblem")
	return A1;
      else if (GetProblemDescriptor()->GetName() == "Dual  Problem")
	return A2;
      else if (GetProblemDescriptor()->GetName() == "BurgerProblem")
	return A3;
      else if (GetProblemDescriptor()->GetName() == "BurgerDualProblem")
	return A4;
      else abort();
    }
    IluInterface* GetIlu() const
    {
      if (GetProblemDescriptor()->GetName() == "LaplaceTProblem")
	return I1;
      else if (GetProblemDescriptor()->GetName() == "Dual  Problem")
	return I2;
      else if (GetProblemDescriptor()->GetName() == "BurgerProblem")
	return I3;
      else if (GetProblemDescriptor()->GetName() == "BurgerDualProblem")
	return I4;
      else abort();
    }
    MatrixInterface*& GetMatrixPointer() {abort();}
    IluInterface*& GetIluPointer() { abort();}

    void RegisterMatrix()
    {
      const Equation*  EQ = GetProblemDescriptor()->GetEquation();
      assert(EQ);
      int ncomp = EQ->GetNcomp();

      if (A1) delete A1;
      if (A2) delete A2;
      if (A3) delete A3;
      if (A4) delete A4;
      
      if (A1==NULL) A1 = new PointMatrix(ncomp,"node");
      if (A2==NULL) A2 = new PointMatrix(ncomp,"node");
      if (A3==NULL) A3 = new PointMatrix(ncomp,"node");
      if (A4==NULL) A4 = new PointMatrix(ncomp,"node");

      if (I1) delete I1;
      if (I2) delete I2;
      if (I3) delete I3;
      if (I4) delete I4;
      
      if (I1==NULL) I1 = new UmfIlu(A1); 
      if (I2==NULL) I2 = new UmfIlu(A2); 
      if (I3==NULL) I3 = new UmfIlu(A3); 
      if (I4==NULL) I4 = new UmfIlu(A4); 
    }
    void ReInitMatrix() 
    {
      GetDiscretization()->InitFilter(GetPfilter());
      SparseStructure SA;
      GetDiscretization()->Structure(&SA);
      
      if (GetFaceDiscretization())
	GetFaceDiscretization()->Structure(&SA);
      
      AddPeriodicNodes(&SA);
      
      A1->ReInit(&SA);      I1->ReInit(&SA);
      A2->ReInit(&SA);      I2->ReInit(&SA);
      A3->ReInit(&SA);      I3->ReInit(&SA);
      A4->ReInit(&SA);      I4->ReInit(&SA);
    }
  };

  class MyMLSO : public StdMultiLevelSolver
  {
  public:
    SolverInterface* NewSolver(int solverlevel)
    { return new MySO; }
  };

  class Loop : public StdLoop
  {
    int _M;
    
  public:

    void BasicInit(const ParamFile* paramfile,
		   const ProblemContainer* PC,
		   const FunctionalContainer* FC=NULL)
    {
      //      GetMultiLevelSolverPointer() = new MyMLSO;
      StdLoop::BasicInit(paramfile,PC,FC);
    }
      

    string SolvePrimalSingle(VectorInterface& u, VectorInterface& f, string name);
  
    string SolvePrimalBurgerSingle(VectorInterface& v, VectorInterface& f, string name);
    
    void SolvePrimalProblem(vector<GlobalVector> &Utotal,vector<GlobalVector> &Vtotal,nvector<double>& Jtotal, VectorInterface& u, VectorInterface& oldu,VectorInterface& v, VectorInterface& oldv, VectorInterface& f, int ADAITER);
   
    double DWRResidual(VectorInterface& f, VectorInterface& u,
		       string U1, VectorInterface& G1, string U2, VectorInterface& G2, string U3, VectorInterface& G3,
		       string U4, VectorInterface& G4, string U5, VectorInterface& G5, string U6, VectorInterface& G6,
		       string U7, VectorInterface& G7, string U8, VectorInterface& G8, string U9, VectorInterface& G9,
		       double rhswgt, double formwgt,
		       DoubleVector& eta, 
		       GlobalVector& Z, GlobalVector& OLDZ,  GlobalVector& NEXTZ,
		       int MC);
    
    
    double Splitting_Form(VectorInterface& f, VectorInterface& u,
			 string U1, VectorInterface& G1, string U2, VectorInterface& G2, string U3, VectorInterface& G3,
			 string U4, VectorInterface& G4, string U5, VectorInterface& G5, string U6, VectorInterface& G6,
			 string U7, VectorInterface& G7, string U8, VectorInterface& G8, string U9, VectorInterface& G9,
			 double formwgt_H, double formwgt,
			 DoubleVector& eta, 
			 GlobalVector& Z, GlobalVector& OLDZ,  GlobalVector& NEXTZ,
			 int MC);
    
    void ETAProduct(DoubleVector& eta, const GlobalVector& F, const GlobalVector& Z, double wgt, int MC);
  
 
    void SolveDualProblem(vector<GlobalVector>& Ztotal,vector<GlobalVector>& Qtotal,
			  vector<GlobalVector>& Vtotal,vector<GlobalVector>& Htotal,
			  VectorInterface& f,
			  VectorInterface& z,  VectorInterface& nextz,
			  VectorInterface& q,  VectorInterface& nextq, 
			  VectorInterface& v,VectorInterface& oldv,  VectorInterface& nextv, 
			  VectorInterface& h,VectorInterface& oldh,  VectorInterface& nexth,
			  int ADAITER);

    void EstimateDualError(DoubleVector& eta,
			   DoubleVector& eta1,
			   DoubleVector& eta2,
			   DoubleVector& eta3,
			   DoubleVector& eta4,
			   DoubleVector& eta5,
			   DoubleVector& eta6,
			   DoubleVector& eta7,
			   vector<GlobalVector>& Htotal,
			   vector<GlobalVector>& Ztotal,
			   vector<GlobalVector>& Vtotal,
			   vector<GlobalVector>& Qtotal,
			   VectorInterface& h,VectorInterface& oldh,VectorInterface& nexth,
			   VectorInterface& z,VectorInterface& nextz,
			   VectorInterface& v,VectorInterface& oldv,VectorInterface& nextv,
			   VectorInterface& q,VectorInterface& nextq,
			   VectorInterface& f );
 
 
 
    double EstimateDWRprim(DoubleVector& eta, int m, vector<GlobalVector>& U,vector<GlobalVector>& V,
			 vector<GlobalVector>& Z,VectorInterface& u, VectorInterface& oldu,VectorInterface& v, VectorInterface& oldv,VectorInterface& f);
 

    double EstimateDWRdual(DoubleVector& eta2, int m,
			   vector<GlobalVector>& Htotal, vector<GlobalVector>& Vtotal, vector<GlobalVector>& Ztotal, vector<GlobalVector>& Qtotal,
			   VectorInterface& z, VectorInterface& nextz,
			   VectorInterface& v,  VectorInterface& nextv, 
			   VectorInterface& q, VectorInterface& nextq,
            VectorInterface& h, VectorInterface& oldh,
			   VectorInterface& f);
      
 
    double EstimateDWRprimBurger(DoubleVector& eta, int m, vector<GlobalVector>& V,vector<GlobalVector>& U,
				 vector<GlobalVector>& Q,VectorInterface& v, VectorInterface& oldv,
				 VectorInterface& u, VectorInterface& oldu,VectorInterface& f);
    
    
    double  EstimateDWRdualBurger(DoubleVector& eta, int m, 
			       vector<GlobalVector>&Vtotal,
			       vector<GlobalVector>&Utotal,
			       vector<GlobalVector>&Ztotal,
			       vector<GlobalVector>&Qtotal,
			       VectorInterface& v, VectorInterface& oldv,
			       VectorInterface& newv,VectorInterface& u, 
			       VectorInterface& oldu,VectorInterface& newu,
			       VectorInterface& q,VectorInterface& oldq,
			       VectorInterface& z,VectorInterface& oldz,VectorInterface& f);


    void Splittingerror(DoubleVector& eta, int m, vector<GlobalVector>& V,vector<GlobalVector>& U,
			vector<GlobalVector>& Q,VectorInterface& v, VectorInterface& oldv,VectorInterface& u, VectorInterface& oldu,VectorInterface& f) ;
    
    void SplittingErrorDual(DoubleVector& eta, int m,
			      vector<GlobalVector>& V,vector<GlobalVector>& H,
			      vector<GlobalVector>& Q,vector<GlobalVector>& Z,
			      VectorInterface& q, VectorInterface& nextq,
                  VectorInterface& z, VectorInterface& nextz,
			      VectorInterface& v, VectorInterface& nextv,
			      VectorInterface& h,
                  VectorInterface& oldh,VectorInterface& f);
    
    
    void SplittingErrorBurgerDual(DoubleVector& eta, int m,
				    vector<GlobalVector>& Vtotal,vector<GlobalVector>& Htotal,
				    vector<GlobalVector>& Qtotal,vector<GlobalVector>& Ztotal,
				    VectorInterface& q, VectorInterface& nextq,VectorInterface& z,  VectorInterface& v, 
				    VectorInterface& h,    VectorInterface& oldh,VectorInterface& f);
    
    
    
    void run(const std::string& problemlabel); 
    
  };
}


/*----------------------------   loop.h     ---------------------------*/
#endif
/*----------------------------   loop.h     ---------------------------*/
