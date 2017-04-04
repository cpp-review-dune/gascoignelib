/*----------------------------   dwr_discretization.h     ---------------------------*/
/*      $Id: dwr_discretization.h,v 1.8 2009/09/15 16:52:52 richter Exp $                 */
#ifndef __dwr_discretization_H
#define __dwr_discretization_H
/*----------------------------   dwr_discretization.h     ---------------------------*/


#include "q22d.h"
#include "q23d.h"
#include "baseq1patch.h"
#include "baseq22d.h"
#include "baseq13dpatch.h"
#include "transformation2d.h"
#include "transformation3d.h"
#include "finiteelement.h"
#include "dwr_integrator.h"
#include "hnstructureq12d.h"
#include "hnstructureq22d.h"


/*---------------------------------------------------*/

namespace Gascoigne
{
  class DWR_Discretization : public PatchDiscretization
  {
  protected:
    typedef Transformation2d<BaseQ12d>         TransQ1;
    typedef Transformation2d<BaseQ22d>         TransQ2;
    typedef Transformation2d<BaseQ12dPatch>    TransQ1Patch;
        
    HNStructureInterface *HNTest, *HNAnsatz;


    FemInterface* __FEM_ANSATZ;
    FemInterface* __FEM_TRIAL;

    std::string GetName() const {return "DWR Boundary Disc"; }
    
    IntVector GetLocalIndices(int iq) const { abort(); }


    FemInterface*       GetAnsatzFem()       { return __FEM_ANSATZ; }
    const FemInterface* GetAnsatzFem() const { return __FEM_ANSATZ; }
    FemInterface*       GetTrialFem()        { return __FEM_TRIAL; }
    const FemInterface* GetTrialFem() const  { return __FEM_TRIAL; }
    nvector<int>        GetAnsatzIndices(int p) const;
    nvector<int>        GetTrialIndices(int p) const;
    
    void TransformationQ1(FemInterface::Matrix& T, int ip) const;
    void TransformationQ2(FemInterface::Matrix& T, int ip) const;
    void InitElement(int ip, const FemInterface* F) const;
    
    const DWR_Integrator<2>* GetDWRIntegrator() const
    {
      assert(dynamic_cast<const DWR_Integrator<2>*> (GetIntegrator()));
      return dynamic_cast<const DWR_Integrator<2>*> (GetIntegrator());
    }

    
    void GlobalToLocalAnsatz(LocalVector& U, const GlobalVector& u, int iq) const;
    void GlobalToLocalTrial (LocalVector& U, const GlobalVector& u, int iq) const;
    void GlobalToLocalDataAnsatz(int iq) const;
    void GlobalToLocalDataTrial (int iq) const;
    void GlobalToLocalSingleAnsatz(LocalVector& U, const GlobalVector& u, int iq) const;
    void GlobalToLocalSingleTrial (LocalVector& U, const GlobalVector& u, int iq) const;

    void LocalToGlobalAnsatz(GlobalVector& f, const LocalVector& F, int iq, double s) const;
    void LocalToGlobalTrial (GlobalVector& f, const LocalVector& F, int iq, double s) const;
    
  public:
    DWR_Discretization(std::string trans_ansatz, std::string base_ansatz,
		       std::string trans_trial,  std::string base_trial);
    ~DWR_Discretization();
    
    void BasicInit(const ParamFile* paramfile);
    void ReInit(const MeshInterface* MP);

    int n()  const {return GetMesh()->nnodes();}
    int nc() const {return GetMesh()->ncells();}
    

    // 
    void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
    void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const 
    {abort();}
    void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;

    void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
    void BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const
    {abort();}
    
    void MassMatrix(MatrixInterface& M) const {abort();}
    void MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const{abort();}
    

    // Hanging nodes
    void HNAverageAnsatz   (GlobalVector& x) const;
    void HNZeroAnsatz      (GlobalVector& x) const;
    void HNDistributeTest  (GlobalVector& x) const;
    void HNAverage   (GlobalVector& x) const {}
    void HNDistribute(GlobalVector& x) const {}
    void HNZero      (GlobalVector& x) const {}
    bool HNZeroCheck (const GlobalVector& x) const
    {
      return false;
    }
    void HNAverageData() const  {}
    void HNZeroData   () const  {}


  };
}




/*----------------------------   dwr_discretization.h     ---------------------------*/
/* end of #ifndef __dwr_discretization_H */
#endif
/*----------------------------   dwr_discretization.h     ---------------------------*/

  
