/*----------------------------   multilevelsolvers.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __multilevelsolvers_H
#define __multilevelsolvers_H
/*----------------------------   multilevelsolvers.h     ---------------------------*/


#include "stdmultilevelsolver.h"
#include "solvers.h"
#include "vankasolver.h"


using namespace std;

namespace Gascoigne
{

  template<int DIM>
  class FSIMultiLevelSolver : public StdMultiLevelSolver
  {
  private:
  map<std::string,std::vector<SolverInterface*>>  _MSP;
  mutable std::string _SolverLabel;
  mutable VectorInterface _f,_u;
  
  void NewtonMatrixControl_nomonitor(VectorInterface& u, NLInfo& nlinfo);
  public:
  FSIMultiLevelSolver();     
  ~FSIMultiLevelSolver();        
    std::string GetName() const {return "FSI MultiLevelSolver";}
	  
    SolverInterface* NewSolver(int solverlevel)
    { //return new FSISolver<DIM>; 
    return new FSIVankaSolver<DIM>; 
    
    }
    
    void NewSolvers(std::string solverlabel);  
      
    void ReInit(const std::string& problemlabel);
    void NewMgInterpolator();
    void MgInterpolatorToSolver();
    void SetSolverLabel (std::string solverlabel) {_SolverLabel=solverlabel;}
    
    std::vector<SolverInterface*>& GetSolverPointers() 
    {
		if(_MSP.find(_SolverLabel)==_MSP.end())
		{ 
		 cout<<"SolverLabel: " <<_SolverLabel<< "does not exist"<<endl;
		 abort();
		}
		return _MSP.find(_SolverLabel)->second; 
    }
    const std::vector<SolverInterface*>& GetSolverPointers() const 
    { 
        if(_MSP.find(_SolverLabel)==_MSP.end())
    	{ 
    	 cout<<"SolverLabel: " <<_SolverLabel<< "does not exist"<<endl;
    	 abort();
    	}
    	return _MSP.find(_SolverLabel)->second; 
    }
    virtual SolverInterface*& GetSolverPointer(int l) 
    {
        if(_MSP.find(_SolverLabel)==_MSP.end())
    	{ 
    	 cout<<"SolverLabel: " <<_SolverLabel<< "does not exist"<<endl;
    	 abort();
    	}
		assert(l<_MSP.find(_SolverLabel)->second.size()); 
		return _MSP.find(_SolverLabel)->second[l];
    }
    
    SolverInterface* GetSolver(int l) 
    {
        if(_MSP.find(_SolverLabel)==_MSP.end())
    	{ 
    	 cout<<"SolverLabel: " <<_SolverLabel<< "does not exist"<<endl;
    	 abort();
    	}
    	assert(l<_MSP.find(_SolverLabel)->second.size()); 
    	return _MSP.find(_SolverLabel)->second[l];
    }
  	const SolverInterface* GetSolver(int l) const 
  	{
  	    if(_MSP.find(_SolverLabel)==_MSP.end())
    	{ 
    	 cout<<"SolverLabel: " <<_SolverLabel<< "does not exist"<<endl;
    	 abort();
    	}
  		assert(l<_MSP.find(_SolverLabel)->second.size());
  	 	return _MSP.find(_SolverLabel)->second[l];
  	}
  	SolverInterface* GetSolver() 
  	{
  	   	if(_MSP.find(_SolverLabel)==_MSP.end())
    	{ 
    	 cout<<"SolverLabel: " <<_SolverLabel<< "does not exist"<<endl;
    	 abort();
    	}
    	assert(_MSP.find(_SolverLabel)->second.size()==nlevels()); 
    	return _MSP.find(_SolverLabel)->second[FinestLevel()];
    }
  	const SolverInterface* GetSolver() const 
  	{
  	    if(_MSP.find(_SolverLabel)==_MSP.end())
    	{ 
    	 cout<<"SolverLabel: " <<_SolverLabel<< "does not exist"<<endl;
    	 abort();
    	}
  		assert(_MSP.find(_SolverLabel)->second.size()==nlevels()); 
  		return _MSP.find(_SolverLabel)->second[FinestLevel()];
  	}

  	SolverInterface* GetSolver(std::string solverlabel) 
  	{
  	   	if(_MSP.find(solverlabel)==_MSP.end())
    	{ 
    	 cout<<"SolverLabel: " <<solverlabel<< "does not exist"<<endl;
    	 abort();
    	}
    	assert(_MSP.find(solverlabel)->second.size()==nlevels()); 
    	return _MSP.find(solverlabel)->second[FinestLevel()];
    }
  	const SolverInterface* GetSolver(std::string solverlabel) const 
  	{
  	    if(_MSP.find(solverlabel)==_MSP.end())
    	{ 
    	 cout<<"SolverLabel: " <<solverlabel<< "does not exist"<<endl;
    	 abort();
    	}
  		assert(_MSP.find(solverlabel)->second.size()==nlevels()); 
  		return _MSP.find(solverlabel)->second[FinestLevel()];
  	}  
  
    const FSISolver<DIM>* GetFSISolver(int l) const
    {
      assert(dynamic_cast<const FSISolver<DIM>* > (GetSolver(l)));
      return dynamic_cast<const FSISolver<DIM>* > (GetSolver(l));
    }
    FSISolver<DIM>* GetFSISolver(int l) 
    {
      assert(dynamic_cast<FSISolver<DIM>* > (GetSolver(l)));
      return dynamic_cast<FSISolver<DIM>* > (GetSolver(l));
    }

    const FSISolver<DIM>* GetFSISolver() const
    {
      assert(dynamic_cast<const FSISolver<DIM>* > (GetSolver()));
      return dynamic_cast<const FSISolver<DIM>* > (GetSolver());
    }
    FSISolver<DIM>* GetFSISolver() 
    {
      assert(dynamic_cast<FSISolver<DIM>* > (GetSolver()));
      return dynamic_cast<FSISolver<DIM>* > (GetSolver());
    }

      
    double NewtonUpdate(GlobalVector& U_Vec_GV, double& rr, VectorInterface& x, VectorInterface& dx, VectorInterface& r, const VectorInterface& f, NLInfo& nlinfo,NLInfo& info_solid_disp,NLInfo& info_meshmotion);
    void AddNodeVectorinAllSolvers(const string& name, VectorInterface& gq)  ;
	void DeleteNodeVectorinAllSolvers(const string& name);
	  void RegisterVectors();
	  void NewtonUpdateDisplacement(GlobalVector& U_Vec_GV,NLInfo& info_solid_disp,NLInfo& info_meshmotion);
	  void newton(VectorInterface& U_Vec, const VectorInterface& f, VectorInterface& r, VectorInterface& w, NLInfo& info_fsi_reduced ,NLInfo& info_solid_disp,NLInfo& info_meshmotion);
	  void UpdateVelocityandPressureinU(GlobalVector& U_Vec_GV, double factor, VectorInterface& dx);
	  std::string Solve( VectorInterface& x, const VectorInterface& b, NLInfo& info_fsi_reduced ,NLInfo& info_solid_disp,NLInfo& info_meshmotion);
	  void VelocityandPressureUto_u(GlobalVector& U_Vec_GV, VectorInterface& _u);
	  void AssembleMatrix(VectorInterface& u);
	  void AssembleMatrix(VectorInterface& u, NLInfo& nlinfo);
	  void ComputeIlu(VectorInterface& u);
  };

  
}

/*----------------------------   multilevelsolvers.h     ---------------------------*/
/* end of #ifndef __multilevelsolvers_H */
#endif
/*----------------------------   multilevelsolvers.h     ---------------------------*/