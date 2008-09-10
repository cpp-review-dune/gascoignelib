#ifndef  __TimeSolver_h
#define  __TimeSolver_h

#include  "stdsolver.h"
#include  "simplematrix.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{

//////////////////////////////////////////////
///
///@brief
/// Default nonlinear Solver for time-dependent Equations

///
///
//////////////////////////////////////////////


class TimeSolver : public StdSolver
{
private:
  
  TimePattern       _TP;
  MatrixInterface*  _MMP;
  double            theta, dt, time;

protected:

  MatrixInterface*& GetMassMatrixPointer() {return _MMP;}

public:

  TimeSolver() : StdSolver(), _MMP(NULL), theta(0.), dt(0.), time(0.) {}
  ~TimeSolver() { if (_MMP) { delete _MMP; _MMP=NULL;} };

  string GetName() const {  return "TimeSolver";}

  void SetTimeData(double _dt, double _theta, double _time) 
  {
    dt    = _dt; 
    theta = _theta;
    time  = _time;
    assert(dt>0.);
    assert(theta>0.);
    GetProblemDescriptor()->SetTime(time,dt);
  };

  const MatrixInterface* GetMassMatrix() const {return _MMP;}
        MatrixInterface* GetMassMatrix()       {return _MMP;}

  void RegisterMatrix()
  {
    const Equation*  EQ = GetProblemDescriptor()->GetEquation();
    assert(EQ);
    int ncomp = EQ->GetNcomp();
    
    if (GetMassMatrixPointer()==NULL)
      GetMassMatrixPointer() = NewMassMatrix(ncomp,_matrixtype);
    
    StdSolver::RegisterMatrix();
  }

  void SetProblem(const ProblemDescriptorInterface& PDX)
  {
    const Equation* EQ = PDX.GetEquation();
    if (EQ) 
      {
	_TP.reservesize(EQ->GetNcomp(),EQ->GetNcomp(),0.);
	EQ->SetTimePattern(_TP);
      }   
    StdSolver::SetProblem(PDX);
  }

  void ReInitMatrix() 
  {
    GetDiscretization()->InitFilter(_PF);
    SparseStructure SA;
    GetDiscretization()->Structure(&SA);
    
    GetMatrix()->ReInit(&SA);
    GetIlu()   ->ReInit(&SA);
    
    GetMassMatrix()->ReInit(&SA);
    GetMassMatrix()->zero();
    GetDiscretization()->MassMatrix(*GetMassMatrix()); 
  }

  MatrixInterface* NewMassMatrix(int ncomp, const string& matrixtype)
    {
      return new SimpleMatrix;
    }

  void AssembleMatrix(const VectorInterface& gu, double d)
  {
    StdSolver::AssembleMatrix(gu,d);

    double scale = d/(dt*theta);
    GetMatrix()->AddMassWithDifferentStencil(GetMassMatrix(),_TP,scale);
    
    StdSolver::DirichletMatrix();
  }
  
  void Form(VectorInterface& gy, const VectorInterface& gx, double d) const
  {
    StdSolver::Form(gy,gx,d);
   
    double scale = d/(dt*theta);
    MassMatrixVector(gy,gx,scale);
  }

  void MassMatrixVector(VectorInterface& gf, const VectorInterface& gu, double d) const
  {
          GlobalVector& f = GetGV(gf);
    const GlobalVector& u = GetGV(gu);
    GetMassMatrix()->vmult_time(f,u,_TP,d);
  }
};

/*-------------------------------------------------------------*/
}

#endif
