#ifndef  __StdTimeSolver_h
#define  __StdTimeSolver_h

#include  "stdsolver.h"

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

class StdTimeSolver : public virtual StdSolver
{
private:
  
  TimePattern       _TP;
  MatrixInterface*  _MMP;

protected:


  double _dt, _theta, _time;
  fixarray<2,double> _rhs;

  MatrixInterface*& GetMassMatrixPointer() {return _MMP;}

  const TimePattern& GetTimePattern() const {return _TP;}
  TimePattern& GetTimePattern() {return _TP;}

  virtual MatrixInterface* NewMassMatrix(int ncomp, const std::string& matrixtype);
  virtual void InitialCondition(GlobalVector& f, double d=1.) const;
  virtual std::string PrecondCGMass(GlobalVector& u, GlobalVector& f, const TimePattern& TP, double s=1.);

public:
  
  StdTimeSolver();
  ~StdTimeSolver();

  void BasicInit(int level, const ParamFile* pfparamfile, const MeshInterface* MP);

  void RegisterMatrix();
  void ReInitMatrix();

  void SetTimeData(double dt, double theta, double time, double oldrhs = -1., double newrhs = 1.);
  void SetProblem(const ProblemDescriptorInterface& PDX);

  void InitialCondition(BasicGhostVector& f, double d=1.) const;
  virtual void TimeRhsOperator(BasicGhostVector& f, const BasicGhostVector& u) const;
  void TimeRhs(int k, BasicGhostVector& f) const;
  virtual void Form (BasicGhostVector& y, const BasicGhostVector& x, double d) const;
  virtual void AssembleMatrix(const BasicGhostVector& u, double d);
  std::string GetName() const;
  void L2Projection(BasicGhostVector& u);
  
  void SetMassMatrix(MatrixInterface &MM, bool init=false);
  const MatrixInterface* GetMassMatrix() const {return _MMP;}
  MatrixInterface* GetMassMatrix() {return _MMP;}
};
}

#endif
