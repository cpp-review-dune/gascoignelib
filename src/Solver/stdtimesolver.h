#ifndef  __StdTimeSolver_h
#define  __StdTimeSolver_h

#include  "stdsolver.h"

//////////////////////////////////////////////
///
///@brief
/// Default nonlinear Solver for time-dependent Equations

///
///
//////////////////////////////////////////////


/*-------------------------------------------------------------*/

class StdTimeSolver : public virtual StdSolver
{
private:

  TimePattern       _TP;
  MatrixInterface*  _MMP;

protected:


  double dt, theta, time, rhs;

  const MatrixInterface* GetMassMatrix() const {return _MMP;}
  MatrixInterface* GetMassMatrix() {return _MMP;}
  MatrixInterface*& GetMassMatrixPointer() {return _MMP;}

  const TimePattern& GetTimePattern() const {return _TP;}
  TimePattern& GetTimePattern() {return _TP;}

  virtual MatrixInterface* NewMassMatrix(int ncomp, const std::string& matrixtype);

public:
  
  StdTimeSolver();
  ~StdTimeSolver();

  void BasicInit(int level, const Gascoigne::ParamFile* pfparamfile, const MeshInterface* MP);

  void RegisterMatrix();
  void ReInitMatrix();

  void SetTimeData(double d, double th, double ti, double rh = -1);
  void SetProblem(const ProblemDescriptorInterface& PDX);

  void TimeRhs(BasicGhostVector& f, const BasicGhostVector& u) const;
  void Residual (BasicGhostVector& y, const BasicGhostVector& x, double d) const;
  void AssembleMatrix(BasicGhostVector& u, double d);
  std::string GetName() const;
  void RhsL2Projection(BasicGhostVector& f) const;
  void L2Projection(BasicGhostVector& u, const BasicGhostVector& f);
};

#endif
