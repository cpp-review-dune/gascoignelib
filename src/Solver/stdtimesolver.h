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


  double _dt, _theta, _time;
  fixarray<2,double> _rhs;

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

  void SetTimeData(double dt, double theta, double time, double oldrhs = -1., double newrhs = 1.);
  void SetProblem(const ProblemDescriptorInterface& PDX);

  void TimeRhsOperator(BasicGhostVector& f, const BasicGhostVector& u) const;
  void TimeRhs(int k, BasicGhostVector& f) const;
  void Form (BasicGhostVector& y, const BasicGhostVector& x, double d) const;
  void AssembleMatrix(BasicGhostVector& u, double d);
  std::string GetName() const;
  void L2Projection(BasicGhostVector& u);
};

#endif
