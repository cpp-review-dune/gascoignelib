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

class StdTimeSolver : public StdSolver
{
private:

  TimePattern       _TP;
  MatrixInterface*  _MMP;

protected:


  mutable BasicGhostVector Gg, Gr, Gd;
  double dt, theta, time;

  const MatrixInterface* GetMassMatrix() const {return _MMP;}
  MatrixInterface* GetMassMatrix() {return _MMP;}
  MatrixInterface*& GetMassMatrixPointer() {return _MMP;}

  const TimePattern& GetTimePattern() const {return _TP;}
  TimePattern& GetTimePattern() {return _TP;}

  virtual MatrixInterface* NewMassMatrix(int ncomp, const string& matrixtype);
  void MemoryMatrix();
  void _BasicInit();

public:
  
  StdTimeSolver();
  ~StdTimeSolver();

  void BasicInit(int level,const std::string& paramfile, const MeshInterface* MP, const ProblemDescriptorInterface& PDX);

  void SetTimeData(double d, double th, double ti);

  void TimeRhs(BasicGhostVector& f, const BasicGhostVector& u) const;
  void Residual (BasicGhostVector& y, const BasicGhostVector& x, double d) const;
  void AssembleMatrix(BasicGhostVector& u, double d);
  std::string GetName() const;
  void RhsL2Projection(BasicGhostVector& f, const BasicGhostVector& u) const;
  void L2Projection(BasicGhostVector& u, const BasicGhostVector& f) const;
};

#endif
