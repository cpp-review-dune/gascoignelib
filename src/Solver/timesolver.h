#ifndef  __TimeSolver_h
#define  __TimeSolver_h

#include  "stdsolver.h"
#include  "simplematrix.h"
#include  "cginfo.h"

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

  TimeSolver() : StdSolver(), _MMP(NULL), theta(1.), dt(0.), time(0.) {}
  ~TimeSolver() { if (_MMP) { delete _MMP; _MMP=NULL;} };

  std::string GetName() const {  return "TimeSolver";}

  void SetTimeData(double _dt, double _theta, double _time);

  const MatrixInterface* GetMassMatrix() const {return _MMP;}
        MatrixInterface* GetMassMatrix()       {return _MMP;}

  void RegisterMatrix();
  void SetProblem(const ProblemDescriptorInterface& PDX);

  void ReInitMatrix();

  MatrixInterface* NewMassMatrix(int ncomp, const std::string& matrixtype)
    {
      return new SimpleMatrix;
    }

  void AssembleMatrix(const VectorInterface& gu, double d);
  void Form(VectorInterface& gy, const VectorInterface& gx, double d) const;
  void MassMatrixVector(VectorInterface& gf, const VectorInterface& gu, double d) const;
  void InverseMassMatrix(VectorInterface& u, const VectorInterface& f, CGInfo& info);
  void precondition(VectorInterface& u, const VectorInterface& f);
  void cgvmult(VectorInterface& y, const VectorInterface& x, double d) const;
};

/*-------------------------------------------------------------*/
}

#endif
