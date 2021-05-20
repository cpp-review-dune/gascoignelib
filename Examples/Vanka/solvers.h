/*----------------------------   solvers.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __solvers_H
#define __solvers_H
/*----------------------------   solvers.h     ---------------------------*/

#include "stdmultilevelsolver.h"
#include "stdsolver.h"

using namespace std;

namespace Gascoigne {

template<int DIM>
class FSISolver : public StdSolver
{
private:
public:
  void smooth(int niter,
              VectorInterface& x,
              const VectorInterface& y,
              VectorInterface& h) const;
};

template<int DIM>
class FSIMultiLevelSolver : public StdMultiLevelSolver
{
public:
  std::string GetName() const { return "FSI MultiLevelSolver"; }

  SolverInterface* NewSolver(int solverlevel) { return new FSISolver<DIM>; }

  const FSISolver<DIM>* GetFSISolver(int l) const
  {
    assert(dynamic_cast<const FSISolver<DIM>*>(GetSolver(l)));
    return dynamic_cast<const FSISolver<DIM>*>(GetSolver(l));
  }
};

} // namespace Gascoigne

/*----------------------------   solvers.h     ---------------------------*/
/* end of #ifndef __solvers_H */
#endif
/*----------------------------   solvers.h     ---------------------------*/
