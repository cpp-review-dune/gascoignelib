/*----------------------------   solvers.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __solvers_H
#define __solvers_H
/*----------------------------   solvers.h     ---------------------------*/

#include "stdmultilevelsolver.h"
#include "stdsolver.h"
#include "alebasediscretization.h"
#include "umfilu.h"

using namespace std;

namespace Gascoigne {
template <int DIM>
class FSISolver : public StdSolver {
private:
  DiscretizationInterface* NewDiscretization(int dimension, const std::string& discname);
  void inline DeleteSolidPressure_DIV(VectorInterface& gf, bool zero) const;
  void inline AssembleMatrix_DIV(const VectorInterface& gu, double d);

public:
  // >>>>>>>>>>>>>>>>> ILU STUFF
  std::string GetName() const {
    return "FSI Solver";
  }

  AleBaseDiscretization* GetAleDiscretization() {
    assert(dynamic_cast<AleBaseDiscretization*>(GetDiscretization()));
    return dynamic_cast<AleBaseDiscretization*>(GetDiscretization());
  }
  const AleBaseDiscretization* GetAleDiscretization() const {
    assert(dynamic_cast<const AleBaseDiscretization*>(GetDiscretization()));
    return dynamic_cast<const AleBaseDiscretization*>(GetDiscretization());
  }

  void ReInitInterface(AleBaseDiscretization* ALEDISC);
  void reinit_element(int                             en,
                      const nvector<int>&             indices,
                      HASHMAP<int, std::vector<int>>& solid_interface_cells,
                      HASHMAP<int, std::vector<int>>& fluid_interface_cells,
                      HASHSET<int>&                   fluid_cells,
                      HASHSET<int>&                   solid_cells,
                      HASHSET<int>&                   interface_nodes,
                      set<int>&                       fluid_nodes,
                      set<int>&                       solid_nodes);

  void NewMesh(const MeshInterface* mp);
  void SetBoundaryVectorZero(VectorInterface& gf) const;
  void SetBoundaryVector(VectorInterface& gf) const;
  void DeleteSolidPressure(VectorInterface& gf, bool zero) const;

  void AssembleMatrix(const VectorInterface& gu, double d);

  void PointVisu(const string& name, const GlobalVector& u, int iter) const;

  //      void ComputeIlu(const VectorInterface& gu) const;
};

template <int DIM>
class FSIMultiLevelSolver : public StdMultiLevelSolver {
public:
  std::string GetName() const {
    return "FSI MultiLevelSolver";
  }

  SolverInterface* NewSolver(int solverlevel) {
    return new FSISolver<DIM>;
  }

  const FSISolver<DIM>* GetFSISolver(int l) const {
    assert(dynamic_cast<const FSISolver<DIM>*>(GetSolver(l)));
    return dynamic_cast<const FSISolver<DIM>*>(GetSolver(l));
  }
};

}  // namespace Gascoigne

/*----------------------------   solvers.h     ---------------------------*/
/* end of #ifndef __solvers_H */
#endif
/*----------------------------   solvers.h     ---------------------------*/
