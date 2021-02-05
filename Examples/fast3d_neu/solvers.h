/*----------------------------   solvers.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __solvers_H
#define __solvers_H
/*----------------------------   solvers.h     ---------------------------*/

#include "alebasediscretization.h"
#include "sparse_umf.h"
#include "stdmultilevelsolver.h"
#include "stdsolver.h"
//#include "umfilu.h"
//#include "mumpsilu.h"
#include "fmatrixblock.h"
#include "lpsequation.h"

using namespace std;

namespace Gascoigne {

template <int DIM> class FSISolver : public StdSolver {
private:
  DiscretizationInterface *NewDiscretization(int dimension,
                                             const std::string &discname);

  IluInterface *__IF, *__IS;

  SparseUmf<FMatrixBlock<DIM + 1>> _LAP_M;
  SparseBlockMatrix<FMatrixBlock<DIM + 1>> _LAP_A;

  // MumpsIlu<DIM+1>                         _LAP_M;

public:
  vector<vector<double>> _precond;

  void do_precondition(VectorInterface gx) const;
  void undo_precondition(VectorInterface gx) const;

  FSISolver();

  // >>>>>>>>>>>>>>>>> ILU STUFF
  void ReInitMatrix();
  void ReInitExtensionMatrix();

  void SolveExtension(VectorInterface &x);

  void ComputeIlu(const VectorInterface &gu) const;
  void modify_ilu(IluInterface &I, int ncomp) const;

  void smooth(int niter, VectorInterface &x, const VectorInterface &y,
              VectorInterface &h) const;
  void smooth_exact(VectorInterface &x, const VectorInterface &y,
                    VectorInterface &help) const;

  std::string GetName() const { return "FSI Solver"; }

  AleBaseDiscretization *GetAleDiscretization() {
    assert(dynamic_cast<AleBaseDiscretization *>(GetDiscretization()));
    return dynamic_cast<AleBaseDiscretization *>(GetDiscretization());
  }
  const AleBaseDiscretization *GetAleDiscretization() const {
    assert(dynamic_cast<const AleBaseDiscretization *>(GetDiscretization()));
    return dynamic_cast<const AleBaseDiscretization *>(GetDiscretization());
  }

  void ReInitInterface(AleBaseDiscretization *ALEDISC);
  void reinit_element(int en, const nvector<int> &indices,
                      HASHSET<int> &fluid_cells, HASHSET<int> &solid_cells,
                      set<int> &fluid_nodes, set<int> &solid_nodes,
                      int material);

  void reinit_interface_element(
      int en, const nvector<int> &indices,
      HASHMAP<int, std::vector<int>> &solid_interface_cells,
      HASHMAP<int, std::vector<int>> &fluid_interface_cells,
      HASHSET<int> &interface_nodes, int material);

  void Form(VectorInterface &y, const VectorInterface &x, double d) const;

  void NewMesh(const MeshInterface *mp);
  void SetBoundaryVectorZero(VectorInterface &gf) const;
  void SetBoundaryVector(VectorInterface &gf) const;
  void DeleteSolidPressure(VectorInterface &gf) const;

  void AssembleMatrix(const VectorInterface &gu, double d);

  void PointVisu(const string &name, const GlobalVector &u, int iter) const;

  //      void ComputeIlu(const VectorInterface& gu) const;

  MatrixInterface *NewMatrix(int ncomp, const std::string &matrixtype);
  IluInterface *NewIlu(int ncomp, const std::string &matrixtype);
};

template <int DIM> class FSIMultiLevelSolver : public StdMultiLevelSolver {
public:
  std::string GetName() const { return "FSI MultiLevelSolver"; }

  SolverInterface *NewSolver(int solverlevel) { return new FSISolver<DIM>; }

  const FSISolver<DIM> *GetFSISolver(int l) const {
    assert(dynamic_cast<const FSISolver<DIM> *>(GetSolver(l)));
    return dynamic_cast<const FSISolver<DIM> *>(GetSolver(l));
  }
  FSISolver<DIM> *GetFSISolver(int l) {
    assert(dynamic_cast<FSISolver<DIM> *>(GetSolver(l)));
    return dynamic_cast<FSISolver<DIM> *>(GetSolver(l));
  }

  const FSISolver<DIM> *GetFSISolver() const {
    assert(dynamic_cast<const FSISolver<DIM> *>(GetSolver()));
    return dynamic_cast<const FSISolver<DIM> *>(GetSolver());
  }
  FSISolver<DIM> *GetFSISolver() {
    assert(dynamic_cast<FSISolver<DIM> *>(GetSolver()));
    return dynamic_cast<FSISolver<DIM> *>(GetSolver());
  }

  void UpdateDeformation(VectorInterface &x);

  double NewtonUpdate(double &rr, VectorInterface &x, VectorInterface &dx,
                      VectorInterface &r, const VectorInterface &f,
                      NLInfo &nlinfo);
};

/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Laplace-matrix-stuff
template <int DIM> // laplace eq for extension
class ExtensionEquation : public LpsEquation {
  mutable double ext;
  mutable int domain;
  Chi chi;

public:
  std::string GetName() const { return "Extension"; }
  int GetNcomp() const { return DIM + 1; }

  void point(double h, const FemFunction &U, const Vertex<DIM> &v) const {
    double dx = std::max(0.0, fabs(v.x() - 0.4) - 0.4);
    double dy = std::max(0.0, fabs(v.y() - 0.2) - 0.01);
    double dist = sqrt(dx * dx + dy * dy);
    ext = 1.0 / (1.e-2 + dist);
    // domain = chi(v);
  }

  void point_cell(int material) const {
    if (material == 1)
      domain = 1;
    if (material == 2)
      domain = -1;
  }

  void Form(VectorIterator b, const FemFunction &U,
            const TestFunction &N) const {
    if (domain > 0)
      return;

    for (int i = 0; i < DIM; ++i)
      for (int j = 0; j < DIM; ++j)
        b[i] += ext * U[i][j + 1] * N[j + 1];
  }
  void Matrix(EntryMatrix &A, const FemFunction &U, const TestFunction &M,
              const TestFunction &N) const {
    if (domain > 0)
      return;

    for (int i = 0; i < DIM; ++i)
      for (int j = 0; j < DIM; ++j)
        A(i, i) += ext * M[j + 1] * N[j + 1];
  }

  void StabForm(VectorIterator b, const FemFunction &U, const FemFunction &UP,
                const TestFunction &N) const {}
  void StabMatrix(EntryMatrix &A, const FemFunction &U, const TestFunction &Np,
                  const TestFunction &Mp) const {}
  void lpspoint(double h, const FemFunction &U, const Vertex<DIM> &v) const {}
};

} // namespace Gascoigne

/*----------------------------   solvers.h     ---------------------------*/
/* end of #ifndef __solvers_H */
#endif
/*----------------------------   solvers.h     ---------------------------*/
