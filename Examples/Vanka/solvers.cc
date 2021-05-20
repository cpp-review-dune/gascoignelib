#include "solvers.h"

#include "columndiagstencil.h"
#include "fmatrixblock.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "sparseblockmatrix.h"
#include <eigen3/Eigen/Dense>

using namespace std;

namespace Gascoigne {

template<int DIM>
void
FSISolver<DIM>::smooth(int niter,
                       VectorInterface& x,
                       const VectorInterface& y,
                       VectorInterface& h) const
{
  if (GetSolverData().GetLinearSmooth() == "vanka") {
    double omega = GetSolverData().GetOmega();
    const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*>(GetMesh());
    assert(M);
    const SparseBlockMatrix<FMatrixBlock<DIM>>* A =
      dynamic_cast<const SparseBlockMatrix<FMatrixBlock<DIM>>*>(GetMatrix());
    assert(A);
    const ColumnDiagStencil* SA =
      dynamic_cast<const ColumnDiagStencil*>(A->GetStencil());
    assert(SA);

    std::vector<int> smooth_weight(GetGV(x).n(), 0);
    GlobalVector smooth_zwisch;
    smooth_zwisch.ncomp() = GetGV(x).ncomp();
    smooth_zwisch.reservesize(GetGV(x).n());
    smooth_zwisch.zero();

    for (int iter = 0; iter < niter; iter++) {
      MatrixResidual(h, x, y);
      for (int p = 0; p < M->npatches(); ++p) {
        const nvector<int>& iop = *(M->IndicesOfPatch(p));
        HASHMAP<int, int> inP;
        for (int i = 0; i < iop.size(); ++i)
          inP[iop[i]] = i;

        int N = iop.size();
        assert(N == 9);
        int ncomp = GetGV(x).ncomp();

        // copy matrix & Vector
        Eigen::MatrixXd P(N * ncomp, N * ncomp);
        Eigen::VectorXd H(N * ncomp);

        for (int r = 0; r < N; ++r) {
          int row = iop[r];
          // matrix
          for (int pos = SA->start(row); pos < SA->stop(row); ++pos) {
            int col = SA->col(pos);
            if (inP.find(col) == inP.end())
              continue;
            int c = inP[col];
            const FMatrixBlock<DIM>& B = (*A->mat(pos));
            for (int cr = 0; cr < ncomp; ++cr)
              for (int cc = 0; cc < ncomp; ++cc)
                P(ncomp * r + cr, ncomp * c + cc) = B(cr, cc);
          }

          // vector residuum
          for (int cr = 0; cr < ncomp; ++cr)
            H(ncomp * r + cr, 0) = GetGV(h)(row, cr);
        }
        // local solve
        Eigen::VectorXd X = P.inverse() * H;

        // Update
        // for (int r=0;r<N;++r)
        //	{
        //				for (int cr=0;cr<ncomp;++cr)
        //				{
        //					GetGV(x)(iop[r],cr) += omega *
        // X(r*ncomp+cr,0);
        //				}
        //	}
        // Update Jacobi
        //	for (int r=0;r<N;++r)
        //	{
        //		if(smooth_weight[iop[r]]==0)
        //				for (int cr=0;cr<ncomp;++cr)
        //				{
        //
        //						GetGV(x)(iop[r],cr) +=
        // omega
        //* X(r*ncomp+cr,0);
        // smooth_weight[iop[r]]+=1;
        //				}
        //	}

        // update mean Jacobi
        for (int r = 0; r < N; ++r) {
          for (int cr = 0; cr < ncomp; ++cr) {
            smooth_zwisch(iop[r], cr) += omega * X(r * ncomp + cr, 0);
          }
          smooth_weight[iop[r]]++;
        }
      }

      for (int i = 0; i < GetGV(x).n(); ++i)
        for (int ii = 0; ii < GetGV(x).ncomp(); ++ii) {
          GetGV(x)(i, ii) += omega / smooth_weight[i] * smooth_zwisch(i, ii);
        }
    }
  } else
    StdSolver::smooth(niter, x, y, h);
}

template class FSISolver<2>;
template class FSISolver<3>;
template class FSIMultiLevelSolver<2>;
template class FSIMultiLevelSolver<3>;
} // namespace Gascoigne
