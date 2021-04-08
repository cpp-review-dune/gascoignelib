#include "fsisparseblockilu.h"
#include "fmatrixblock.h"

namespace Gascoigne {
template <class B>
void FSISparseBlockIlu<B>::copy_entries(const HASHMAP<int, int> &G2L,
                                        const HASHSET<int> &INT,
                                        const MatrixInterface *A) {
  const SparseBlockMatrix<B> *MM =
      dynamic_cast<const SparseBlockMatrix<B> *>(A);
  assert(MM);
  const ColumnDiagStencil *AS =
      dynamic_cast<const ColumnDiagStencil *>(MM->GetStencil());
  assert(AS);

  B diagall;
  diagall.zero();
  for (int i = 0; i < diagall.ncomp(); ++i)
    diagall.uno_diag(i);
  B diagvel;
  diagvel.zero();
  for (int i = 1; i < diagvel.ncomp(); ++i)
    diagvel.uno_diag(i);

  for (int i = 0; i < SparseBlockMatrix<B>::US.n(); i++) {
    int pi = this->p[i];

    // row in domain? row on interface
    auto it_row = G2L.find(pi);
    bool row_found = (it_row != G2L.end());
    bool on_int = (INT.find(pi) != INT.end());

    for (int posA = AS->start(pi); posA < AS->stop(pi); posA++) {
      int j = AS->col(posA);

      // column in domain?
      auto it_col = G2L.find(j);
      bool col_found = (it_col != G2L.end());

      int pj = this->q[j];
      bool found = 0;
      for (int pos = SparseBlockMatrix<B>::US.start(i);
           pos < SparseBlockMatrix<B>::US.stop(i); pos++) {
        int k = SparseBlockMatrix<B>::US.col(pos);
        if (k == pj) {
          found = 1;
          if ((col_found) && (row_found) && !((on_int))) {
            *SparseBlockMatrix<B>::mat(pos) += *MM->mat(posA);
            break;
          }
          if (on_int)
            if (pi == pj) {
              *SparseBlockMatrix<B>::mat(pos) += diagvel;
              break;
            }
          if (!(row_found))
            if (pi == pj) {
              *SparseBlockMatrix<B>::mat(pos) += diagall;
              break;
            }
          break;
        }
      }
      if (!found) {
        std::cout << "not found " << std::endl;
        std::cout << *MM->mat(posA) << std::endl;

        //	      *SparseBlockMatrix<B>::mat(SparseBlockMatrix<B>::US.diag(i))
        //+= *MM.mat(posA);
      }
    }
  }
}

template class FSISparseBlockIlu<FMatrixBlock<3>>;
template class FSISparseBlockIlu<FMatrixBlock<4>>;

} // namespace Gascoigne
