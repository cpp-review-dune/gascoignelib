#ifndef CUSPARSEHELPER_H
#define CUSPARSEHELPER_H

#include <vector>

#include "columndiagstencil.h"
#include "fmatrixblock.h"

#include "check_cuda.h"

namespace Gascoigne {

// This is important for alignment issues.
#ifdef __MATRIX_SINGLE_PRECISION__
typedef uint32_t uint;
#else
typedef uint64_t uint;
#endif

/// Is not templated and runs the cuda part.
void
invert_host(size_t n_entries,
            size_t n_comp,
            const MatrixEntryType* src_vals,
            const IndexType* src_diag,
            MatrixEntryType* dest_vals,
            int32_t* dest_rows,
            int32_t* dest_cols);

} // namespace Gascoigne
#endif