#include "umfilu.h"
#include <fstream>

#ifdef __WITH_UMFPACK__

using namespace std;

#define UMFPACK_OK 0
#define UMFPACK_INFO 90
#define UMFPACK_CONTROL 20

#define UMFPACK_A (0)  /* Ax=b		*/
#define UMFPACK_At (1) /* A'x=b	*/

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/

namespace Gascoigne {
extern "C" long
umfpack_dl_symbolic(long n,
                    long m,
                    const long Ap[],
                    const long Ai[],
                    const double Ax[],
                    void** Symbolic,
                    const double Control[UMFPACK_CONTROL],
                    double Info[UMFPACK_INFO]);
extern "C" long
umfpack_dl_numeric(const long Ap[],
                   const long Ai[],
                   const double Ax[],
                   void* Symbolic,
                   void** Numeric,
                   const double Control[UMFPACK_CONTROL],
                   double Info[UMFPACK_INFO]);
extern "C" long
umfpack_dl_solve(long sys,
                 const long Ap[],
                 const long Ai[],
                 const double Ax[],
                 double X[],
                 const double B[],
                 void* Numeric,
                 const double Control[UMFPACK_CONTROL],
                 double Info[UMFPACK_INFO]);
extern "C" void
umfpack_dl_free_symbolic(void** Symbolic);
extern "C" void
umfpack_dl_free_numeric(void** Numeric);

extern "C" long
umfpack_dl_triplet_to_col(long n,
                          long nz,
                          const long Ti[],
                          const long Tj[],
                          const double Tx[],
                          long Bp[],
                          long Bi[],
                          double Bx[]);

extern "C" void
umfpack_dl_report_status(const double Control[UMFPACK_CONTROL], long status);
extern "C" void
umfpack_dl_report_info(const double Control[UMFPACK_CONTROL],
                       const double Info[UMFPACK_INFO]);
extern "C" void
umfpack_dl_report_control(const double Control[UMFPACK_CONTROL]);
extern "C" long
umfpack_dl_report_symbolic(const char name[],
                           void* Symbolic,
                           const double Control[UMFPACK_CONTROL]);
extern "C" long
umfpack_dl_report_numeric(const char name[],
                          void* Numeric,
                          const double Control[UMFPACK_CONTROL]);
extern "C" void
umfpack_dl_report_control(const double Control[UMFPACK_CONTROL]);
extern "C" void
umfpack_dl_defaults(const double Control[UMFPACK_CONTROL]);

/* ----------------------------------------- */

UmfIluLong::UmfIluLong(const MatrixInterface* A)
  : SimpleMatrix()
  , Control(NULL)
  , Info(NULL)
  , Symbolic(NULL)
  , Numeric(NULL)
{
  AP = dynamic_cast<const SimpleMatrix*>(A);
  assert(AP);

  Control = new double[UMFPACK_CONTROL];
  umfpack_dl_defaults(Control);
  Control[0] = 2;
}

/* ----------------------------------------- */

UmfIluLong::~UmfIluLong()
{
  umfpack_dl_free_symbolic(&Symbolic);
  umfpack_dl_free_numeric(&Numeric);

  if (Control)
    delete[] Control;
  Control = NULL;
  if (Info)
    delete[] Info;
  Info = NULL;
}

/*-------------------------------------------------------------*/

void
UmfIluLong::ReInit(const SparseStructureInterface* SS)
{
  const ColumnStencil* SA =
    dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);
  SimpleMatrix::ReInit(SA->n(), SA->nentries());

  umfpack_dl_free_symbolic(&Symbolic);

  long n = SA->n();

  IndexVector start = SA->start();
  nvector<long> start_long;

  start_long.resize(start.size());
  for (size_t i = 0; i < start.size(); i++) {
    start_long[i] = static_cast<long>(start[i]);
  }

  IndexVector col = SA->col();
  nvector<long> col_long;

  col_long.resize(col.size());
  for (size_t i = 0; i < col.size(); i++) {
    col_long[i] = static_cast<long>(col[i]);
  }

  const long* sb = &(*start_long.begin());
  const long* cb = &(*col_long.begin());

  long status =
    umfpack_dl_symbolic(n, n, sb, cb, NULL, &Symbolic, Control, Info);

  if (status != UMFPACK_OK) {
    umfpack_dl_report_info(Control, Info);
    umfpack_dl_report_status(Control, status);
    ofstream file("MATRIX");
    //       AP->Write(file);
    cerr << "umfpack_symbolic failed\n";
    exit(1);
  }
}

/*-------------------------------------------------*/

void
UmfIluLong::copy_entries(const MatrixInterface& A)
{}

/*-------------------------------------------------------------*/

void
UmfIluLong::ConstructStructure(const IndexVector& perm,
                               const MatrixInterface& A)
{}

/*-----------------------------------------*/

void
UmfIluLong::Factorize()
{
  //
  // baue LU auf
  //

  umfpack_dl_free_numeric(&Numeric);

  const ColumnStencil* SA =
    dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);

  IndexVector start = SA->start();
  nvector<long> start_long;

  start_long.resize(start.size());
  for (size_t i = 0; i < start.size(); i++) {
    start_long[i] = static_cast<long>(start[i]);
  }

  IndexVector col = SA->col();
  nvector<long> col_long;

  col_long.resize(col.size());
  for (size_t i = 0; i < col.size(); i++) {
    col_long[i] = static_cast<long>(col[i]);
  }

  const long* sb = &(*start_long.begin());
  const long* cb = &(*col_long.begin());
  const double* mb = &AP->GetValue(0);
  long status =
    umfpack_dl_numeric(sb, cb, mb, Symbolic, &Numeric, Control, Info);
  if (status != UMFPACK_OK) {
    umfpack_dl_report_info(Control, Info);
    umfpack_dl_report_status(Control, status);
    cerr << "umfpack_numeric failed\n";
    exit(1);
  }
  //   umfpack_report_numeric("LU von A\n",Numeric,Control);
}

/*-----------------------------------------*/

void
UmfIluLong::Solve(DoubleVector& x, const DoubleVector& b) const
{
  const ColumnStencil* SA =
    dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);

  IndexVector start = SA->start();
  nvector<long> start_long;

  start_long.resize(start.size());
  for (size_t i = 0; i < start.size(); i++) {
    start_long[i] = static_cast<long>(start[i]);
  }

  IndexVector col = SA->col();
  nvector<long> col_long;

  col_long.resize(col.size());
  for (size_t i = 0; i < col.size(); i++) {
    col_long[i] = static_cast<long>(col[i]);
  }

  const long* sb = &(*start_long.begin());
  const long* cb = &(*col_long.begin());
  const double* mb = &AP->GetValue(0);
  double* xb = &(*x.begin());
  const double* bb = &(*b.begin());
  long status =
    umfpack_dl_solve(UMFPACK_At, sb, cb, mb, xb, bb, Numeric, Control, Info);

  if (status != UMFPACK_OK) {
    umfpack_dl_report_info(Control, Info);
    umfpack_dl_report_status(Control, status);
    cerr << "umfpack_dl_solve failed\n";
    exit(1);
  }
}

/*-----------------------------------------*/

void
UmfIluLong::SolveTranspose(DoubleVector& x, const DoubleVector& b)
{
  const ColumnStencil* SA =
    dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);

  IndexVector start = SA->start();
  nvector<long> start_long;

  start_long.resize(start.size());
  for (size_t i = 0; i < start.size(); i++) {
    start_long[i] = static_cast<long>(start[i]);
  }

  IndexVector col = SA->col();
  nvector<long> col_long;

  col_long.resize(col.size());
  for (size_t i = 0; i < col.size(); i++) {
    col_long[i] = static_cast<long>(col[i]);
  }

  const long* sb = &(*start_long.begin());
  const long* cb = &(*col_long.begin());
  const double* mb = &AP->GetValue(0);
  double* xb = &(*x.begin());
  const double* bb = &(*b.begin());
  long status =
    umfpack_dl_solve(UMFPACK_A, sb, cb, mb, xb, bb, Numeric, Control, Info);

  if (status != UMFPACK_OK) {
    umfpack_dl_report_info(Control, Info);
    umfpack_dl_report_status(Control, status);
    cerr << "umfpack_dl_solve failed\n";
    exit(1);
  }
}
} // namespace Gascoigne

#undef UMFPACK_OK
#undef UMFPACK_INFO
#undef UMFPACK_CONTROL

#undef UMFPACK_A
#undef UMFPACK_At

#endif
