#include  "umfilu.h"
#include  <fstream>


using namespace std;
 
#define UMFPACK_OK       0
#define UMFPACK_INFO    90
#define UMFPACK_CONTROL 20

#define UMFPACK_A	(0)	/* Ax=b		*/
#define UMFPACK_At	(1)	/* A'x=b	*/

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/

extern "C" int umfpack_di_symbolic
(
    int n,
    int m,
    const int Ap [ ],
    const int Ai [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;
extern "C" int umfpack_di_numeric
(
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    void *Symbolic,
    void **Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;
extern "C" int umfpack_di_solve
(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    double X [ ],
    const double B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;
extern "C" void umfpack_di_free_symbolic
(
    void **Symbolic
) ;
extern "C" void umfpack_di_free_numeric
(
    void **Numeric
) ;

extern "C" int umfpack_triplet_to_col
(
    int n,
    int nz,
    const int Ti [ ],
    const int Tj [ ],
    const double Tx [ ],
    int Bp [ ],
    int Bi [ ],
    double Bx [ ]
) ;

extern "C" void umfpack_di_report_status
(
    const double Control [UMFPACK_CONTROL],
    int status
) ;
extern "C" void umfpack_di_report_info
(
    const double Control [UMFPACK_CONTROL],
    const double Info [UMFPACK_INFO]
) ;
extern "C" void umfpack_di_report_control
(
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" int umfpack_di_report_symbolic
(
    const char name [ ],
    void *Symbolic,
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" int umfpack_di_report_numeric
(
    const char name [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" void umfpack_di_report_control
(
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" void umfpack_di_defaults
(
    const double Control [UMFPACK_CONTROL]
) ;

/* ----------------------------------------- */

UmfIlu::UmfIlu(const MatrixInterface* A) 
  : SimpleMatrix(), Symbolic(NULL), Numeric(NULL), Control(NULL), Info(NULL)
{
  AP = dynamic_cast<const SimpleMatrix*>(A);
  assert(AP);

  Control = new double[UMFPACK_CONTROL];
  umfpack_di_defaults(Control);
  Control[0] = 2;
}

/* ----------------------------------------- */

UmfIlu::~UmfIlu()
{
  if(Control) delete Control; Control=NULL;
  if(Info) delete Info; Info=NULL;
}

/*-------------------------------------------------------------*/

void UmfIlu::ReInit(const SparseStructureInterface* SS)
{
  const ColumnStencil* SA = dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);
  SimpleMatrix::ReInit(SA->n(),SA->nentries());

//   const UmfStructure& SA = *GetStructure();
  
  umfpack_di_free_symbolic (&Symbolic) ;

  int n = SA->n();
  const int* sb = &(*SA->start().begin());
  const int* cb = &(*SA->col().begin());
  
  int status = umfpack_di_symbolic(n, n, sb, cb, &Symbolic, Control, Info);
  status = umfpack_di_symbolic(n, n, sb, cb, &Symbolic, Control, Info);
  
  if(status != UMFPACK_OK)
    {
      umfpack_di_report_info(Control,Info);
      umfpack_di_report_status(Control,status);
      ofstream file("MATRIX");
//       AP->Write(file);
      cerr << "umfpack_symbolic failed\n"; exit(1);
    }
}

/*-------------------------------------------------*/

void UmfIlu::copy_entries(const MatrixInterface&  A)
{
}

/*-------------------------------------------------------------*/

void UmfIlu::ConstructStructure(const nvector<int>& perm, const MatrixInterface& A)
{
}

/*-----------------------------------------*/

void UmfIlu::Factorize()
{
  //
  // baue LU auf
  //

  umfpack_di_free_numeric (&Numeric) ;

  const ColumnStencil* SA = dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);

  const int* sb = &(*SA->start().begin());
  const int* cb = &(*SA->col().begin());
  const double* mb = &AP->GetValue(0);
  int status = umfpack_di_numeric(sb, cb, mb, Symbolic, &Numeric, Control, Info) ;
  if(status != UMFPACK_OK)
    {
      umfpack_di_report_info(Control,Info);
      umfpack_di_report_status(Control,status);
      cerr << "umfpack_numeric failed\n"; exit(1);
    }
  //   umfpack_report_numeric("LU von A\n",Numeric,Control);
}

/*-----------------------------------------*/

void UmfIlu::Solve(nvector<double>& x, const nvector<double>& b)
{
  const ColumnStencil* SA = dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);

  const int* sb = &(*SA->start().begin());
  const int* cb = &(*SA->col().begin());
  const double* mb = &AP->GetValue(0);
  double* xb = &(*x.begin());
  const double* bb = &(*b.begin());
  int status = umfpack_di_solve (UMFPACK_At, sb, cb, mb, xb, bb, Numeric, Control, Info) ;

    if(status != UMFPACK_OK)
      {
	umfpack_di_report_info(Control,Info);
	umfpack_di_report_status(Control,status);
	cerr << "umfpack_di_solve failed\n"; exit(1);
      }
}

/*-----------------------------------------*/

void UmfIlu::SolveTranspose(nvector<double>& x, const nvector<double>& b)
{
  const ColumnStencil* SA = dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);

  const int* sb = &(*SA->start().begin());
  const int* cb = &(*SA->col().begin());
  const double* mb = &AP->GetValue(0);
  double* xb = &(*x.begin());
  const double* bb = &(*b.begin());
  int status = umfpack_di_solve (UMFPACK_A, sb, cb, mb, xb, bb, Numeric, Control, Info) ;

    if(status != UMFPACK_OK)
      {
	umfpack_di_report_info(Control,Info);
	umfpack_di_report_status(Control,status);
	cerr << "umfpack_di_solve failed\n"; exit(1);
      }
}

 
#undef UMFPACK_OK     
#undef UMFPACK_INFO   
#undef UMFPACK_CONTROL

#undef UMFPACK_A	
#undef UMFPACK_At	
