#include "boundaryfsi.h"
#include "filescanner.h"

extern double __DT;
extern double __THETA;
extern double __TIME;

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne {

//////////////////////////////////////////////////

////////////////////////////////////////////////// BOUNDARY
template<int DIM>
BoundaryFSI<DIM>::BoundaryFSI(const ParamFile* pf)
{
  DataFormatHandler DFH;
  DFH.insert("nu_f", &__nu_f, 0.0);
  DFH.insert("rho_f", &__rho_f);
  FileScanner FS(DFH, pf, "Equation");
  p_2 = 2.266 * 1.0e4;
  p_4 = 2.286 * 1.0e4;
  cout << "%%%%%%%%%%Fluid_Stat%%%%%%%%%%" << endl;
  cout << "Boundary 2/4 -- do-nothing with p=" << p_2
       << "g/cm/s^2=" << p_2 / 1333.22 << "mmHg and p=" << p_4
       << "g/cm/s^2=" << p_4 / 1333.22 << "mmHg" << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
}

template<int DIM>
void
BoundaryFSI<DIM>::Form(VectorIterator b,
                       const FemFunction& U_Dummy,
                       const TestFunction& N,
                       int col) const
{}

template<int DIM>
void
BoundaryFSI<DIM>::Matrix(EntryMatrix& A,
                         const FemFunction& U_Dummy,
                         const TestFunction& M,
                         const TestFunction& N,
                         int col) const
{}

template<int DIM>
void
BoundaryFSI<DIM>::pointboundary(double h,
                                const FemFunction& U_Dummy,
                                const Vertex<DIM>& v,
                                const Vertex<DIM>& n) const
{

  __n[0] = n[0];
  __n[1] = n[1];
  if (DIM == 3)
    __n[2] = n[2];
}

template class BoundaryFSI<2>;
template class BoundaryFSI<3>;

} // namespace Gascoigne

/*-----------------------------------------*/
