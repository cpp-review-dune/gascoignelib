#include  "laplace3d.h"
#include  "filescanner.h"

/*-----------------------------------------*/

Laplace3d::Laplace3d(const std::string& filename) : 
  Laplace(filename)
{
  DataFormatHandler DFH;
  DFH.insert("betax",&betax,0.);
  DFH.insert("betay",&betay,0.);
  DFH.insert("betaz",&betaz,0.);
  FileScanner FS(DFH,filename,"Equation");
}

/*-----------------------------------------*/

void Laplace3d::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
  b[0] += visc* (U[0].x()*N.x()+U[0].y()*N.y()+U[0].z()*N.z());
  b[0] += (betax * U[0].x() + betay * U[0].y() + betaz * U[0].z()) * N.m();
}

/*-----------------------------------------*/

void Laplace3d::Matrix(EntryMatrix& A, const FemFunction& U, const DerivativeVector& M, const TestFunction& N) const
{
  A(0,0) += visc* (M.x()*N.x()+M.y()*N.y()+M.z()*N.z());
  A(0,0) += (betax * M.x() + betay * M.y() + betaz * M.z()) * N.m();
}

/*-----------------------------------------*/
