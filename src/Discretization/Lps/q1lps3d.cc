#include  "q1lps3d.h"

namespace Gascoigne
{
/* ----------------------------------------- */

void Q1Lps3d::BasicInit(const ParamFile* paramfile)
{
  Q13d::BasicInit(paramfile);
  S    .BasicInit(paramfile,HN);
}

/* ----------------------------------------- */

void Q1Lps3d::ReInit(const MeshInterface* M)
{
  Q13d::ReInit(M);
  S    .ReInit(M);
}

/* ----------------------------------------- */

void Q1Lps3d::Structure(SparseStructureInterface* SI) const
{
  S.Structure(SI);
}

/* ----------------------------------------- */

void Q1Lps3d::StabForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  S.Form(f,u,EQ,d);
}

/* ----------------------------------------- */

void Q1Lps3d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  Q13d::Form(f,u,EQ,d);
  S    .Form(f,u,EQ,d);
}

/* ----------------------------------------- */

void Q1Lps3d::Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const
{
  Q13d::Matrix(A,u,EQ,d);
  S    .Matrix(A,u,EQ,d);
}
}
