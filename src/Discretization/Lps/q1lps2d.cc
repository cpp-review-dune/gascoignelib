#include  "q1lps2d.h"

namespace Gascoigne
{
/* ----------------------------------------- */

void Q1Lps2d::BasicInit(const ParamFile* paramfile)
{
  Q12d::BasicInit(paramfile);
  S    .BasicInit(paramfile,HN);
}

/* ----------------------------------------- */

void Q1Lps2d::ReInit(const MeshInterface* M)
{
  Q12d::ReInit(M);
  S    .ReInit(M);
}

/* ----------------------------------------- */

void Q1Lps2d::Structure(SparseStructureInterface* SI) const
{
  S.Structure(SI);
}


/* ----------------------------------------- */

void Q1Lps2d::StabForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  S.Form(f,u,EQ,d);
}

/* ----------------------------------------- */

void Q1Lps2d::Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const
{
  Q12d::Form(f,u,EQ,d);
  S    .Form(f,u,EQ,d);
}

/* ----------------------------------------- */

void Q1Lps2d::Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const
{
  Q12d::Matrix(A,u,EQ,d);
  S    .Matrix(A,u,EQ,d);
}
}
