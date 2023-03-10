/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2020 by the
 *Gascoigne 3D authors
 *
 * This file is part of Gascoigne 3D
 *
 * Gascoigne 3D is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version.
 *
 * Gascoigne 3D is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * Please refer to the file LICENSE.TXT for further information
 * on this license.
 *
 **/

#include <iomanip>
#include <list>

#include "stdsolver.h"
#include <numeric>

#include "pointilu.h"
#include "pointmatrix.h"

#include "dynamicblockilu.h"
#include "dynamicblockmatrix.h"

#include "vankasmoother.h"

/*--------------------------------*/
#include "cfdblock3d.h"
#include "fmatrixblock.h"
#include "sparseblockilu.h"

/*--------------------------------*/
#ifdef __WITH_UMFPACK__
#include "sparse_umf.h"
#include "umfilu.h"
#endif
/*--------------------------------*/

#include "backup.h"
#include "cuthillmckee.h"
#include "ilupermutate.h"
#include "pi.h"
#include "stopwatch.h"
#include "visu_eps.h"

#include "diracrighthandside.h"

//////////////////// Discretization
#include "baseq12d.h"
#include "baseq13d.h"
#include "baseq1patch.h"
#include "baseq22d.h"
#include "baseq23d.h"
#include "cgdisc.h"
#include "cgmixeddisc.h"
#include "lagrangedisc.h"

#include "elementintegrator.h"
#include "elementlpsintegrator.h"
#include "finiteelement.h"
#include "patchintegrationformula.h"
#include "transformation2d.h"
#include "transformation3d.h"

#include "glsequation.h"
#include "lpsequation.h"

using namespace std;

namespace Gascoigne {
extern Timer GlobalTimer;

StdSolver::StdSolver()
  : _MP(NULL)
  , _HM(NULL)
  , _ZP(NULL)
  , _PDX(NULL)
  //      , _NI(NULL)
  , _distribute(true)
  , _ndirect(1000)
  , _directsolver(0)
  , _discname("none")
  , _matrixtype("point_node")
  , _PrimalSolve(1)
  , _useUMFPACK(true)
// , omega_domain(0.)
{
}

/*-----------------------------------------*/

StdSolver::~StdSolver()
{
  if (_ZP)
    delete _ZP;
  _ZP = NULL;
}

/*-------------------------------------------------------*/

void
StdSolver::_check_consistency(const Equation* EQ,
                              const DiscretizationInterface* DI) const
{
  return;

  string eq = DI->GetName();

  bool glseq = false, glsdi = false;

  if (dynamic_cast<const GlsEquation*>(EQ)) {
    glseq = true;
  }
  if (eq == "Q1Gls2d" || eq == "Q2Gls2d" || eq == "Q1Gls3d" ||
      eq == "Q2Gls3d") {
    glsdi = true;
  }
  if (glseq && !glsdi) {
    cerr << "Warning: Discretization \"" << eq
         << "\" doesn't go with type of given Equation!" << endl;
  } else if (!glseq && glsdi) {
    cerr << "Error: Discretization \"" << eq
         << "\" doesn't go with type of given Equation!" << endl;
    abort();
  }

  bool lpseq = false, lpsdi = false;

  if (dynamic_cast<const LpsEquation*>(EQ)) {
    lpseq = true;
  }
  if (eq == "Q1Lps2d" || eq == "Q2Lps2d" || eq == "Q1Lps3d" ||
      eq == "Q2Lps3d") {
    lpsdi = true;
  }

  if (lpseq && !lpsdi) {
    cerr << "Warning: Discretization \"" << eq
         << "\" doesn't go with type of given Equation!" << endl;
  } else if (!lpseq && lpsdi) {
    cerr << "Error: Discretization \"" << eq
         << "\" doesn't go with type of given Equation!" << endl;
    abort();
  }
}

/*-------------------------------------------------------*/

void
StdSolver::MatrixZero(Matrix& A) const
{
  GetMatrix(A).zero();
}

/*-------------------------------------------------------*/

void
StdSolver::OutputSettings() const
{
  cout << "==================================================" << endl;
  cout << "Solver:                   " << GetName() << endl;
  cout << "Discretization:           " << GetDiscretization()->GetName()
       << endl;
  GetProblemDescriptor()->OutputSettings(cout);
  cout << "==================================================" << endl;
}

/*-------------------------------------------------------*/

void
StdSolver::SetProblem(const ProblemDescriptorInterface& PDX)
{
  _PDX = &PDX;
  assert(_PDX);

  const Equation* EQ = GetProblemDescriptor()->GetEquation();

  if (EQ)
    _check_consistency(EQ, GetDiscretization());
  _PF.SetComponents(GetSolverData().GetPfilter());
}

/*-------------------------------------------------------*/

void
StdSolver::SetDiscretization(DiscretizationInterface& DI, bool init)
{
  if (init) {
    DI.ReInit(GetMesh());
    DI.SetDataContainer(GetDiscretization()->GetDataContainer());
  }

  GetDiscretizationPointer() = &DI;
}

/*-------------------------------------------------------*/

/*-------------------------------------------------------*/

void
StdSolver::NewMesh(const GascoigneMesh* mp)
{
  _MP = mp;
  assert(_MP);

  if (_MP->nnodes() < _ndirect) {
    _directsolver = 1;
  } else {
    _directsolver = 0;
  }
  GetDiscretization()->ReInit(_MP);

  // reinit the pressure filter
  GetDiscretization()->InitFilter(GetPfilter());
}

/*-----------------------------------------*/

void
StdSolver::SetDefaultValues(string discname, string matrixtype, int ndirect)
{
  _discname = discname;
  _matrixtype = matrixtype;
  _ndirect = ndirect;
}

/*-------------------------------------------------------*/

void
StdSolver::BasicInit(const ParamFile& paramfile, const int dimension)
//                            const NumericInterface *NI)
{
  _paramfile = paramfile;
  //    _NI = NI;

  string xxx;

  DataFormatHandler DFH;
  DFH.insert("matrixtype", &_matrixtype);
  DFH.insert("ndirect", &_ndirect);
  DFH.insert("useUMFPACK", &_useUMFPACK);
  DFH.insert("discname", &_discname);
  DFH.insert("disc", &xxx, "void");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile, "Solver");

  if (xxx != "void") {
    cout << "Expression 'disc' in ParamFile not longer valid !" << endl;
    abort();
  }
  if (GetDiscretizationPointer() == NULL)
    GetDiscretizationPointer() = NewDiscretization(dimension, _discname);
  assert(_ZP);

  GetDiscretization()->BasicInit(_paramfile);
}

/*-------------------------------------------------------*/

DiscretizationInterface*
StdSolver::NewDiscretization(int dimension, const string& discname)
{
  // if (_NI)
  // {
  //   return _NI->NewDiscretization();
  // }
  // else
  // {
  // }
  std::cout << discname << std::endl;

  if (dimension == 2) {
    if (discname == "CGQ1")
      return new CGDiscQ12d;
    else if (discname == "CGQ2")
      return new CGDiscQ22d;

    else if (discname == "CGQ1Patch")
      return new CGDiscQ12dPatch;

    else if (discname == "CGQ1Mixed")
      return new CGMixedDiscQ12dPatch;

    else if (discname == "CGP1")
      return new CGDiscP12d;

    else if (discname == "CGQ1Lps")
      return new CGDiscQ12dLps;
    else if (discname == "CGQ2Lps")
      return new CGDiscQ22dLps;

    else if (discname == "LQ1")
      return new LagrangeDiscQ12d;
    else if (discname == "LQ2")
      return new LagrangeDiscQ22d;
    else if (discname == "LQ4")
      return new LagrangeDiscQ42d;

    else {
      cerr << " Solver::NewDiscretization()\tunknown discname=" << discname
           << endl;
      abort();
    }
  } else if (dimension == 3) {
    if (discname == "CGQ1")
      return new CGDiscQ13d;
    else if (discname == "CGQ2")
      return new CGDiscQ23d;
    else if (discname == "CGQ1Lps")
      return new CGDiscQ13dLps;
    else if (discname == "CGQ2Lps")
      return new CGDiscQ23dLps;

    else if (discname == "LQ1")
      return new LagrangeDiscQ13d;
    else if (discname == "LQ2")
      return new LagrangeDiscQ23d;
    else if (discname == "LQ4")
      return new LagrangeDiscQ43d;

    else {
      cerr << " Solver::NewDiscretization()\tunknown discname=" << discname
           << endl;
      abort();
    }
  } else {
    cerr << " Solver::NewDiscretization()\tdimension must either 2 or 3"
         << endl;
    abort();
  }
}

/*-------------------------------------------------------------*/

MatrixInterface*
StdSolver::NewMatrix(int ncomp, const string& matrixtype)
{
  if (_directsolver || matrixtype == "point_node") {
    return new PointMatrix(ncomp, "node");
  } else if ((matrixtype == "block") || (matrixtype == "sparseumf") ||
             (matrixtype == "vanka")) {
    if (ncomp == 1)
      return new SparseBlockMatrix<FMatrixBlock<1>>;
    else if (ncomp == 2)
      return new SparseBlockMatrix<FMatrixBlock<2>>;
    else if (ncomp == 3)
      return new SparseBlockMatrix<FMatrixBlock<3>>;
    else if (ncomp == 4)
      return new SparseBlockMatrix<FMatrixBlock<4>>;
    else if (ncomp == 5)
      return new SparseBlockMatrix<FMatrixBlock<5>>;
    else if (ncomp == 6)
      return new SparseBlockMatrix<FMatrixBlock<6>>;
    else {
      cerr << "No SparseBlockMatrix for " << ncomp << "components." << endl;
      abort();
    }
  } else if (matrixtype == "dynamic") {
    if (ncomp == 1)
      return new DynamicBlockMatrix<FMatrixBlock<1>>;
    else if (ncomp == 2)
      return new DynamicBlockMatrix<FMatrixBlock<2>>;
    else if (ncomp == 3)
      return new DynamicBlockMatrix<FMatrixBlock<3>>;
    else if (ncomp == 4)
      return new DynamicBlockMatrix<FMatrixBlock<4>>;
    else {
      cerr << "No SparseBlockMatrix for " << ncomp << "components." << endl;
      abort();
    }
  } else if (matrixtype == "component") {
    return new PointMatrix(ncomp, "component");
  } else if (matrixtype == "cfd") {
    if (ncomp == 4) {
      abort();
      //        return new SparseBlockMatrix<CFDBlock3d>;
    } else {
      cerr << "No SparseBlockMatrix for " << ncomp << "components." << endl;
      abort();
    }
  } else {
    cerr << "No such matrix type \"" << matrixtype << "\"." << endl;
    abort();
  }
}

/*-------------------------------------------------------------*/

IluInterface*
StdSolver::NewIlu(const Matrix& A, int ncomp, const string& matrixtype)
{
#ifdef __WITH_UMFPACK__
  if (_directsolver && _useUMFPACK) {
#ifdef __WITH_UMFPACK_LONG__
    return new UmfIluLong(&GetMatrix(A));
#else
    return new UmfIlu(&GetMatrix(A));
#endif
  }
#endif

  // analog zu NewMatrix muss hier auch _directsolver eingehen,
  // sonst gibts aerger nachher beim
  // GetIlu()->copy_entries(GetMatrix());
  if (_directsolver || matrixtype == "point_node") {
    return new PointIlu(ncomp, "node");
  }

  else if (matrixtype == "block") {
    if (ncomp == 1)
      return new SparseBlockIlu<FMatrixBlock<1>>;
    else if (ncomp == 2)
      return new SparseBlockIlu<FMatrixBlock<2>>;
    else if (ncomp == 3)
      return new SparseBlockIlu<FMatrixBlock<3>>;
    else if (ncomp == 4)
      return new SparseBlockIlu<FMatrixBlock<4>>;
    else if (ncomp == 5)
      return new SparseBlockIlu<FMatrixBlock<5>>;
    else if (ncomp == 6)
      return new SparseBlockIlu<FMatrixBlock<6>>;
    else {
      cerr << "No SparseBlockIlu for " << ncomp << "components." << endl;
      abort();
    }
  } else if (matrixtype == "vanka") {
    return new VankaSmoother;
  } else if (matrixtype == "sparseumf") {
#ifdef __WITH_UMFPACK__
    if (ncomp == 1)
      return new SparseUmf<FMatrixBlock<1>>(&GetMatrix(A));
    else if (ncomp == 2)
      return new SparseUmf<FMatrixBlock<2>>(&GetMatrix(A));
    else if (ncomp == 3)
      return new SparseUmf<FMatrixBlock<3>>(&GetMatrix(A));
    else if (ncomp == 4)
      return new SparseUmf<FMatrixBlock<4>>(&GetMatrix(A));
    else if (ncomp == 5)
      return new SparseUmf<FMatrixBlock<5>>(&GetMatrix(A));
    else
#endif
    {
      cerr << "No SparseBlockIlu for " << ncomp << "components." << endl;
      abort();
    }
  } else if (matrixtype == "dynamic") {
    if (ncomp == 1)
      return new DynamicBlockIlu<FMatrixBlock<1>>;
    else if (ncomp == 2)
      return new DynamicBlockIlu<FMatrixBlock<2>>;
    else if (ncomp == 3)
      return new DynamicBlockIlu<FMatrixBlock<3>>;
    else if (ncomp == 4)
      return new DynamicBlockIlu<FMatrixBlock<4>>;
    else {
      cerr << "No DynamicBlockIlu for " << ncomp << "components." << endl;
      abort();
    }
  } else if (matrixtype == "component") {
    return new PointIlu(ncomp, "component");
  } else if (matrixtype == "cfd") {
    if (ncomp == 4) {
      abort();
      //        return new SparseBlockIlu<CFDBlock3d>;
    }
  }
  cerr << "No such matrix type \"" << matrixtype << "and ncomp \"." << ncomp
       << endl;
  abort();
}

/*-------------------------------------------------------*/

void
StdSolver::AddPeriodicNodes(SparseStructure* SA)
{
  /*-------------------------------------------------------
  | Vereinigt die Nachbarschaften von zusammengehoerenden
  | Knoten auf den periodischen Raendern, so dass es
  | moeglich wird, die Kopplungen des einen Knoten zum
  | anderen zu schieben.
  -------------------------------------------------------*/

  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntVector& iv_PeriodicColors = BM->GetPeriodicDataColors();

  const GascoigneMesh* p_mesh = GetMesh();
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(p_mesh);
  assert(GMP);

  map<int, map<int, int>> mm_PeriodicPairs =
    GMP->GetBoundaryIndexHandler().GetPeriodicPairs();

  for (IntVector::const_iterator p_col = iv_PeriodicColors.begin();
       p_col != iv_PeriodicColors.end();
       p_col += 2) {
    int col = *p_col;

    IntSet is_neighbours1;
    IntSet is_neighbours2;

    for (map<int, int>::const_iterator p_pair = mm_PeriodicPairs[col].begin();
         p_pair != mm_PeriodicPairs[col].end();
         p_pair++) {
      // beide raender abgrasen und die kopplungen in columns1 und columns2
      // eintragen
      is_neighbours1 = SA->row(p_pair->first);
      is_neighbours2 = SA->row(p_pair->second);

      for (IntSet::const_iterator p_neighbour1 = is_neighbours1.begin();
           p_neighbour1 != is_neighbours1.end();
           p_neighbour1++) {
        SA->build_add(p_pair->second, *p_neighbour1);
        SA->build_add(*p_neighbour1, p_pair->second);
      }

      for (IntSet::const_iterator p_neighbour2 = is_neighbours2.begin();
           p_neighbour2 != is_neighbours2.end();
           p_neighbour2++) {
        SA->build_add(p_pair->first, *p_neighbour2);
        SA->build_add(*p_neighbour2, p_pair->first);
      }
    }
  }
  SA->build_end();
}

/*-------------------------------------------------------*/

// Access to Vector & Matrix Data
GlobalVector&
StdSolver::GetGV(Vector& u) const
{
  return vector_agent(u);
}
const GlobalVector&
StdSolver::GetGV(const Vector& u) const
{
  return vector_agent(u);
}

/*-------------------------------------------------------*/

void
StdSolver::ReInitVector(Vector& dst)
{
  IndexType size = 1;
  if (dst.GetType() == "node") {
    size = GetDiscretization()->ndofs();
  } else if (dst.GetType() == "cell") {
    size = GetDiscretization()->nelements();
  } else if (dst.GetType() == "parameter") {
    size = 1;
  } else {
    throw std::runtime_error(std::string("Can not initialize vector of type ") +
                             dst.GetType());
  }

  IndexType comp = GetProblemDescriptor()->GetNcomp();
  auto& p = vector_agent[dst];

  // If not allready in Vector Agent reserve new
  if (p == NULL) {
    p = new GlobalVector(comp, size, 0);
  } else {
    p->ncomp() = comp;
    p->reservesize(size);
    p->zero();
  }
}

/*---------------------------------------------------*/

void
StdSolver::DeleteVector(Vector& p) const
{
  vector_agent.Delete(p);
}

/*-------------------------------------------------------*/

void
StdSolver::ReInitMatrix(const Matrix& A)
{
  GlobalTimer.start("---> matrix init");

  IndexType ncomp = GetProblemDescriptor()->GetNcomp();

  auto& matrix = matrix_agent[A];

  // check if matrix is an umfpack object. If yes, delete it
#ifdef __WITH_UMFPACK__
  if (_useUMFPACK && matrix != NULL) {
    SimpleMatrix* SM = dynamic_cast<SimpleMatrix*>(matrix);
    if ((SM && !_directsolver && _matrixtype != "point_node") ||
        (!SM && _directsolver)) {
      delete matrix;
      matrix = NULL;
    }
  }
#endif

  if (matrix == NULL) {
    matrix = NewMatrix(ncomp, _matrixtype);
  }

  auto& ilu = ilu_agent[A];

  // if umfpack is used and matrix is umfpack object, delete it
#ifdef __WITH_UMFPACK__
  if (_useUMFPACK && ilu != NULL) {
#ifdef __WITH_UMFPACK_LONG__
    UmfIluLong* UM = dynamic_cast<UmfIluLong*>(ilu);
#else
    UmfIlu* UM = dynamic_cast<UmfIlu*>(ilu);
#endif

    if ((UM && !_directsolver) || (!UM && _directsolver)) {
      delete ilu;
      ilu = NULL;
    }
  }
#endif

  if (ilu == NULL)
    ilu = NewIlu(A, ncomp, _matrixtype);

  // if Vankasmoother is used, attach dofhandler
  if (!_directsolver && (_matrixtype == "vanka")) {
    assert(dynamic_cast<const VankaSmoother*>(ilu));
    dynamic_cast<const VankaSmoother*>(ilu)->SetDofHandler(GetMesh());
  }

  ////////// setup the stencil
  // is it a drawback that the stencil cannot be reused for multiple matrices?

  SparseStructure SA;
  GetDiscretization()->Structure(&SA);
  AddPeriodicNodes(&SA);

  matrix->ReInit(&SA);

  if (ilu != NULL)
    ilu->ReInit(&SA);

  GlobalTimer.stop("---> matrix init");
}

/*-------------------------------------------------------*/

void
StdSolver::Zero(Vector& dst) const
{
  GetGV(dst).zero();
}

/*-----------------------------------------*/

double
StdSolver::NewtonNorm(const Vector& u) const
{
  return GetGV(u).norm_l8();
}

/*-----------------------------------------*/

void
StdSolver::HNAverageData() const
{
  GetDiscretization()->HNAverageData();
}
void
StdSolver::HNZeroData() const
{
  GetDiscretization()->HNZeroData();
}
void
StdSolver::HNAverage(Vector& x) const
{
  GetDiscretization()->HNAverage(GetGV(x));
}
void
StdSolver::HNZero(Vector& x) const
{
  GetDiscretization()->HNZero(GetGV(x));
}
void
StdSolver::HNDistribute(Vector& x) const
{
  if (GetDistribute()) {
    GetDiscretization()->HNDistribute(GetGV(x));
  }
}

/*-------------------------------------------------------*/

void
StdSolver::InterpolateSolution(Vector& gu, const GlobalVector& uold) const
{
  GlobalVector& u = GetGV(gu);

  u.zero();
  GetDiscretization()->InterpolateSolution(u, uold);
  SubtractMean(gu);
}

/*-----------------------------------------*/

void
StdSolver::residualgmres(const Matrix& A,
                         Vector& gy,
                         const Vector& gx,
                         const Vector& gb) const
{
  GlobalVector& y = GetGV(gy);
  const GlobalVector& b = GetGV(gb);

  vmulteq(A, gy, gx, 1.);
  y.sadd(-1., 1., b);
  SetBoundaryVectorZero(gy);
}

/*-----------------------------------------*/

void
StdSolver::vmult(const Matrix& A, Vector& gy, const Vector& gx, double d) const
{
  GlobalTimer.start("---> vmult");
  GetMatrix(A).vmult(GetGV(gy), GetGV(gx), d);
  GlobalTimer.stop("---> vmult");
}

/*-----------------------------------------*/

void
StdSolver::vmulteq(const Matrix& A,
                   Vector& gy,
                   const Vector& gx,
                   double d) const
{
  GlobalTimer.start("---> vmult");
  Zero(gy);
  GlobalTimer.stop("---> vmult");
  vmult(A, gy, gx, d);
}

/*-----------------------------------------*/

void
StdSolver::MatrixResidual(const Matrix& A,
                          Vector& gy,
                          const Vector& gx,
                          const Vector& gb) const
{
  Equ(gy, 1., gb);
  vmult(A, gy, gx, -1.);
  SubtractMeanAlgebraic(gy);
}

void
StdSolver::Jacobi(const Matrix& A, Vector& y) const
{
  GetMatrix(A).Jacobi(GetGV(y));
}

/*-------------------------------------------------------*/

void
StdSolver::SetBoundaryVectorZero(Vector& gf) const
{
  GetDiscretization()->StrongDirichletVectorZero(GetGV(gf),
                                                 *GetProblemDescriptor());

  // GlobalVector &f = GetGV(gf);

  // const BoundaryManager *BM = GetProblemDescriptor()->GetBoundaryManager();
  // const IntSet &Colors = BM->GetDirichletDataColors();

  // for (IntSet::const_iterator p = Colors.begin(); p != Colors.end(); p++)
  // {
  //   int col = *p;
  //   const IntVector &comp = BM->GetDirichletDataComponents(col);
  //   GetDiscretization()->StrongDirichletVectorZero(f, col, comp);
  //    }
}

/*-------------------------------------------------------*/

void
StdSolver::SetBoundaryVector(Vector& gf) const
{
  GetDiscretization()->StrongDirichletVector(
    GetGV(gf), GetProblemDescriptor()->GetDirichletData());
  // const BoundaryManager *BM = GetProblemDescriptor()->GetBoundaryManager();
  // const DirichletData *DD = GetProblemDescriptor()->GetDirichletData();
  // if (DD == NULL)
  // {
  //   if (BM->GetDirichletDataColors().size() != 0)
  //   {
  //     cerr << "No DirichetData given but DirichetColors in ParamFile!"
  //          << endl;
  //     abort();
  //   }
  //   return;
  // }
  // SetBoundaryVectorStrong(gf, *BM, *DD);
}

/*-------------------------------------------------------*/

void
StdSolver::SetPeriodicVector(Vector& gf) const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const PeriodicData* PD = GetProblemDescriptor()->GetPeriodicData();

  if (PD == NULL) {
    if (BM->GetPeriodicDataColors().size() != 0) {
      cerr << "No PeriodicData given but PeriodicColors in ParamFile!" << endl;
      abort();
    }
    return;
  }

  SetPeriodicVectorStrong(gf, *BM, *PD);
}

/*-------------------------------------------------------*/

void
StdSolver::SetBoundaryVectorStrong(Vector& gf,
                                   const BoundaryManager& BM,
                                   const DirichletData& DD,
                                   double d) const
{
  cerr << "StdSolver::SetBoundaryVectorStrong(...). Old Interface. Not used "
          "any more"
       << endl
       << "\t New interface has colors in Dirichlet Data. Directly called in "
          "Discretization"
       << endl;
  abort();
  // GlobalVector &f = GetGV(gf);

  // IntSet PrefCol = DD.preferred_colors();
  // list<int> colors(BM.GetDirichletDataColors().begin(),
  //                  BM.GetDirichletDataColors().end());

  // for (IntSet::const_iterator p = PrefCol.begin(); p != PrefCol.end(); p++)
  // {
  //   int col = *p;
  //   colors.remove(col);
  //   colors.push_back(col);
  // }
  // for (list<int>::const_iterator p = colors.begin(); p != colors.end(); p++)
  // {
  //   int col = *p;
  //   const IntVector &comp = BM.GetDirichletDataComponents(col);
  //   GetDiscretization()->StrongDirichletVector(f, DD, col, comp, d);
  // }
}

/*-------------------------------------------------------*/

void
StdSolver::SetPeriodicVectorStrong(Vector& gf,
                                   const BoundaryManager& BM,
                                   const PeriodicData& PD,
                                   double d) const
{
  GlobalVector& f = GetGV(gf);

  list<int> periodic_cols(BM.GetPeriodicDataColors().begin(),
                          BM.GetPeriodicDataColors().end());
  for (list<int>::const_iterator p = periodic_cols.begin();
       p != periodic_cols.end();
       p++) {
    int col = *p;
    const IntVector& comp = BM.GetPeriodicDataComponents(col);
    GetDiscretization()->StrongPeriodicVector(f, PD, col, comp, d);
  }
}

/*-------------------------------------------------------*/

void
StdSolver::SetPeriodicVectorZero(Vector& gf) const
{
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntVector& iv_PeriodicColors = BM->GetPeriodicDataColors();

  GlobalVector& f = GetGV(gf);
  const GascoigneMesh* p_mesh = GetMesh();
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(p_mesh);
  assert(GMP);

  map<int, map<int, int>> mm_PeriodicPairs =
    GMP->GetBoundaryIndexHandler().GetPeriodicPairs();

  for (IntVector::const_iterator p_col = iv_PeriodicColors.begin();
       p_col != iv_PeriodicColors.end();) {
    int col = *p_col++;
    *p_col++;

    const IntVector& iv_PeriodicComponents = BM->GetPeriodicDataComponents(col);

    for (map<int, int>::const_iterator p_pair = mm_PeriodicPairs[col].begin();
         p_pair != mm_PeriodicPairs[col].end();
         p_pair++) {
      for (IntVector::const_iterator p_comp = iv_PeriodicComponents.begin();
           p_comp != iv_PeriodicComponents.end();
           p_comp++) {
        f(p_pair->second, *p_comp) =
          .5 * (f(p_pair->second, *p_comp) + f(p_pair->first, *p_comp));
        f(p_pair->first, *p_comp) = f(p_pair->second, *p_comp);
      }
    }
  }
}

/*-------------------------------------------------------*/

void
StdSolver::smooth(int niter,
                  const Matrix& A,
                  Vector& x,
                  const Vector& y,
                  Vector& h) const
{

  GlobalTimer.start("---> smooth");
  double omega = GetSolverData().GetOmega();

  for (int iter = 0; iter < niter; iter++) {
    if (GetSolverData().GetLinearSmooth() == "ilu") {
      GlobalTimer.stop("---> smooth");
      MatrixResidual(A, h, x, y);
      GlobalTimer.start("---> smooth");
      GetIlu(A).solve(GetGV(h));
      Add(x, omega, h);
    } else if (GetSolverData().GetLinearSmooth() == "jacobi") {
      GlobalTimer.stop("---> smooth");
      MatrixResidual(A, h, x, y);
      GlobalTimer.start("---> smooth");
      Jacobi(A, h);
      Add(x, omega, h);
    } else if (GetSolverData().GetLinearSmooth() == "richardson") {
      MatrixResidual(A, h, x, y);
      Add(x, omega, h);
    } else if (GetSolverData().GetLinearSmooth() == "none") {
    } else {
      cerr << "Smoother: " << GetSolverData().GetLinearSmooth()
           << " not valid!\n";
      abort();
    }
    SubtractMean(x);
  }
  GlobalTimer.stop("---> smooth");
}

/*-------------------------------------------------------*/

void
StdSolver::smooth_pre(const Matrix& A,
                      Vector& x,
                      const Vector& y,
                      Vector& help) const
{
  int niter = GetSolverData().GetIterPre();
  smooth(niter, A, x, y, help);
}

/*-------------------------------------------------------*/

void
StdSolver::smooth_exact(const Matrix& A,
                        Vector& x,
                        const Vector& y,
                        Vector& help) const
{
#ifdef __WITH_UMFPACK__
  if (_directsolver && _useUMFPACK) {
    GlobalTimer.start("---> smooth (ex)");
#ifdef __WITH_UMFPACK_LONG__
    const UmfIluLong* UM = dynamic_cast<const UmfIluLong*>(&GetIlu(A));
#else
    const UmfIlu* UM = dynamic_cast<const UmfIlu*>(&GetIlu(A));
#endif
    assert(UM);
    UM->Solve(GetGV(x), GetGV(y));
    GlobalTimer.stop("---> smooth (ex)");
  } else
#endif
  {
    int niter = GetSolverData().GetIterExact();
    smooth(niter, A, x, y, help);
  }
}

/*-------------------------------------------------------*/

void
StdSolver::smooth_post(const Matrix& A,
                       Vector& x,
                       const Vector& y,
                       Vector& help) const
{
  int niter = GetSolverData().GetIterPost();
  smooth(niter, A, x, y, help);
}

/*-------------------------------------------------------*/

void
StdSolver::Form(Vector& gy, Vector& gx, double d) const
{
  GlobalTimer.start("---> form");

  HNAverage(gx);
  HNAverageData();

  //////////// Form in elements
  GetDiscretization()->Form(
    GetGV(gy), GetGV(gx), *GetProblemDescriptor()->GetEquation(), d);

  //////////// Boundary
  GetDiscretization()->BoundaryForm(
    GetGV(gy), GetGV(gx), *GetProblemDescriptor(), d);

  HNZero(gx);
  HNZeroData();
  HNDistribute(gy);
  SubtractMeanAlgebraic(gy);

  GlobalTimer.stop("---> form");
}

/*-------------------------------------------------------*/

void
StdSolver::AdjointForm(Vector& gy, Vector& gx, double d) const
{
  GlobalVector& y = GetGV(gy);
  const GlobalVector& x = GetGV(gx);

  HNAverage(gx);
  HNAverageData();

  const Equation* EQ = GetProblemDescriptor()->GetEquation();
  assert(EQ);
  GetDiscretization()->AdjointForm(y, x, *EQ, d);

  abort();
  const BoundaryEquation* BE = GetProblemDescriptor()->GetBoundaryEquation();
  if (BE) {
    GetDiscretization()->BoundaryForm(y, x, *GetProblemDescriptor(), d);
  }

  HNZero(gx);
  HNZeroData();
  HNDistribute(gy);
  SubtractMeanAlgebraic(gy);
}

/*-------------------------------------------------------*/

void
StdSolver::BoundaryInit(Vector& Gu) const
{
  GlobalVector& u = GetGV(Gu);
  const DirichletData* DD = GetProblemDescriptor()->GetDirichletData();
  if (DD == NULL) {
    u.zero();
    cerr << "StdSolver::BoundaryInit():\t No DirichletData given but required "
            "for "
            "boundary init!!\n";
    cerr << "\t\t\t\t Alternative: use initial solution" << endl;
    cerr << "\t\t\t\t or: use Dirichlet data without dirichlet_colors" << endl;
    abort();
  }

  int color = 0;
  int ncomp = GetProblemDescriptor()->GetNcomp();
  DoubleVector y(ncomp);

  int dim = GetMesh()->dimension();

  for (int ind = 0; ind < GetDiscretization()->ndofs(); ind++) {
    if (dim == 2)
      (*DD)(y, GetDiscretization()->vertex2d(ind), color);
    else
      (*DD)(y, GetDiscretization()->vertex3d(ind), color);

    for (int c = 0; c < u.ncomp(); c++) {
      u(ind, c) = y[c];
    }
  }
}

/*-------------------------------------------------------*/

void
StdSolver::SolutionInit(Vector& Gu) const
{
  GlobalVector& u = GetGV(Gu);
  const DomainInitialCondition* u0 = dynamic_cast<const DomainRightHandSide*>(
    GetProblemDescriptor()->GetInitialCondition());
  if (u0 == NULL) {
    u.zero();
    return;
  }

  assert(u.ncomp() == GetProblemDescriptor()->GetNcomp());

  assert(0);
  for (int ind = 0; ind < GetMesh()->nnodes(); ind++) {
    if (GetMesh()->dimension() == 2) {
      for (int c = 0; c < u.ncomp(); c++) {
        u(ind, c) = (*u0)(c, GetMesh()->vertex2d(ind));
      }
    } else {
      for (int c = 0; c < u.ncomp(); c++) {
        u(ind, c) = (*u0)(c, GetMesh()->vertex3d(ind));
      }
    }
  }
}

/*-------------------------------------------------------*/

void
StdSolver::ComputeError(Vector& u, GlobalVector& err) const
{
  if (GetProblemDescriptor()->GetExactSolution() == NULL)
    return;
  HNAverage(u);
  GetDiscretization()->ComputeError(
    GetGV(u), err, GetProblemDescriptor()->GetExactSolution());
  HNZero(u);
}

/*-------------------------------------------------------*/

void
StdSolver::AssembleError(GlobalVector& eta, Vector& u, GlobalVector& err) const
{
  if (GetProblemDescriptor()->GetExactSolution() == NULL)
    return;
  HNAverage(u);
  GetDiscretization()->AssembleError(
    eta, GetGV(u), err, GetProblemDescriptor()->GetExactSolution());
  HNZero(u);
}

/*-------------------------------------------------------*/

double
StdSolver::ComputeFunctional(Vector& gf, Vector& gu, const Functional* FP)
{
  Vector gh("XXX");
  ReInitVector(gh);
  double val = 0.;

  const DomainFunctional* DFP = dynamic_cast<const DomainFunctional*>(FP);
  const BoundaryFunctional* BFP = dynamic_cast<const BoundaryFunctional*>(FP);
  const ResidualFunctional* RFP = dynamic_cast<const ResidualFunctional*>(FP);
  const PointFunctional* NPFP = dynamic_cast<const PointFunctional*>(FP);

  if (DFP)
    val = ComputeDomainFunctional(gu, DFP);
  else if (BFP)
    val = ComputeBoundaryFunctional(gf, gu, gh, BFP);
  else if (RFP)
    val = ComputeResidualFunctional(gf, gu, gh, RFP);
  else if (NPFP)
    val = ComputePointFunctional(gf, gu, gh, NPFP);
  else {
    cerr << "Functional must be either of type DomainFunctional, "
            "BoundaryFunctional or PointFunctional!!!"
         << endl;
    abort();
  }
  DeleteVector(gh);
  return val;
}

/*-------------------------------------------------------*/

double
StdSolver::ComputeBoundaryFunctional(Vector& gf,
                                     Vector& gu,
                                     Vector& gz,
                                     const BoundaryFunctional* FP) const
{
  HNAverage(gu);
  HNAverageData();
  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  double J = GetDiscretization()->ComputeBoundaryFunctional(
    GetGV(gu), BM->GetBoundaryFunctionalColors(), *FP);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double
StdSolver::ComputeDomainFunctional(Vector& gu, const DomainFunctional* FP) const
{
  HNAverage(gu);
  HNAverageData();
  double J = GetDiscretization()->ComputeDomainFunctional(GetGV(gu), *FP);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double
StdSolver::ComputePointFunctional(Vector& gf,
                                  Vector& gu,
                                  Vector& gz,
                                  const PointFunctional* FP) const
{
  HNAverage(gu);
  HNAverageData();
  double J = GetDiscretization()->ComputePointFunctional(GetGV(gu), *FP);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

double
StdSolver::ComputeResidualFunctional(Vector& gf,
                                     Vector& gu,
                                     Vector& gz,
                                     const ResidualFunctional* FP) const
{
  // cerr << "Aenderung in ResidualFunctional" << endl;

  // //    const BoundaryManager *BM =
  // GetProblemDescriptor()->GetBoundaryManager();

  // HNAverage(gu);
  // HNAverageData();
  // Zero(gf);
  // Rhs(gf);
  // // Do not include Boundary Equation!
  // //    Form(gf, gu, -1.); >>>>>>>>
  // const Equation *EQ = GetProblemDescriptor()->GetEquation();
  // GetDiscretization()->Form(GetGV(gf), GetGV(gu), *EQ, -1.0);
  // //    Form(gf, gu, -1.); <<<<<<<<

  // const DirichletData *ABD = FP->GetDirichletData();
  // assert(ABD);

  // Zero(gz);

  // // modification needed when funcitonal is not on dirichlet boundary
  // //    SetBoundaryVectorStrong(gz, *BM, *ABD); >>>>>

  // IntSet       PrefCol = ABD->preferred_colors();
  // nvector<int> comps   = FP->GetComps();

  // for (auto it : PrefCol)
  //   GetDiscretization()->StrongDirichletVector(GetGV(gz), *ABD, it,
  //   comps,1.0);
  // //    SetBoundaryVectorStrong(gz, *BM, *ABD); <<<<
  // HNAverage(gz);

  // double J = ScalarProduct(gz, gf);
  // HNZero(gu);
  // HNZeroData();
  // return J;

  HNAverage(gu);
  HNAverageData();
  Zero(gf);
  Rhs(gf);
  Form(gf, gu, -1.);

  const DirichletData* ABD = FP->GetDirichletData();
  assert(ABD);

  Zero(gz);
  GetDiscretization()->StrongDirichletVector(GetGV(gz), ABD);

  HNAverage(gz);

  double J = ScalarProduct(gz, gf);
  HNZero(gu);
  HNZeroData();
  return J;
}

/*-------------------------------------------------------*/

void
StdSolver::EvaluateCellRightHandSide(Vector& f,
                                     const DomainRightHandSide& CF,
                                     double d) const
{
  assert(f.GetType() == "cell");
  HNAverageData();

  GetDiscretization()->EvaluateCellRightHandSide(GetGV(f), CF, d);

  HNZeroData();
}

/*-------------------------------------------------------*/

void
StdSolver::EvaluateBoundaryCellRightHandSide(Vector& f,
                                             const BoundaryRightHandSide& CF,
                                             const BoundaryManager& BM,
                                             double d) const
{
  assert(f.GetType() == "cell");
  HNAverageData();

  GetDiscretization()->EvaluateBoundaryCellRightHandSide(
    GetGV(f), BM.GetBoundaryRightHandSideColors(), CF, d);

  HNZeroData();
}

/*-------------------------------------------------------*/

void
StdSolver::EvaluateParameterRightHandSide(Vector& f,
                                          const DomainRightHandSide& CF,
                                          double d) const
{
  assert(f.GetType() == "parameter");
  HNAverageData();

  GetDiscretization()->EvaluateParameterRightHandSide(GetGV(f), CF, d);

  HNZeroData();
}

/*-------------------------------------------------------*/

void
StdSolver::EvaluateBoundaryParameterRightHandSide(
  Vector& f,
  const BoundaryRightHandSide& CF,
  const BoundaryManager& BM,
  double d) const
{
  assert(f.GetType() == "parameter");
  HNAverageData();

  GetDiscretization()->EvaluateBoundaryParameterRightHandSide(
    GetGV(f), BM.GetBoundaryRightHandSideColors(), CF, d);

  HNZeroData();
}

/*-------------------------------------------------------*/

void
StdSolver::InterpolateDomainFunction(Vector& f, const DomainFunction& DF) const
{
  HNAverageData();

  if (f.GetType() == "node") {
    GetDiscretization()->InterpolateDomainFunction(GetGV(f), DF);
  } else if (f.GetType() == "cell") {
    GetDiscretization()->InterpolateCellDomainFunction(GetGV(f), DF);
  } else {
    cerr << "No such vector type: " << f.GetType() << endl;
    abort();
  }

  HNZeroData();
}

/*-------------------------------------------------------*/

void
StdSolver::Rhs(Vector& gf, double d) const
{
  const auto* RHS = GetProblemDescriptor()->GetRightHandSide();
  const auto* DRHS = GetProblemDescriptor()->GetDiracRightHandSide();
  const auto* BRHS = GetProblemDescriptor()->GetBoundaryRightHandSide();

  if ((RHS == NULL) && (BRHS == NULL) && (DRHS == NULL))
    return;

  GlobalVector& f = GetGV(gf);
  HNAverageData();

  if (DRHS) {
    GetDiscretization()->DiracRhs(f, *DRHS, d);
  }

  if (RHS) {
    bool done = false;
    const DomainRightHandSide* DRHS =
      dynamic_cast<const DomainRightHandSide*>(RHS);
    if (DRHS) {
      GetDiscretization()->Rhs(f, *DRHS, d);
      done = true;
    }
    const DiracRightHandSide* NDRHS =
      dynamic_cast<const DiracRightHandSide*>(RHS);
    if (NDRHS) {
      abort();
      GetDiscretization()->DiracRhs(f, *NDRHS, d);
      done = true;
    }
    if (!done) {
      cerr << "RightHandSide should be either of type DomainRightHandSide or "
              "DiracRightHandSide!!!"
           << endl;
      abort();
    }
  }

  if (BRHS) {
    assert(BRHS->GetNcomp() == f.ncomp());
    const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
    GetDiscretization()->BoundaryRhs(
      f, BM->GetBoundaryRightHandSideColors(), *BRHS, d);
  }

  HNZeroData();
  HNDistribute(gf);
}

/*-------------------------------------------------------*/

void
StdSolver::AssembleMatrix(Matrix& A, Vector& gu, double d) const
{
  GlobalTimer.start("---> matrix");

  const GlobalVector& u = GetGV(gu);
  HNAverage(gu);
  HNAverageData();

  //////////// Elements
  GetDiscretization()->Matrix(
    GetMatrix(A), u, *GetProblemDescriptor()->GetEquation(), d);

  //////////// Boundary
  GetDiscretization()->BoundaryMatrix(
    GetMatrix(A), u, *GetProblemDescriptor(), d);

  PeriodicMatrix(A);
  DirichletMatrix(A);
  HNZero(gu);
  HNZeroData();

  GlobalTimer.stop("---> matrix");
}

/*-------------------------------------------------------*/

void
StdSolver::DirichletMatrix(Matrix& A) const
{
  GetDiscretization()->StrongDirichletMatrix(GetMatrix(A),
                                             *GetProblemDescriptor());
}

/* -------------------------------------------------------*/

void
StdSolver::DirichletMatrixOnlyRow(Matrix& A) const
{
  GetDiscretization()->StrongDirichletMatrixOnlyRow(GetMatrix(A),
                                                    *GetProblemDescriptor());
}

/* -------------------------------------------------------*/

void
StdSolver::PeriodicMatrix(Matrix& A) const
{
  /*-------------------------------------------------------
  | Modifiziert die Systemmatrix, um den periodischen
  | Raendern Rechnung zu tragen.
  | Vgl. DirichletMatrix bzw. DirichletMatrixOnlyRow.
  | Ruft dazu Matrix->periodic() auf.
  -------------------------------------------------------*/

  const BoundaryManager* BM = GetProblemDescriptor()->GetBoundaryManager();
  const IntVector& iv_PeriodicColors = BM->GetPeriodicDataColors();

  const GascoigneMesh* p_mesh = GetMesh();
  const GascoigneMesh* GMP = dynamic_cast<const GascoigneMesh*>(p_mesh);
  assert(GMP);

  map<int, map<int, int>> mm_PeriodicPairs =
    GMP->GetBoundaryIndexHandler().GetPeriodicPairs();

  for (IntVector::const_iterator p_col = iv_PeriodicColors.begin();
       p_col != iv_PeriodicColors.end();) {
    int col = *p_col++;
    *p_col++;

    const IntVector iv_PeriodicComponents = BM->GetPeriodicDataComponents(col);

    GetMatrix(A).periodic(mm_PeriodicPairs[col], iv_PeriodicComponents);
  }
}

/* -------------------------------------------------------*/

void
StdSolver::ComputeIlu(Matrix& A, const Vector& gu) const
{
#ifdef __WITH_UMFPACK__
  if (_directsolver && _useUMFPACK) {
    GlobalTimer.start("---> direct");
#ifdef __WITH_UMFPACK_LONG__
    UmfIluLong* UM = dynamic_cast<UmfIluLong*>(&GetIlu(A));
#else
    UmfIlu* UM = dynamic_cast<UmfIlu*>(&GetIlu(A));
#endif
    assert(UM);
    //       if(PrimalSolve==0) return;
    UM->Factorize();
    GlobalTimer.stop("---> direct");
  } else
#endif
    if (GetSolverData().GetLinearSmooth() == "ilu") {
    GlobalTimer.start("---> ilu");
    int ncomp = GetProblemDescriptor()->GetNcomp();
    PermutateIlu(A, gu);
    GetIlu(A).zero();
    GetIlu(A).copy_entries(GetMatrix(A));
    modify_ilu(GetIlu(A), ncomp);
    GetIlu(A).compute_ilu();
    GlobalTimer.stop("---> ilu");
  }
}

/*-------------------------------------------------------*/

void
StdSolver::modify_ilu(IluInterface& I, int ncomp) const
{
  if (GetSolverData().GetIluModify().size() == 0)
    return;
  if (GetSolverData().GetIluModify().size() != ncomp) {
    cerr << "ERROR: GetSolverData().GetIluModify().size()="
         << GetSolverData().GetIluModify().size() << " and ";
    cerr << "ncomp=" << ncomp << endl;
    abort();
    // assert(GetSolverData().GetIluModify().size()==ncomp);
  }

  for (int c = 0; c < ncomp; c++) {
    double s = GetSolverData().GetIluModify(c);
    I.modify(c, s);
  }
}

/* -------------------------------------------------------*/

void
StdSolver::PermutateIlu(Matrix& A, const Vector& gu) const
{
  const GlobalVector& u = GetGV(gu);

  int n = GetMatrix(A).GetStencil()->n();
  IntVector perm(n);

  iota(perm.begin(), perm.end(), 0);
  if (_matrixtype != "vanka") // no need to sort for vanka
  {
    if (GetSolverData().GetIluSort() == "cuthillmckee") {
      CuthillMcKee cmc(GetMatrix(A).GetStencil());
      cmc.Permutate(perm);
    } else if (GetSolverData().GetIluSort() == "streamdirection") {
      const Equation* EQ = GetProblemDescriptor()->GetEquation();
      (void)EQ;
      assert(EQ);
      assert(GetSolverData().GetStreamDirection().size() <= EQ->GetNcomp());
      StreamDirection sd(GetMesh(), GetMatrix(A).GetStencil(), u);
      sd.Permutate(perm, GetSolverData().GetStreamDirection());
    } else if (GetSolverData().GetIluSort() == "vectordirection") {
      VecDirection vd(GetMesh());
      vd.Permutate(perm, GetSolverData().GetVectorDirection());
    }
  }
  GetIlu(A).ConstructStructure(perm, GetMatrix(A));
}

/* -------------------------------------------------------*/

void
StdSolver::Visu(const string& name, const Vector& gu, int i) const
{
  GlobalTimer.start("---> visu");
  if (gu.GetType() == "node") {
    PointVisu(name, GetGV(gu), i);
  } else if (gu.GetType() == "cell") {
    CellVisu(name, GetGV(gu), i);
  } else {
    cerr << "No such vector type: " << gu.GetType() << endl;
    abort();
  }
  GlobalTimer.stop("---> visu");
}

/* -------------------------------------------------------*/

void
StdSolver::PointVisu(const string& name, const GlobalVector& u, int i) const
{
  GetDiscretization()->VisuVtk(
    GetProblemDescriptor()->GetComponentInformation(), _paramfile, name, u, i);
}

/* -------------------------------------------------------*/

void
StdSolver::CellVisu(const string& name, const GlobalVector& u, int i) const
{
  GascoigneVisualization Visu;

  Visu.SetMesh(GetMesh());

  const ComponentInformation* CI =
    GetProblemDescriptor()->GetComponentInformation();
  if (CI) {
    Visu.AddCellVector(CI, &u);
  } else {
    Visu.AddCellVector(&u);
  }

  Visu.read_parameters(_paramfile);
  Visu.set_name(name);
  Visu.step(i);
  Visu.write();
}

/* -------------------------------------------------------*/

void
StdSolver::VisuGrid(const string& name, int i) const
{
  assert(GetMesh());

  if (GetMesh()->dimension() == 2) {
    VisuEPS eps(_paramfile);
    //  eps.SetOption(VisuEPS::LINEWIDTH,0.1);
    if (_discname[1] == '2') {
      eps.SetOption(VisuEPS::WRITE_PATCH, 1);
    }
    eps.SetMesh(*GetMesh());
    eps.WriteGrid(name, i);
  }
}

/*-------------------------------------------------------*/

void
StdSolver::Read(Vector& gu, const string& filename) const
{
  GlobalVector& u = GetGV(gu);
  u.zero();
  ReadBackUp(u, filename);
}

/*-------------------------------------------------------*/

void
StdSolver::Write(const Vector& gu, const string& filename) const
{
  const GlobalVector& u = GetGV(gu);
  WriteBackUp(u, filename);
}

/*-------------------------------------------------------*/

void
StdSolver::ConstructInterpolator(MgInterpolatorInterface* I,
                                 const MeshTransferInterface* MT)
{
  GetDiscretization()->ConstructInterpolator(I, MT);
}

/* -------------------------------------------------------*/

DoubleVector
StdSolver::IntegrateSolutionVector(Vector& gu) const
{
  HNAverage(gu);
  DoubleVector dst = GetPfilter().IntegrateVector(GetGV(gu));
  HNZero(gu);
  return dst;
}

/* -------------------------------------------------------*/

void
StdSolver::SubtractMean(Vector& gx) const
{
  GlobalVector& x = GetGV(gx);
  // In each nonlinear step: applied to Newton correction,
  // in each smoothing step
  //
  if (GetPfilter().Active()) {
    GetDiscretization()->HNZeroCheck(x);
    GetPfilter().SubtractMean(x);
    HNZero(gx);
  }
}

/* -------------------------------------------------------*/

void
StdSolver::SubtractMeanAlgebraic(Vector& gx) const
{
  GlobalVector& x = GetGV(gx);

  // applies to residuals
  if (GetPfilter().Active()) {
    GetDiscretization()->HNZeroCheck(x);
    GetPfilter().SubtractMeanAlgebraic(x);
    HNZero(gx);
  }
}

/*-----------------------------------------*/

void
StdSolver::Equ(Vector& dst, double s, const Vector& src) const
{
  GetGV(dst).equ(s, GetGV(src));
}

/*-----------------------------------------*/

void
StdSolver::Add(Vector& dst, double s, const Vector& src) const
{
  GetGV(dst).add(s, GetGV(src));
}

/*-----------------------------------------*/

void
StdSolver::SAdd(double s1, Vector& dst, double s2, const Vector& src) const
{
  GetGV(dst).sadd(s1, s2, GetGV(src));
}

/*-----------------------------------------*/

double
StdSolver::Norm(const Vector& dst) const
{
  return GetGV(dst).norm();
}

/*-----------------------------------------*/

double
StdSolver::ScalarProduct(const Vector& y, const Vector& x) const
{
  return GetGV(y) * GetGV(x);
}

/*---------------------------------------------------*/

void
StdSolver::AssembleDualMatrix(Matrix& A, Vector& gu, double d)
{
  GlobalTimer.start("---> matrix");

  MatrixInterface& M = GetMatrix(A);

  HNAverage(gu);

  M.zero();
  GetDiscretization()->Matrix(
    M, GetGV(gu), *GetProblemDescriptor()->GetEquation(), d);
  M.transpose();

  // PeriodicMatrix() hier nicht getestet!
  PeriodicMatrix(A);
  DirichletMatrixOnlyRow(A);
  HNZero(gu);

  GlobalTimer.stop("---> matrix");
}

/*---------------------------------------------------*/

void
StdSolver::RhsCurve(Vector& f, const Curve& C, int comp, int N) const
{
  HNAverageData();

  GetDiscretization()->RhsCurve(GetGV(f), C, comp, N);

  HNZeroData();
  HNDistribute(f);
}
/*--------------------------------------------------------*/

double
StdSolver::ScalarProductWithFluctuations(DoubleVector& eta,
                                         const Vector& gf,
                                         const Vector& gz) const
{
  const GlobalVector& f = GetGV(gf);
  const GlobalVector& z = GetGV(gz);

  GlobalVector dz(f.ncomp(), f.n());

  dz.zero();
  Pi pi;
  pi.Init(GetMesh());
  pi.vmult(dz, z);

  for (int i = 0; i < z.n(); i++) {
    for (int c = 0; c < z.ncomp(); c++) {
      eta[i] += fabs(f(i, c) * dz(i, c));
    }
  }
  return dz * f;
}
} // namespace Gascoigne
