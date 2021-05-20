#include "solvers.h"
#include <algorithm>

#include "alediscretization.h"
#include "cuthillmckee.h"
#include "fmatrixblock.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "mginterpolatornested.h"
#include "sparse_umf.h"
#include "umfilu.h"

using namespace std;
using namespace Eigen;

extern double __DT, __THETA;

namespace Gascoigne {
template<int DIM>
FSISolver<DIM>::FSISolver()
  : StdSolver()
{}

template<int DIM>
void
FSISolver<DIM>::Form(Vector& y, const Vector& x, double d) const
{
  StdSolver::Form(y, x, d);
}

template<int DIM>
void
FSISolver<DIM>::ReInitExtensionMatrix()
{
  cout << "Initialize extension matrix. Rand noch nicht sinnvoll" << endl;

  // set up structure
  SparseStructure SA;
  GetDiscretization()->Structure(&SA);
  _LAP_A.ReInit(&SA);
  _LAP_M.SetMatrix(&_LAP_A);
  _LAP_M.ReInit(&SA);

  // ILU for partitioning
  if (!_directsolver) {
    //	_ILU_F.SetMatrix(GetMatrix());
    _ILU_S.SetMatrix(GetMatrix());

    _ILU_F.ReInit(&SA);
    _ILU_S.ReInit(&SA);
  }

  _LAP_A.zero();
  GlobalVector dummy(DIM + 1, GetMesh()->nnodes());

  // discretization for laplace, 1component, delete interface-test functions
  DiscretizationInterface* Q;
  vector<int> delf, dels;
  for (int i = 0; i < DIM; ++i)
    delf.push_back(i);
  if (DIM == 2) {
    Q = new AleQ2Lps2d;
    dynamic_cast<AleQ2Lps2d*>(Q)->InitInterfaceComponents(delf, dels);
  } else if (DIM == 3) {
    Q = new AleQ2Lps3d;
    dynamic_cast<AleQ2Lps3d*>(Q)->InitInterfaceComponents(delf, dels);
  } else
    abort();
  assert(Q);
  Q->BasicInit(_paramfile);
  Q->ReInit(GetMesh());
  ExtensionEquation<DIM> exteq;
  Q->Matrix(_LAP_A, dummy, exteq, 1.0);

  // Read dirichlet boundary data for extension problem
  IntSet dirichletextension;
  std::map<int, IntVector> dirichletextensioncomp;
  DataFormatHandler DF;
  DF.insert("dirichletextension", &dirichletextension);
  DF.insert("dirichletextensioncomp", &dirichletextensioncomp);
  FileScanner FS(DF, _paramfile, "BoundaryManager");

  for (auto color : dirichletextension) {
    auto extensioncomp = dirichletextensioncomp.find(color);
    if (extensioncomp == dirichletextensioncomp.end())
      continue;
    const auto& comps = extensioncomp->second;
    Q->StrongDirichletMatrix(_LAP_A, color, comps);
  }

  // modify in solid
  const vector<int>& s_nodes = GetAleDiscretization()->GetSolidL2G();
  vector<int> cv;
  for (int i = 0; i < DIM + 1; ++i)
    cv.push_back(i);

  // for (auto node : s_nodes)
  //_LAP_A.dirichlet_only_row(node,cv);
  for (vector<int>::const_iterator it = s_nodes.begin(); it != s_nodes.end();
       ++it)
    _LAP_A.dirichlet_only_row(*it, cv);
  cv.clear();
  cv.push_back(DIM);
  for (int i = 0; i < GetMesh()->nnodes(); ++i)
    _LAP_A.dirichlet(i, cv);

  // copy to mumps
  _LAP_M.ConstructStructure(vector<int>(0), _LAP_A);
  _LAP_M.copy_entries(&_LAP_A);
  _LAP_M.compute_ilu();
}

template<int DIM>
void
FSISolver<DIM>::ReInitMatrix()
{
  ReInitExtensionMatrix();

  GetDiscretization()->InitFilter(GetPfilter());
  SparseStructure SA;
  GetDiscretization()->Structure(&SA);
  _P.ReInit(&SA);
  _P.zero();

  StdSolver::ReInitMatrix();
}

template<int DIM>
void
FSISolver<DIM>::ComputeSAI() const
{
  if (_directsolver)
    return;

  _P.zero();
  StopWatch SW;
  SW.start();

  const SparseBlockMatrix<FMatrixBlock<DIM + 1>>& A =
    dynamic_cast<const SparseBlockMatrix<FMatrixBlock<DIM + 1>>&>(*GetMatrix());
  const ColumnDiagStencil& ST =
    dynamic_cast<const ColumnDiagStencil&>(*A.GetStencil());
  // const vector<FMatrixBlock<y> >& value = A.mat();

  int N = ST.n();

  //#pragma omp parallel for

  for (int k = 0; k < ST.n(); ++k) {
    // Spalte k von P
    // nicht-nullen in Spalte k genau bei nicht-Nullen in Zeile k
    vector<int> indexcol;
    for (int i = ST.start(k); i < ST.stop(k); ++i)
      indexcol.push_back(ST.col(i));
    set<int> indexrowset;
    for (int j = 0; j < indexcol.size(); ++j) {
      int ri = indexcol[j];
      for (int i = ST.start(ri); i < ST.stop(ri); ++i)
        indexrowset.insert(ST.col(i));
    }
    vector<int> indexrow;
    for (set<int>::const_iterator it = indexrowset.begin();
         it != indexrowset.end();
         ++it) {
      indexrow.push_back(*it);
    }

    // gleiche daten zum Suchen
    map<int, int> indexcolset;
    for (int i = 0; i < indexcol.size(); ++i)
      indexcolset[indexcol[i]] = i;

    const vector<FMatrixBlock<DIM + 1>>& value = A.mat();

    // Initialisierung Matrix
    Matrix<double, Dynamic, Dynamic> Aklein(indexrow.size() * (DIM + 1),
                                            indexcol.size() * (DIM + 1));
    Aklein.setZero();

    // zeilen durchlaufen
    for (int rAklein = 0; rAklein < indexrow.size();
         ++rAklein) // Initialisierung Aklein
    {
      int rA = indexrow[rAklein];
      for (int pA = ST.start(rA); pA < ST.stop(rA); ++pA) {
        int cA = ST.col(pA);
        if (indexcolset.find(cA) != indexcolset.end()) {
          for (int iF = rAklein * (DIM + 1); iF < (rAklein + 1) * (DIM + 1);
               ++iF) {
            int cP = indexcolset[cA];
            for (int jF = cP * (DIM + 1); jF < (cP + 1) * (DIM + 1); ++jF)
              Aklein(iF, jF) =
                value[pA](iF - rAklein * (DIM + 1), jF - cP * (DIM + 1));
          }
        }
      }
    }
    // QR-Zerlegung berechnen

    ColPivHouseholderQR<Matrix<double, Dynamic, Dynamic>> qr(Aklein.rows(),
                                                             Aklein.cols());
    qr.compute(Aklein);

    // mk berechnen
    for (int mS = 0; mS < (DIM + 1); ++mS) // ueber die Spalten von mk iterieren
    {
      Matrix<double, Dynamic, 1> ek(indexrow.size() * (DIM + 1));
      for (int rAklein = 0; rAklein < indexrow.size();
           ++rAklein) // Initialisierung Aklein
      {
        int rA = indexrow[rAklein];
        if (rA == k) // Initialisierung rechte Seite ebenfalls spaltenweise
        {
          for (int iE = (rAklein * (DIM + 1)); iE < (rAklein + 1) * (DIM + 1);
               ++iE) {
            if (iE - (rAklein * (DIM + 1)) == mS) // ?????
              ek(iE) = 1;
            else
              ek(iE) = 0;
          }
        } else {
          for (int iE = rAklein * (DIM + 1); iE < (rAklein + 1) * (DIM + 1);
               ++iE) {
            ek(iE) = 0;
          }
        }
      }
      Matrix<double, Dynamic, 1> mk = qr.solve(ek);
      // cout << Aklein << "\n" << endl;
      // cout << mk << "\n" << endl;
      for (int cPP = 0; cPP < indexcol.size(); ++cPP) {
        int rP = indexcol[cPP];
        // P[rP,k] = mk(cPP,0);
        for (int iP = 0; iP < (DIM + 1); ++iP) {
          if (fabs(mk((cPP * (DIM + 1)) + iP) < 1.e-14)) {
            _P.mat(ST.Find(rP, k))->value(iP, mS) = 0;
          } else
            _P.mat(ST.Find(rP, k))->value(iP, mS) = mk((cPP * (DIM + 1)) + iP);
        }
      }
    }
  }

  SW.stop();
  cout << "Precond: " << N << "\t" << SW.read() << endl;
}

template<int DIM>
void
FSISolver<DIM>::ComputeIlu(const Vector& gu) const
{
  ComputeSAI();

  if ((_matrixtype == "sparseumf") && (!_directsolver)) {
    IntVector perm;
    GetIlu()->ConstructStructure(perm, *GetMatrix());
    GetIlu()->copy_entries(GetMatrix());
    GetIlu()->compute_ilu();
  } else

    // if (!_directsolver)
    // {
    //   int ncomp = GetProblemDescriptor()->GetEquation()->GetNcomp();
    //   PermutateIlu(gu);

    //   GetIlu()->zero();
    //   GetIlu()->copy_entries(GetMatrix());

    //   modify_ilu(*GetIlu(), ncomp);
    //   GetIlu()->compute_ilu();

    //   // FOR FSI ILU
    //   IntVector perm(GetMatrix()->GetStencil()->n());
    //   iota(perm.begin(),perm.end(),0);
    //   CuthillMcKee    cmc(GetMatrix()->GetStencil());
    //   cmc.Permutate      (perm);

    //   _ILU_F.ConstructStructure(perm,
    // 				*GetMatrix());
    //   _ILU_S.ConstructStructure(vector<int>(0),
    // 				*GetMatrix());
    //   HASHSET<int> nullset;
    //   _ILU_F.copy_entries(GetAleDiscretization()->GetFluidG2L(),
    //   GetAleDiscretization()->GetInterfaceNodes(), GetMatrix());
    //   //      _ILU_F.copy_entries(GetAleDiscretization()->GetFluidG2L(),
    //   nullset, GetMatrix());
    //   _ILU_S.copy_entries(GetAleDiscretization()->GetSolidG2L(), nullset,
    //   GetMatrix()); _ILU_F.compute_ilu(); _ILU_S.compute_ilu();
    // }
    // else
    StdSolver::ComputeIlu(gu);
}

template<int DIM>
void
FSISolver<DIM>::modify_ilu(IluInterface& I, int ncomp) const
{
  if (GetSolverData().GetIluModify().size() == 0)
    return;
  if (GetSolverData().GetIluModify().size() != ncomp) {
    cerr << "ERROR: GetSolverData().GetIluModify().size()="
         << GetSolverData().GetIluModify().size() << " and ";
    cerr << "ncomp=" << ncomp << endl;
    abort();
  }

  for (int c = 0; c < ncomp; c++) {
    double s = GetSolverData().GetIluModify(c);
    I.modify(c, s);
  }
}

template<int DIM>
void
FSISolver<DIM>::smooth(int niter, Vector& x, const Vector& y, Vector& h) const
{
  // const std::vector<int>& fl2g = GetAleDiscretization()->GetFluidL2G();
  // const std::vector<int>& sl2g = GetAleDiscretization()->GetSolidL2G();
  // const HASHSET<int>& interface =
  // GetAleDiscretization()->GetInterfaceNodes(); if ((_matrixtype ==
  // "splitumf") && (!_directsolver))
  // {
  //   abort();
  // }
  // else if (GetSolverData().GetLinearSmooth() == "ilu")
  // {
  //   assert(!_directsolver);
  //   double omega = GetSolverData().GetOmega();
  //   for (int iter = 0; iter < niter; iter++)
  //   {
  // 	// FLUID PROBLEM
  //     MatrixResidual(h, x, y);
  // 	// Set Dirichlet-Zero within the Solid domain for velocity & Pressure
  // 	// on the interface: only velocity to zero
  // 	for (auto it : sl2g)
  // 	  if (interface.find(it)==interface.end())
  // 	    GetGV(h).zero_node(it);
  // 	for (auto it : interface)
  // 	  for (int c=0;c<DIM;++c)
  // 	    GetGV(h)(it,c+1)=0.0;
  // 	_ILU_F.solve(GetGV(h));
  // 	Add(x, omega, h);

  // 	// SOLID PROBLEM
  //     MatrixResidual(h, x, y);
  // 	for (auto it : fl2g)
  // 	  if (interface.find(it)==interface.end())
  // 	    GetGV(h).zero_node(it);
  // 	_ILU_S.solve(GetGV(h));
  // 	for (auto it : fl2g)
  // 	  if (interface.find(it)==interface.end())
  // 	    GetGV(h).zero_node(it);

  // 	Add(x, omega, h);
  //   }
  // }
  // else

  if (GetSolverData().GetLinearSmooth() == "sai") {
    double omega = GetSolverData().GetOmega();

    for (int iter = 0; iter < niter; iter++) {
      MatrixResidual(h, x, y); // h = y - A x
      //  	    GetIlu()->solve(GetGV(h)); // h = (LR)^{-1} h
      //  	    Add(x,omega,h);            // x = x + omega * h
      _P.vmult(GetGV(x), GetGV(h), omega);
      SubtractMean(x);
    }
  } else
    StdSolver::smooth(niter, x, y, h);
}

template<int DIM>
void
FSISolver<DIM>::smooth_exact(Vector& x, const Vector& y, Vector& help) const
{
  StdSolver::smooth_exact(x, y, help);
}

/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

template<int DIM>
MatrixInterface*
FSISolver<DIM>::NewMatrix(int ncomp, const std::string& matrixtype)
{

#ifdef __WITH_UMFPACK__
  if (_directsolver && _useUMFPACK)
    return StdSolver::NewMatrix(ncomp, matrixtype);
#endif

  if (matrixtype == "sparseumf")
    return StdSolver::NewMatrix(ncomp, "block");
  else
    return StdSolver::NewMatrix(ncomp, matrixtype);
}

template<int DIM>
IluInterface*
FSISolver<DIM>::NewIlu(int ncomp, const string& matrixtype)
{
#ifdef __WITH_UMFPACK__
  // if(_directsolver && _useUMFPACK)             return new
  // UmfIluLong(GetMatrix());
  if (_directsolver && _useUMFPACK)
    return new UmfIlu(GetMatrix());
//  if(_directsolver && _useUMFPACK) {cout<<"keine ahnung was hier
//  passiert"<<endl; abort();}
#endif
  if (matrixtype == "block") {
    if (ncomp == 1)
      return new SparseUmf<FMatrixBlock<1>>(GetMatrix());
    else if (ncomp == 2)
      return new SparseUmf<FMatrixBlock<2>>(GetMatrix());
    else if (ncomp == 3)
      return new SparseUmf<FMatrixBlock<3>>(GetMatrix());
    else if (ncomp == 4)
      return new SparseUmf<FMatrixBlock<4>>(GetMatrix());
    else
      abort();
  } else
    return StdSolver::NewIlu(ncomp, matrixtype);

  // SplittingIlu<DIM>* TMP = new SplittingIlu<DIM>;
  // TMP->SetInterface(GetAleDiscretization()->GetFluidL2G(),
  // 		      GetAleDiscretization()->GetSolidL2G(),
  // 		      GetAleDiscretization()->GetFluidG2L(),
  // 		      GetAleDiscretization()->GetSolidG2L(),
  // 		      GetAleDiscretization()->GetInterfaceNodes());

  // return TMP;
}

/*-------------------------------------------------------*/

template<int DIM>
void
FSISolver<DIM>::DeleteSolidPressure(Vector& gf) const
{
  ////////////////////

  const HASHSET<int>& i_nodes = GetAleDiscretization()->GetInterfaceNodes();
  const vector<int>& s_nodes = GetAleDiscretization()->GetSolidL2G();
  /*
    for (auto node : s_nodes)
    if (i_nodes.find(node)==i_nodes.end())
    {
    GetGV(gf)(node,0) = 0.0;
    }
  */
  for (vector<int>::const_iterator it = s_nodes.begin(); it != s_nodes.end();
       ++it)
    if (i_nodes.find(*it) == i_nodes.end()) {
      GetGV(gf)(*it, 0) = 0.0;
    }
  // for (auto node : s_nodes)
  //   for (int c=0;c<DIM;++c)
  // 	GetGV(gf)(node,c+1+DIM) = 0.0;
}

template<int DIM>
void
FSISolver<DIM>::SetBoundaryVectorZero(Vector& gf) const
{
  StdSolver::SetBoundaryVectorZero(gf);

  DeleteSolidPressure(gf);
}

template<int DIM>
void
FSISolver<DIM>::SetBoundaryVector(Vector& gf) const
{
  StdSolver::SetBoundaryVector(gf);

  DeleteSolidPressure(gf);
}

template<int DIM>
void
FSISolver<DIM>::AssembleMatrix(const Vector& gu, double d)
{

  StdSolver::AssembleMatrix(gu, d);

  if ((_directsolver) || (_matrixtype == "block") ||
      (_matrixtype == "sparseumf")) {

    // Modify for pressure zero in Solid-Part
    const HASHSET<int>& i_nodes = GetAleDiscretization()->GetInterfaceNodes();
    const vector<int>& s_nodes = GetAleDiscretization()->GetSolidL2G();
    vector<int> cv;
    cv.push_back(0);
    for (int i = 0; i < s_nodes.size(); ++i)
      if (i_nodes.find(s_nodes[i]) == i_nodes.end())
        GetMatrix()->dirichlet(s_nodes[i], cv);
  }
}

template<int DIM>
DiscretizationInterface*
FSISolver<DIM>::NewDiscretization(int dimension, const string& discname)
{
  if (dimension == 2) {
    if (discname == "AleQ1") {
      abort();
      return new AleQ12d;
    }
    //	else if (discname=="AleQ2")               return new AleQ22d;
    else if (discname == "AleQ1Lps") {
      abort();
      return new AleQ1Lps2d;
    }

    else if (discname == "AleQ2Lps") {
      DiscretizationInterface* X = new AleQ2Lps2d;
      vector<int> delf, dels;
      dels.push_back(0);
      dynamic_cast<AleQ2Lps2d*>(X)->InitInterfaceComponents(delf, dels);
      return X;
    }

    else
      return StdSolver::NewDiscretization(dimension, discname);
  } else if (dimension == 3) {
    // if      (discname=="AleQ1")               return new AleQ13d;
    if (discname == "AleQ1Lps")
      return new AleQ1Lps3d;
    else if (discname == "AleQ2Lps") {
      DiscretizationInterface* X = new AleQ2Lps3d;
      vector<int> delf, dels;
      //	    for (int i=0;i<DIM;++i) delf.push_back(i+1+DIM);
      dels.push_back(0);
      dynamic_cast<AleQ2Lps3d*>(X)->InitInterfaceComponents(delf, dels);
      return X;
    }

    // else if (discname=="AleQ2Lps")            return new AleQ2Lps3d;
    // else
    return StdSolver::NewDiscretization(dimension, discname);
  } else
    abort();
}

template<int DIM>
void
FSISolver<DIM>::reinit_interface_element(
  int en,
  const nvector<int>& indices,
  HASHMAP<int, std::vector<int>>& solid_interface_cells,
  HASHMAP<int, std::vector<int>>& fluid_interface_cells,
  HASHSET<int>& interface_nodes,
  int material)
{
  vector<int> ni;
  for (int i = 0; i < indices.size(); ++i) {
    if (interface_nodes.find(indices[i]) != interface_nodes.end()) {
      ni.push_back(i);
    }
  }

  if (ni.size() > 0) {
    if (material == 1) { // solid cell
      solid_interface_cells[en] = ni;
    } else if (material == 2) { // fluid cell
      fluid_interface_cells[en] = ni;
    } else {
      cout << "	Fluid Cells need to have the material value 2 and solid "
              "cells the material value 1!"
           << endl;
      abort();
    }
  }
}

template<int DIM>
void
FSISolver<DIM>::reinit_element(int en,
                               const nvector<int>& indices,
                               HASHSET<int>& fluid_cells,
                               HASHSET<int>& solid_cells,
                               set<int>& fluid_nodes,
                               set<int>& solid_nodes,
                               int material)
{

  if (material == 1) { // solid cell
    for (int i = 0; i < indices.size(); ++i) {
      solid_nodes.insert(indices[i]);
      solid_cells.insert(en);
    }
  } else if (material == 2) { // fluid cell
    for (int i = 0; i < indices.size(); ++i) {
      fluid_nodes.insert(indices[i]);
      fluid_cells.insert(en);
    }
  } else {
    cout << "	Fluid Cells need to have the material value 2 and solid "
            "cells the material value 1!"
         << endl;
    abort();
  }
}

template<int DIM>
void
FSISolver<DIM>::ReInitInterface(AleBaseDiscretization* ALEDISC)
{
  HASHMAP<int, std::vector<int>>& solid_interface_cells =
    ALEDISC->GetSolidInterfaceCells();
  HASHMAP<int, std::vector<int>>& fluid_interface_cells =
    ALEDISC->GetFluidInterfaceCells();
  HASHSET<int>& interface_nodes = ALEDISC->GetInterfaceNodes();
  HASHSET<int>& fluid_cells = ALEDISC->GetFluidCells();
  HASHSET<int>& solid_cells = ALEDISC->GetSolidCells();
  vector<int>& fluid_l2g = ALEDISC->GetFluidL2G();
  vector<int>& solid_l2g = ALEDISC->GetSolidL2G();
  HASHMAP<int, int>& fluid_g2l = ALEDISC->GetFluidG2L();
  HASHMAP<int, int>& solid_g2l = ALEDISC->GetSolidG2L();

  set<int> fluid_nodes, solid_nodes;

  solid_interface_cells.clear();
  fluid_interface_cells.clear();
  interface_nodes.clear();
  fluid_cells.clear();
  solid_cells.clear();

  int dim = GetMesh()->dimension();

  if (dim == 2) {
    const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*>(GetMesh());
    assert(M);
    if ((GetDiscretization()->GetName() == "Q1 Ale 2d Lps") ||
        (GetDiscretization()->GetName() == "Q1 Ale 2d")) {
      // Einsortieren der Solid und Fluid Nodes
      for (int c = 0; c < M->ncells(); ++c) {
        reinit_element(c,
                       M->IndicesOfCell(c),
                       fluid_cells,
                       solid_cells,
                       fluid_nodes,
                       solid_nodes,
                       M->material(c));
      }
      // Interface Node: Sowohl Fluid als auch Solid Node
      // Kann erst aufgerufen werden wenn man durch alle Zellen einmal durch
      // ist!
      for (set<int>::const_iterator it = fluid_nodes.begin();
           it != fluid_nodes.end();
           ++it)
        if (solid_nodes.find(*it) != solid_nodes.end())
          interface_nodes.insert(*it);
      // Interface Cells und Interfacenodes on InterfaceCells abspeichern
      for (int c = 0; c < M->ncells(); ++c) {
        reinit_interface_element(c,
                                 M->IndicesOfCell(c),
                                 solid_interface_cells,
                                 fluid_interface_cells,
                                 interface_nodes,
                                 M->material(c));
      }
    } else if ((GetDiscretization()->GetName() == "Q2 Ale 2d Lps") ||
               (GetDiscretization()->GetName() == "Q2 Ale 2d")) {
      // Einsortieren der Solid und Fluid Nodes
      for (int c = 0; c < M->npatches(); ++c) {
        reinit_element(c,
                       *(M->IndicesOfPatch(c)),
                       fluid_cells,
                       solid_cells,
                       fluid_nodes,
                       solid_nodes,
                       M->material_patch(c));
      }
      // Interface Node: Sowohl Fluid als auch Solid Node
      // Kann erst aufgerufen werden wenn man durch alle Zellen einmal durch
      // ist!
      for (set<int>::const_iterator it = fluid_nodes.begin();
           it != fluid_nodes.end();
           ++it)
        if (solid_nodes.find(*it) != solid_nodes.end())
          interface_nodes.insert(*it);
      // Interface Cells und Interfacenodes on InterfaceCells abspeichern
      for (int c = 0; c < M->npatches(); ++c) {
        reinit_interface_element(c,
                                 *(M->IndicesOfPatch(c)),
                                 solid_interface_cells,
                                 fluid_interface_cells,
                                 interface_nodes,
                                 M->material_patch(c));
      }
    } else
      abort();
  } else if (dim == 3) {
    const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*>(GetMesh());
    assert(M);

    if (GetDiscretization()->GetName() == "Q1 Ale 3d Lps") {
      // Einsortieren der Solid und Fluid Nodes
      for (int c = 0; c < M->ncells(); ++c) {
        reinit_element(c,
                       M->IndicesOfCell(c),
                       fluid_cells,
                       solid_cells,
                       fluid_nodes,
                       solid_nodes,
                       M->material(c));
      }
      // Interface Node: Sowohl Fluid als auch Solid Node
      // Kann erst aufgerufen werden wenn man durch alle Zellen einmal durch
      // ist!
      for (set<int>::const_iterator it = fluid_nodes.begin();
           it != fluid_nodes.end();
           ++it)
        if (solid_nodes.find(*it) != solid_nodes.end())
          interface_nodes.insert(*it);
      // Interface Cells und Interfacenodes on InterfaceCells abspeichern
      for (int c = 0; c < M->ncells(); ++c) {
        reinit_interface_element(c,
                                 M->IndicesOfCell(c),
                                 solid_interface_cells,
                                 fluid_interface_cells,
                                 interface_nodes,
                                 M->material(c));
      }
    } else if (GetDiscretization()->GetName() == "Q2 Ale 3d Lps") {
      // Einsortieren der Solid und Fluid Nodes
      for (int c = 0; c < M->npatches(); ++c) {
        reinit_element(c,
                       *(M->IndicesOfPatch(c)),
                       fluid_cells,
                       solid_cells,
                       fluid_nodes,
                       solid_nodes,
                       M->material_patch(c));
      }
      // Interface Node: Sowohl Fluid als auch Solid Node
      // Kann erst aufgerufen werden wenn man durch alle Zellen einmal durch
      // ist!
      for (set<int>::const_iterator it = fluid_nodes.begin();
           it != fluid_nodes.end();
           ++it)
        if (solid_nodes.find(*it) != solid_nodes.end())
          interface_nodes.insert(*it);
      // Interface Cells und Interfacenodes on InterfaceCells abspeichern
      for (int c = 0; c < M->npatches(); ++c) {
        reinit_interface_element(c,
                                 *(M->IndicesOfPatch(c)),
                                 solid_interface_cells,
                                 fluid_interface_cells,
                                 interface_nodes,
                                 M->material_patch(c));
      }
    } else {
      std::cout << GetDiscretization()->GetName() << std::endl;

      abort();
    }
  } else
    abort();

  // Nodes Fluid & Solid,  local <-> global (fluid nodes include interface and
  // same for solid)
  fluid_l2g.clear();
  solid_l2g.clear();
  fluid_g2l.clear();
  solid_g2l.clear();

  // l2g
  for (set<int>::const_iterator it = fluid_nodes.begin();
       it != fluid_nodes.end();
       ++it)
    fluid_l2g.push_back(*it);
  for (set<int>::const_iterator it = solid_nodes.begin();
       it != solid_nodes.end();
       ++it)
    solid_l2g.push_back(*it);

  // g2l
  for (int i = 0; i < fluid_l2g.size(); ++i)
    fluid_g2l[fluid_l2g[i]] = i;
  for (int i = 0; i < solid_l2g.size(); ++i)
    solid_g2l[solid_l2g[i]] = i;
}

// --------------------------------------------------

template<int DIM>
void
FSISolver<DIM>::NewMesh(const MeshInterface* mp)
{
  StdSolver::NewMesh(mp);

  ReInitInterface(GetAleDiscretization());
}

template<>
void
FSISolver<2>::PointVisu(const string& name,
                        const GlobalVector& u,
                        int iter) const
{
  Vector def("def");
  const GlobalVector& DEF = GetGV(def);

  GlobalVector U;
  U.ncomp() = 2 * u.ncomp();
  assert(2 * 2 + 2 == U.ncomp());
  U.resize(u.n());

  const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*>(GetMesh());
  assert(M);

  for (int i = 0; i < u.n(); ++i) {
    const Vertex2d v = M->vertex2d(i);

    int domain;

    if (find(GetAleDiscretization()->GetFluidL2G().begin(),
             GetAleDiscretization()->GetFluidL2G().end(),
             i) != GetAleDiscretization()->GetFluidL2G().end()) {
      if (GetAleDiscretization()->GetInterfaceNodes().find(i) !=
          GetAleDiscretization()->GetInterfaceNodes().end())
        domain = 0;
      else
        domain = -1;
    } else {
      if (find(GetAleDiscretization()->GetSolidL2G().begin(),
               GetAleDiscretization()->GetSolidL2G().end(),
               i) != GetAleDiscretization()->GetSolidL2G().end()) {
        if (GetAleDiscretization()->GetInterfaceNodes().find(i) !=
            GetAleDiscretization()->GetInterfaceNodes().end())
          domain = 0;
        else
          domain = 1;
      } else {
        cout << "error writing mesh. Node neither fluid nor solid" << endl;
        abort();
      }
    }

    for (int c = 0; c < u.ncomp(); ++c)
      U(i, c) = u(i, c);
    for (int c = 0; c < 2; ++c)
      U(i, c + 2 + 1) = DEF(i, c);
    U(i, U.ncomp() - 1) = domain;
  }
  StdSolver::PointVisu(name, U, iter);
}
template<>
void
FSISolver<3>::PointVisu(const string& name,
                        const GlobalVector& u,
                        int iter) const
{
  Vector def("def");
  const GlobalVector& DEF = GetGV(def);
  GlobalVector U;
  U.ncomp() = 2 * u.ncomp();
  // 1 pressure ,3 fluid ,3 solid ,1 domain
  assert(2 * 3 + 2 == U.ncomp());
  U.resize(u.n());

  const GascoigneMesh3d* M = dynamic_cast<const GascoigneMesh3d*>(GetMesh());
  assert(M);

  for (int i = 0; i < u.n(); ++i) {
    const Vertex3d v = M->vertex3d(i);
    int domain;
    if (find(GetAleDiscretization()->GetFluidL2G().begin(),
             GetAleDiscretization()->GetFluidL2G().end(),
             i) != GetAleDiscretization()->GetFluidL2G().end()) {
      if (GetAleDiscretization()->GetInterfaceNodes().find(i) !=
          GetAleDiscretization()->GetInterfaceNodes().end())
        domain = 0;
      else
        domain = -1;
    } else {
      if (find(GetAleDiscretization()->GetSolidL2G().begin(),
               GetAleDiscretization()->GetSolidL2G().end(),
               i) != GetAleDiscretization()->GetSolidL2G().end()) {
        if (GetAleDiscretization()->GetInterfaceNodes().find(i) !=
            GetAleDiscretization()->GetInterfaceNodes().end())
          domain = 0;
        else
          domain = 1;
      } else {
        cout << "error writing mesh. Node neither fluid nor solid" << endl;
        abort();
      }
    }

    for (int c = 0; c < u.ncomp(); ++c)
      U(i, c) = u(i, c);
    for (int c = 0; c < 3; ++c)
      U(i, c + 3 + 1) = DEF(i, c);
    U(i, U.ncomp() - 1) = domain;
  }
  StdSolver::PointVisu(name, U, iter);
}
//////////////////////////////////////////////////

//////// deformation update
template<int DIM>
void
FSISolver<DIM>::SolveExtension(Vector& x)
{

  Vector def("def");
  GlobalVector& DEF = GetGV(def);
  for (int c = 0; c < DIM; ++c) {
    //	DEF.CompEq(c,1.0,c+DIM+1,X);
    const vector<int>& f_nodes = GetAleDiscretization()->GetFluidL2G();
    const HASHSET<int>& i_nodes = GetAleDiscretization()->GetInterfaceNodes();
    // for (auto node : f_nodes)
    //  if (i_nodes.find(node)==i_nodes.end())
    //    DEF(node,c)=0.0;
    for (vector<int>::const_iterator it = f_nodes.begin(); it != f_nodes.end();
         ++it)
      if (i_nodes.find(*it) == i_nodes.end())
        DEF(*it, c) = 0.0;
  }
  DEF.zero_comp(DIM);

  _LAP_M.solve(DEF);
}

template<int DIM>
void
FSIMultiLevelSolver<DIM>::UpdateDeformation(Vector& x)
{

  // update solid deformation
  GlobalVector& X = GetSolver()->GetGV(x);

  Vector def("def");
  GlobalVector& DEF = GetSolver()->GetGV(def);
  Vector defold("defold");
  GlobalVector& DEFOLD = GetSolver()->GetGV(defold);
  Vector old("old");
  const GlobalVector& OLD = GetSolver()->GetGV(old);

  const vector<int>& s_nodes =
    GetFSISolver()->GetAleDiscretization()->GetSolidL2G();
  /*for (auto node : s_nodes)
    {
    for (int c=0;c<DIM;++c)
    {
    DEF(node,c) = DEFOLD(node,c)
    + __DT * (1.0-__THETA) * OLD(node,c+1)
    + __DT * __THETA       * X(node,c+1);
    }
    }
  */
  for (vector<int>::const_iterator it = s_nodes.begin(); it != s_nodes.end();
       ++it) {
    for (int c = 0; c < DIM; ++c) {
      DEF(*it, c) = DEFOLD(*it, c) + __DT * (1.0 - __THETA) * OLD(*it, c + 1) +
                    __DT * __THETA * X(*it, c + 1);
    }
  }

  // update deformation in fluid
  GetFSISolver()->SolveExtension(x);
}

template<int DIM>
double
FSIMultiLevelSolver<DIM>::NewtonUpdate(double& rr,
                                       Vector& x,
                                       Vector& dx,
                                       Vector& r,
                                       const Vector& f,
                                       NLInfo& nlinfo)
{

  const CGInfo& linfo = nlinfo.GetLinearInfo();
  bool lex = linfo.control().status() == "exploded";

  double nn = NewtonNorm(dx);
  double nr = GetSolver(ComputeLevel)->Norm(r);

  if (nn > 1.e30)
    lex = 1;
  if (!(nn >= 0.))
    lex = 1;
  if (nr > 1.e30)
    lex = 1;
  if (!(nr >= 0.))
    lex = 1;

  if (lex) {
    nlinfo.control().status() = "diverged";
    cerr << "linear : " << linfo.control().status() << endl;
    cerr << "nonlinear : " << nn << endl;
    return NewtonNorm(dx);
  }

  double omega = 0.7;
  double relax = 1.;

  GetSolver(ComputeLevel)->SetPeriodicVectorZero(dx);

  Vector XXX("old");
  GlobalVector OLD = GetSolver()->GetGV(XXX);
  GetSolver(ComputeLevel)->Add(x, relax, dx);
  UpdateDeformation(x);

  NewtonResidual(r, x, f);
  rr = NewtonNorm(r);

  string message = "";
  for (int iter = 0; iter < nlinfo.user().maxrelax(); iter++) {
    message = nlinfo.check_damping(iter, rr);

    if (message == "ok")
      break;
    if (message == "continue") {
      GetSolver(ComputeLevel)->Add(x, relax * (omega - 1.), dx);
      UpdateDeformation(x);

      NewtonResidual(r, x, f);
      rr = NewtonNorm(r);
      relax *= omega;
      continue;
    }
    if (message == "exploded") {
      GetSolver(ComputeLevel)->Add(x, -relax, dx);
      UpdateDeformation(x);

      relax = 0.;
      cout << "Damping exploded !!!!!" << endl;
      nlinfo.control().status() = "diverged";
      break;
    }
  }
  return NewtonNorm(dx);
}

template class FSISolver<2>;
template class FSISolver<3>;
template class FSIMultiLevelSolver<2>;
template class FSIMultiLevelSolver<3>;
} // namespace Gascoigne
