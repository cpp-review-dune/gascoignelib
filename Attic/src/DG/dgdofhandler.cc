#include "dgdofhandler.h"
#include "dgbase.h"
#include "gascoignehash.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "sparsestructure.h"

namespace Gascoigne {
struct TmpEdgeHash {
  std::size_t operator()(const std::array<size_t, 2> &a) const {
    std::size_t h = 0;
    for (auto e : a)
      h += std::hash<int>{}(e);
    return h;
  }
};

template <class BASE>
void DGDofHandler<BASE>::InitFromGascoigneMesh(const GascoigneMesh *M) {
  assert(M);
  const GascoigneMesh2d *M2 = dynamic_cast<const GascoigneMesh2d *>(M);
  const GascoigneMesh3d *M3 = dynamic_cast<const GascoigneMesh3d *>(M);
  if (!(M2 || M3)) {
    std::cerr << "cannot cast Mesh to GascoigneMesh2d/3d" << std::endl;
    abort();
  }
  ////////////////////////////////////////////////// ELEMENTS
  _ndofs = M->ncells() * BASE::N;
  _elements.resize(M->ncells());
  int gi = 0;
  for (int i = 0; i < M->ncells(); ++i)
    for (int d = 0; d < BASE::N; ++d, ++gi)
      _elements[i][d] = gi;
  assert(_ndofs == gi);

  ////////////////////////////////////////////////// EDGES
  assert(M2); // only 2d
  // tmp-structure to identify edges, point nodes + edgeno.
  HASHMAP<std::array<size_t, 2>, size_t, TmpEdgeHash> edgetmp;
  int ei = 0;
  int itoj[4] = {0, 1, 3, 2};
  for (int cell = 0; cell < M->ncells(); ++cell) {
    const IntVector ioc = M->IndicesOfCell(cell);
    for (int e = 0; e < 4; ++e) {
      std::array<size_t, 2> EE;
      EE[0] = std::min(ioc[itoj[e]], ioc[itoj[(e + 1) % 4]]);
      EE[1] = std::max(ioc[itoj[e]], ioc[itoj[(e + 1) % 4]]);

      if (edgetmp.find(EE) == edgetmp.end()) {
        edgetmp[EE] = ei;
        ++ei;
      }
    }
  }
  assert(ei == edgetmp.size());

  // setup edges
  EdgeType noedge;
  noedge.fill(-1);
  _edges.clear();
  _edges.resize(ei, noedge);

  for (int cell = 0; cell < M->ncells(); ++cell) {
    const IntVector ioc = M->IndicesOfCell(cell);
    for (int e = 0; e < 4; ++e) {
      std::array<size_t, 2> EE;
      EE[0] = std::min(ioc[itoj[e]], ioc[itoj[(e + 1) % 4]]);
      EE[1] = std::max(ioc[itoj[e]], ioc[itoj[(e + 1) % 4]]);
      auto IT = edgetmp.find(EE);
      assert(IT != edgetmp.end());
      int index = IT->second;

      assert(index < _edges.size());
      EdgeType &E = _edges[index];

      if (E[0] == -1) {
        E[0] = cell; // cell number
        E[2] = e;    // local index
      } else if (E[1] == -1) {
        E[1] = cell; // cell number
        E[3] = e;    // local index
      } else {
        std::cerr << "Construction of Edge failed " << E[0] << " " << E[1]
                  << " " << E[2] << " " << E[3] << "\t" << ioc << std::endl;
        abort();
      }
    }
  }

  std::cout << "DG Dofs with " << _elements.size() << " elements and "
            << _edges.size() << " edges" << std::endl;
}

//////////////////////////////////////////// Global/Local

template <class BASE>
void DGDofHandler<BASE>::GlobalToLocal(LocalVector &U, const GlobalVector &u,
                                       int iq) const {
  assert(iq < nelements());
  const ElementType &LI = getelement(iq);
  assert(BASE::N == LI.size());
  U.ReInit(u.ncomp(), BASE::N);

  for (int ii = 0; ii < LI.size(); ii++) {
    int i = LI[ii];
    U.equ_node(ii, i, u);
  }
}

template <class BASE>
void DGDofHandler<BASE>::LocalToGlobal(GlobalVector &f, const LocalVector &F,
                                       int iq, double s) const {
  assert(iq < nelements());
  const ElementType &LI = getelement(iq);
  assert(BASE::N == LI.size());

  assert(F.ncomp() == f.ncomp());
  assert(F.n() == BASE::N);

  for (int ii = 0; ii < LI.size(); ii++) {
    int i = LI[ii];
    f.add_node(i, s, ii, F);
  }
}

//////////////////////////////////////////////////

template <class BASE>
void DGDofHandler<BASE>::LocalToGlobalMatrix(MatrixInterface &A, EntryMatrix &E,
                                             int iq, double s) const {
  assert(iq < nelements());
  const ElementType &LI = getelement(iq);

  /// temp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IntVector indices;
  for (auto i : LI)
    indices.push_back(i);

  auto start = indices.begin();
  auto stop = indices.end();
  A.entry(start, stop, E, s);
}

//////////////////////////////////////////////////

template <class BASE>
void DGDofHandler<BASE>::LocalToGlobalMatrix(MatrixInterface &A, EntryMatrix &E,
                                             int q1, int q2, double s) const {
  assert(q1 < nelements());
  assert(q2 < nelements());
  const ElementType &LM = getelement(q1);
  const ElementType &LS = getelement(q2);

  /// temp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IntVector iq1;
  for (auto i : LM)
    iq1.push_back(i);
  IntVector iq2;
  for (auto i : LS)
    iq2.push_back(i);

  auto startq1 = iq1.begin();
  auto stopq1 = iq1.end();
  auto startq2 = iq2.begin();
  auto stopq2 = iq2.end();
  A.entry(startq1, stopq1, startq2, stopq2, E, s);
}

//// Matrix Structure
template <class BASE>
void DGDofHandler<BASE>::Structure(SparseStructureInterface *SI) const {
  SparseStructure *S = dynamic_cast<SparseStructure *>(SI);
  assert(S);

  S->build_begin(ndofs());
  for (auto e : _elements) {
    S->build_add(e.begin(), e.end());
  }
  for (auto e : _edges) {
    if (e[1] == -1)
      continue;
    std::vector<size_t> tmp;
    tmp.insert(tmp.begin(), _elements[e[0]].begin(), _elements[e[0]].end());
    tmp.insert(tmp.begin(), _elements[e[1]].begin(), _elements[e[1]].end());
    S->build_add(tmp.begin(), tmp.end());
  }

  S->build_end();
}

template class DGDofHandler<BASEQ12D>;
template class DGDofHandler<BASEQ22D>;
} // namespace Gascoigne
