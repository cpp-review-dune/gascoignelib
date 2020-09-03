/*----------------------------   alebasediscretization.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __alebasediscretization_H
#define __alebasediscretization_H
/*----------------------------   alebasediscretization.h     ---------------------------*/

#include "q12d.h"
#include "q1lps2d.h"
#include "q13d.h"
#include "q1lps3d.h"
#include "chi.h"
#include "q2lps2d.h"
#include "q2lps3d.h"

#include "gascoignemesh2d.h"

extern bool __ADJOINT;

namespace Gascoigne {

class AleBaseDiscretization {
public:
  std::vector<int> __DEL_F, __DEL_S;

  HASHMAP<int, std::vector<int>> __solid_interface_cells, __fluid_interface_cells;
  HASHSET<int>                   __fluid_cells, __solid_cells, __interface_nodes;

  std::vector<int>  __fluid_l2g, __solid_l2g;
  HASHMAP<int, int> __fluid_g2l, __solid_g2l;

public:
  HASHMAP<int, std::vector<int>>& GetSolidInterfaceCells() {
    return __solid_interface_cells;
  }
  HASHMAP<int, std::vector<int>>& GetFluidInterfaceCells() {
    return __fluid_interface_cells;
  }
  HASHSET<int>& GetFluidCells() {
    return __fluid_cells;
  }
  HASHSET<int>& GetSolidCells() {
    return __solid_cells;
  }
  HASHSET<int>& GetInterfaceNodes() {
    return __interface_nodes;
  }
  HASHMAP<int, int>& GetFluidG2L() {
    return __fluid_g2l;
  }
  HASHMAP<int, int>& GetSolidG2L() {
    return __solid_g2l;
  }
  std::vector<int>& GetFluidL2G() {
    return __fluid_l2g;
  }
  std::vector<int>& GetSolidL2G() {
    return __solid_l2g;
  }

  const HASHMAP<int, std::vector<int>>& GetSolidInterfaceCells() const {
    return __solid_interface_cells;
  }
  const HASHMAP<int, std::vector<int>>& GetFluidInterfaceCells() const {
    return __fluid_interface_cells;
  }
  const HASHSET<int>& GetFluidCells() const {
    return __fluid_cells;
  }
  const HASHSET<int>& GetSolidCells() const {
    return __solid_cells;
  }
  const HASHSET<int>& GetInterfaceNodes() const {
    return __interface_nodes;
  }
  const HASHMAP<int, int>& GetFluidG2L() const {
    return __fluid_g2l;
  }
  const HASHMAP<int, int>& GetSolidG2L() const {
    return __solid_g2l;
  }
  const std::vector<int>& GetFluidL2G() const {
    return __fluid_l2g;
  }
  const std::vector<int>& GetSolidL2G() const {
    return __solid_l2g;
  }

  bool SolidInterfaceCell(int iq) const {
    return (__solid_interface_cells.find(iq) != __solid_interface_cells.end());
  }
  bool FluidInterfaceCell(int iq) const {
    return (__fluid_interface_cells.find(iq) != __fluid_interface_cells.end());
  }

  bool SolidCell(int iq) const {
    return (__solid_cells.find(iq) != __solid_cells.end());
  }
  bool FluidCell(int iq) const {
    return (__fluid_cells.find(iq) != __fluid_cells.end());
  }
  bool InterfaceNode(int iq) const {
    return (__interface_nodes.find(iq) != __interface_nodes.end());
  }

  bool SolidNode(int i) const {
    return (__solid_g2l.find(i) != __solid_g2l.end());
  }
  bool FluidNode(int i) const {
    return (__fluid_g2l.find(i) != __fluid_g2l.end());
  }

  const std::vector<int> SolidInterfaceNodes(int iq) const {
    assert(SolidInterfaceCell(iq));
    return __solid_interface_cells.find(iq)->second;
  }
  const std::vector<int> FluidInterfaceNodes(int iq) const {
    assert(FluidInterfaceCell(iq));
    return __fluid_interface_cells.find(iq)->second;
  }

  void DeleteTestFunctionsVector(const LocalVector&      F,
                                 const std::vector<int>& indices,
                                 int                     comp) const {
    LocalVector& __F= const_cast<LocalVector&>(F);
    for (int i= 0; i < indices.size(); ++i)
      __F(indices[i], comp)= 0.0;
  }

  void DeleteTestFunctionMatrix(EntryMatrix& E, int index, int comp) const {
    for (int j= 0; j < E.Mdof(); ++j)
      E(index, j, comp, comp)= 0.0;
  }
  void DeleteAnsatzFunctionMatrix(EntryMatrix& E, int index, int comp) const {
    for (int j= 0; j < E.Mdof(); ++j)
      E(j, index, comp, comp)= 0.0;
  }

  void DeleteTestFunctionsMatrix(EntryMatrix&            E,
                                 const std::vector<int>& indices,
                                 int                     comp) const {
    for (int i= 0; i < indices.size(); ++i)
      DeleteTestFunctionMatrix(E, indices[i], comp);
  }
  void DeleteAnsatzFunctionsMatrix(EntryMatrix&            E,
                                   const std::vector<int>& indices,
                                   int                     comp) const {
    for (int i= 0; i < indices.size(); ++i)
      DeleteAnsatzFunctionMatrix(E, indices[i], comp);
  }
};

template <class BASE>
class AleDiscretization : public AleBaseDiscretization, virtual public BASE {
  void LocalToGlobal(GlobalVector& f, const LocalVector& F, int iq, double s) const {
    if (FluidInterfaceCell(iq))
      for (int i= 0; i < __DEL_F.size(); ++i)
        DeleteTestFunctionsVector(F, FluidInterfaceNodes(iq), __DEL_F[i]);
    if (SolidInterfaceCell(iq))
      for (int i= 0; i < __DEL_S.size(); ++i)
        DeleteTestFunctionsVector(F, SolidInterfaceNodes(iq), __DEL_S[i]);

    BasicDiscretization::LocalToGlobal(f, F, iq, s);
  }

  void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const {
    if (SolidInterfaceCell(iq))  // In Primal, delete Test-functions
      for (int i= 0; i < __DEL_S.size(); ++i)
        DeleteTestFunctionsMatrix(E, SolidInterfaceNodes(iq), __DEL_S[i]);

    if (FluidInterfaceCell(iq))  // In Primal, delete Test-functions
      for (int i= 0; i < __DEL_F.size(); ++i)
        DeleteTestFunctionsMatrix(E, FluidInterfaceNodes(iq), __DEL_F[i]);

    BASE::LocalToGlobal(A, E, iq, s);
  }

  void GlobalToLocal(LocalVector& U, const GlobalVector& u, int iq) const {
    this->GlobalToLocalFirst(U, u, iq);
    this->GlobalToLocalData(iq);
  }

  void GlobalToLocalFirst(LocalVector& U, const GlobalVector& u, int iq) const {
    BASE::GlobalToLocalSingle(U, u, iq);

    /* if (FluidInterfaceCell(iq)) */
    /* 	for (int i=0;i<__DEL_F.size();++i) */
    /* 	  DeleteTestFunctionsVector(U,FluidInterfaceNodes(iq),__DEL_F[i]); */
    /* if (SolidInterfaceCell(iq)) */
    /* 	for (int i=0;i<__DEL_S.size();++i) */
    /* 	  DeleteTestFunctionsVector(U,SolidInterfaceNodes(iq),__DEL_S[i]); */
  }
};

}  // namespace Gascoigne

/*----------------------------   alebasediscretization.h     ---------------------------*/
/* end of #ifndef __alebasediscretization_H */
#endif
/*----------------------------   alebasediscretization.h     ---------------------------*/
