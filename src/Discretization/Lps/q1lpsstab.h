#ifndef  __Q1LpsStab_h
#define  __Q1LpsStab_h

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1LpsStab

////
////
/////////////////////////////////////////////

#include  "patchmeshinterpretor.h"
#include  "hnstructureinterface.h"

namespace Gascoigne
{

/*----------------------------------------------*/

class Q1LpsStab : public PatchMeshInterpretor
{
 protected:

  const HNStructureInterface* HN;

  nvector<int> GetLocalIndices(int iq) const {
    return *GetPatchMesh()->IndicesOfPatch(iq);
  }
  void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;

 public:

  Q1LpsStab() : PatchMeshInterpretor() {};
  int n() const { return GetMesh()->nnodes();}
  void BasicInit(const Gascoigne::ParamFile* paramfile, const HNStructureInterface*);
};

/*----------------------------------------------*/

class Q1LpsStab2d : public Q1LpsStab
{
 protected:

  void Transformation(FemInterface::Matrix& T, int iq) const;

 public:

  Q1LpsStab2d() : Q1LpsStab() {};
  void BasicInit(const Gascoigne::ParamFile* paramfile, const HNStructureInterface*);
};

/*----------------------------------------------*/

class Q1LpsStab3d : public Q1LpsStab
{
 protected:

  void Transformation(FemInterface::Matrix& T, int iq) const;

 public:

  Q1LpsStab3d() : Q1LpsStab() {};
  void BasicInit(const Gascoigne::ParamFile* paramfile, const HNStructureInterface*);
};

/*----------------------------------------------*/
}
#endif
