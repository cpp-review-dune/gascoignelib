#include  "basicmeshinterpretor.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

BasicMeshInterpretor::BasicMeshInterpretor() : MeshInterpretorInterface(), __FEM(NULL), __INT(NULL)
{
}

/* ----------------------------------------- */

BasicMeshInterpretor::~BasicMeshInterpretor()
{
  if(__FEM==NULL) {delete __FEM; __FEM=NULL;}
  if(__INT==NULL) {delete __INT; __INT=NULL;}
}

/* ----------------------------------------- */

void BasicMeshInterpretor::GlobalToLocalData(int iq) const
{
  const GlobalNodeData& gd = GetGlobalData().GetNodeData();
  GlobalNodeData::const_iterator p=gd.begin();
  __Q.resize(gd.size());
  int i=0;
  assert(gd.size()==__Q.size());
  for(; p!=gd.end(); p++)
    {
      const GlobalVector& q=**p;
      GlobalToLocalSingle(__Q[i++],q,iq);
    }
}

/* ----------------------------------------- */

void BasicMeshInterpretor::GlobalToLocalSingle(LocalVector& U, const GlobalVector& u, int iq) const
{
  nvector<int> indices = GetLocalIndices(iq);
  U.ReInit(u.ncomp(),indices.size());
  for(int ii=0; ii<indices.size(); ii++) 
    {
      int i = indices[ii];
      U.equ_node(ii,i,u);
    }
}

/* ----------------------------------------- */

void BasicMeshInterpretor::LocalToGlobal(GlobalVector& f, const LocalVector& F, int iq, double s) const
{
  nvector<int> indices = GetLocalIndices(iq);
  for(int ii=0; ii<indices.size(); ii++) 
    {
      int i = indices[ii];
      f.add_node(i,s,ii,F);
    }
}

/* ----------------------------------------- */

void BasicMeshInterpretor::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  nvector<int> indices = GetLocalIndices(iq);
  nvector<int>::const_iterator  start = indices.begin();
  nvector<int>::const_iterator  stop  = indices.end();
  A.entry(start,stop,__E,s);
}
