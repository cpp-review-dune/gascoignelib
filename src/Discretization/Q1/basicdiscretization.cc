#include  "basicdiscretization.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
BasicDiscretization::BasicDiscretization() : DiscretizationInterface(), __MP(NULL)
{
}

/* ----------------------------------------- */

BasicDiscretization::~BasicDiscretization()
{
}

/* ----------------------------------------- */

void BasicDiscretization::HNAverageData() const
{
  const GlobalData& gd = GetDataContainer().GetNodeData();
  GlobalData::const_iterator p=gd.begin();
  for(; p!=gd.end(); p++)
    {
      GlobalVector* v = const_cast<GlobalVector*>(p->second);
      HNAverage(*v);
//      cerr << "HNAverage " << p->first << endl;
    }
}

/* ----------------------------------------- */

void BasicDiscretization::HNZeroData() const
{
  const GlobalData& gd = GetDataContainer().GetNodeData();
  GlobalData::const_iterator p=gd.begin();
  for(; p!=gd.end(); p++)
    {
      GlobalVector* v = const_cast<GlobalVector*>(p->second);
      HNZero(*v);
//      cerr << "HNZero " << p->first << endl;
    }
}

/* ----------------------------------------- */

void BasicDiscretization::GlobalToLocalData(int iq) const
{
  const GlobalData& gnd = GetDataContainer().GetNodeData();
  __QN.clear();
  GlobalData::const_iterator p=gnd.begin();
  for(; p!=gnd.end(); p++)
    {
      GlobalToLocalSingle(__QN[p->first],*p->second,iq);
    }

  const GlobalData& gcd = GetDataContainer().GetCellData();
  __QC.clear();
  GlobalData::const_iterator q=gcd.begin();
  for(; q!=gcd.end(); q++)
    {
      GlobalToLocalCell(__QC[q->first],*q->second,iq);
    }
}

/* ----------------------------------------- */

void BasicDiscretization::GlobalToGlobalData() const
{
  const GlobalParameterData& gpd = GetDataContainer().GetParameterData();
  __QP.clear();
  GlobalParameterData::const_iterator p=gpd.begin();
  for(; p!=gpd.end(); p++)
    {
      __QP.insert(make_pair(p->first,*p->second));
    }
}

/* ----------------------------------------- */

void BasicDiscretization::GlobalToLocalSingle(LocalVector& U, const GlobalVector& u, int iq) const
{
  IntVector indices = GetLocalIndices(iq);
  U.ReInit(u.ncomp(),indices.size());
  for(int ii=0; ii<indices.size(); ii++) 
    {
      int i = indices[ii];
      U.equ_node(ii,i,u);
    }
}

/* ----------------------------------------- */

void BasicDiscretization::GlobalToLocalCell(LocalVector& U, const GlobalVector& u, int iq) const
{
  U.ReInit(u.ncomp(),1);
  for(int c=0;c<u.ncomp();++c)
    {
      U(0,c) = u(iq,c);
    }
}

/* ----------------------------------------- */

void BasicDiscretization::LocalToGlobal(GlobalVector& f, const LocalVector& F, int iq, double s) const
{
  IntVector indices = GetLocalIndices(iq);
  for(int ii=0; ii<indices.size(); ii++) 
    {
      int i = indices[ii];
      f.add_node(i,s,ii,F);
    }
}

/* ----------------------------------------- */

void BasicDiscretization::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  IntVector indices = GetLocalIndices(iq);
  IntVector::const_iterator  start = indices.begin();
  IntVector::const_iterator  stop  = indices.end();
  A.entry(start,stop,E,s);
}

/*-----------------------------------------*/

}
