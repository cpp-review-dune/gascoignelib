#include  "basicmeshinterpretor.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
BasicMeshInterpretor::BasicMeshInterpretor() : MeshInterpretorInterface()
{
}

/* ----------------------------------------- */

BasicMeshInterpretor::~BasicMeshInterpretor()
{
}

/* ----------------------------------------- */

void BasicMeshInterpretor::HNAverageData() const
{
  const GlobalNodeData& gd = GetGlobalData().GetNodeData();
  GlobalNodeData::const_iterator p=gd.begin();
  for(; p!=gd.end(); p++)
    {
      GlobalVector* v = const_cast<GlobalVector*>(p->second);
      HNAverage(*v);
//      cerr << "HNAverage " << p->first << endl;
    }
}

/* ----------------------------------------- */

void BasicMeshInterpretor::HNZeroData() const
{
  const GlobalNodeData& gd = GetGlobalData().GetNodeData();
  GlobalNodeData::const_iterator p=gd.begin();
  for(; p!=gd.end(); p++)
    {
      GlobalVector* v = const_cast<GlobalVector*>(p->second);
      HNZero(*v);
//      cerr << "HNZero " << p->first << endl;
    }
}

/* ----------------------------------------- */

void BasicMeshInterpretor::GlobalToLocalData(int iq) const
{
  const GlobalNodeData& gd = GetGlobalData().GetNodeData();
  __Q.clear();
  GlobalNodeData::const_iterator p=gd.begin();
  for(; p!=gd.end(); p++)
    {
      GlobalToLocalSingle(__Q[p->first],*p->second,iq);
    }
}

/* ----------------------------------------- */

void BasicMeshInterpretor::GlobalToGlobalData() const
{
  const GlobalParameterData& gd = GetGlobalData().GetParameterData();
  __qq.clear();
  GlobalParameterData::const_iterator p=gd.begin();
  for(; p!=gd.end(); p++)
    {
      __qq.insert(make_pair(p->first,*p->second));
    }
}

/* ----------------------------------------- */

void BasicMeshInterpretor::GlobalToLocalSingle(LocalVector& U, const GlobalVector& u, int iq) const
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

void BasicMeshInterpretor::LocalToGlobal(GlobalVector& f, const LocalVector& F, int iq, double s) const
{
  IntVector indices = GetLocalIndices(iq);
  for(int ii=0; ii<indices.size(); ii++) 
    {
      int i = indices[ii];
      f.add_node(i,s,ii,F);
    }
}

/* ----------------------------------------- */

void BasicMeshInterpretor::LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
{
  IntVector indices = GetLocalIndices(iq);
  IntVector::const_iterator  start = indices.begin();
  IntVector::const_iterator  stop  = indices.end();
  A.entry(start,stop,__E,s);
}

/* ----------------------------------------- */

void BasicMeshInterpretor::VertexTransformation(const nvector<Vertex2d>& p, const Vertex2d& p0, Vertex2d& tp) const
{
  vector<double> DT(4);

  double x0 = p[0].x();
  double x1 = p[1].x();
  double x2 = p[2].x();
  double x3 = p[3].x();

  double y0 = p[0].y();
  double y1 = p[1].y();
  double y2 = p[2].y();
  double y3 = p[3].y();

  double a1,a2,b1,b2,c1,c2,d1,d2;
  
  
  d1 = x0; a1 = x1-x0; b1 = x3-x0; c1 = x2-x1+x0-x3;
  d2 = y0; b2 = y3-y0; a2 = y1-y0; c2 = y2-y1-y3+y0;
  
  tp.x()=0.5;
  tp.y()=0.5;
  
  Vertex2d dp;
  Vertex2d res;

  res.x() = p0.x() - ( a1*tp.x() + b1*tp.y() + c1*tp.x()*tp.y() + d1);
  res.y() = p0.y() - ( a2*tp.x() + b2*tp.y() + c2*tp.x()*tp.y() + d2);
  
  do 
    {
      DT[0] = a1+c1*tp.y();
      DT[1] = b1+c1*tp.x();

      DT[2] = a2+c2*tp.y();
      DT[3] = b2+c2*tp.x();

      double det = DT[0]*DT[3]-DT[1]*DT[2];
      assert(fabs(det)>=1.e-20);
      
      dp.x()=(DT[3]*res.x()-DT[2]*res.y())/det;
      dp.y()=(DT[0]*res.y()-DT[1]*res.x())/det;

      tp.x()+=dp.x();
      tp.y()+=dp.y();

      res.x() = p0.x() - ( a1*tp.x() + b1*tp.y() + c1*tp.x()*tp.y() + d1);
      res.y() = p0.y() - ( a2*tp.x() + b2*tp.y() + c2*tp.x()*tp.y() + d2);
      assert(tp.norm_l8()<=2);
    } while (res.norm_l8()>1.e-14);
}

/*-----------------------------------------*/

}
