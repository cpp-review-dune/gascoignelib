#include  "coarsehierarchicalmesh3d.h"
#include  "set2vec.h"

using namespace Gascoigne;

/*---------------------------------------------------*/

CoarseHierarchicalMesh3d::CoarseHierarchicalMesh3d
(const HierarchicalMesh3d& HM) :  HierarchicalMesh3d(HM){}

/*---------------------------------------------------*/

void CoarseHierarchicalMesh3d::BasicInit(int depth)
{
  if (depth==1)
    {
      loop(cn2o);
    }
  else if (depth==2)
    {
      IntVector   cn2oA, cn2oB;
      IntVector   co2n2;
      loop(cn2oA);
      co2n2 = co2n;
      loop(cn2oB);
      IntVector   co2n3(co2n2.size());
      for (int i=0; i<co2n2.size(); i++)
	{
	  if (co2n2[i]<0) 
	    {
	      co2n3[i] = co2n2[i];
	    }
	  else
	    {
	      assert(co2n2[i]<co2n.size());
	      co2n3[i] = co2n[co2n2[i]];
	    }
	}
      co2n.resize(co2n3.size());
      co2n = co2n3;

      cn2o.resize(ncells());
      cn2o = -1;
      for (int i=0; i<cn2oB.size(); i++)
	{
	  int j = cn2oB[i];
	  if(j>=0)  
	    {
	      int k = cn2oA[j];
	      if(k>=0)  
		{
		  cn2o[i] = k;
		}
	    }
	}
    }
  else assert(0);
}

/*---------------------------------------------------*/

void CoarseHierarchicalMesh3d::loop(IntVector& dst)
{
  global_coarse3d();
  dst.resize(ncells());
  dst = -1;
  for (int i=0; i<co2n.size(); i++)
    {
      int j = co2n[i];
      assert(j<ncells());
      if (j>=0)	dst[j] = i;
    }
}

/*---------------------------------------------------*/

void CoarseHierarchicalMesh3d::refine
(const IntVector& cell_ref_old, const IntVector& cell_coarse_old)
{
  CellRefList   .clear();
  CellCoarseList.clear();

  IntVector cell_ref(0),cell_coarse(0);

  for(int i=0; i<cell_ref_old.size(); i++) 
    {
      int newc = co2n[cell_ref_old[i]];
      if (newc>=0)
	cell_ref.push_back(newc);
    }

  for(int i=0; i<cell_coarse_old.size(); i++) 
    {
      int newc = co2n[cell_coarse_old[i]];
      if (newc>=0)
	cell_coarse.push_back(newc);
    }

  _refine3d(CellRefList,CellCoarseList,cell_ref,cell_coarse);
}

/*---------------------------------------------------*/

void CoarseHierarchicalMesh3d::GetRefinedList(IntVector& ref)
{
  ref.resize(0);
  IntVector ref2;
  Set2Vec(ref2,CellRefList);
  for (int i=0; i<ref2.size(); i++)
    {
      int j = cn2o[ref2[i]];
      if (j>=0)
	ref.push_back(j);
    }
}

/*---------------------------------------------------*/

void CoarseHierarchicalMesh3d::GetCoarsedList(IntVector& coarse)
{
  coarse.resize(0);
  IntVector coarse2;
  Set2Vec(coarse2,CellCoarseList);
  for (int i=0; i<coarse2.size(); i++)
    {
      int j = cn2o[coarse2[i]];
      if (j>=0)
	coarse.push_back(j);
    }
}

/*---------------------------------------------------*/
