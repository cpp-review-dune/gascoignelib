#include "edgeinfocontainer.h"

using namespace std;
using namespace Gascoigne;

/**********************************************************/

template<int DIM>
EdgeInfoContainer<DIM>::~EdgeInfoContainer<DIM>()
{
  for (int i=0; i<size(); i++)
    {
      if ((*this)[i])
	{
	  delete (*this)[i];
	  (*this)[i] = NULL;
	}
    }
  resize(0);
}

/**********************************************************/

template<int DIM>
void EdgeInfoContainer<DIM>::basicInit(const HierarchicalMesh* HM, int ncomp)
{
  resize(0);
  _HMP   = HM;
  _ncomp = ncomp;

  resize(_HMP->nedges());
  for (int i=0; i<size(); i++)
    {
      (*this)[i]=NULL;
    }
}

/**********************************************************/

template<int DIM>
void EdgeInfoContainer<DIM>::showStatistics() const
{
  for (int i=0; i<size(); i++)
    {
      if ((*this)[i]!=NULL)
	{
	  cout << "Edge " << i << ": ";
	  (*this)[i]->showStatistics();
	}
    }
}

/**********************************************************/

template<int DIM>
const HierarchicalMesh* EdgeInfoContainer<DIM>::getMesh() const
{
  return _HMP;
}

/**********************************************************/

void EdgeInfoContainer<2>::modifyHanging()
{
  const QuadLawAndOrder& QLAO = dynamic_cast<const HierarchicalMesh2d*>(_HMP)->QuadLawOrder();
  fixarray<2,int>        vertexes;
  LocalVector            lu,lul,lur;
  
  lu.ReInit(_ncomp,2);
  lu.zero();
  lul.ReInit(_ncomp,2);
  lul.zero();
  lur.ReInit(_ncomp,2);
  lur.zero();
  
  for (int i=0; i<size(); i++)
    {
      if ((*this)[i]!=NULL && (*this)[i]->getCount()==1 && (*this)[i]->getEdge().slave()!=-1)
	{
	  const Edge& edge = (*this)[i]->getEdge();
	  vertexes = (*this)[i]->vertex();
	  
	  int left  = QLAO.GlobalChildEdge(vertexes,edge.slave(),0);
	  int right = QLAO.GlobalChildEdge(vertexes,edge.slave(),1);
	  
	  for (int c=0; c<_ncomp; c++)
	    {
	      lu(0,c) = ((*this)[right]->getValue())(1,c);
	      lu(1,c) = ((*this)[left]->getValue())(0,c);
	      
	      lul(0,c) = ((*this)[i]->getValue())(1,c);
	      lul(1,c) = 0.5 * (((*this)[i]->getValue())(0,c) + ((*this)[i]->getValue())(1,c));
	      
	      lur(0,c) = 0.5 * (((*this)[i]->getValue())(0,c) + ((*this)[i]->getValue())(1,c));
	      lur(1,c) = ((*this)[i]->getValue())(0,c);
	    }
	  (*this)[i]->addNodes(lu);
	  (*this)[left]->addNodes(lul);
	  (*this)[right]->addNodes(lur);
	}
    }
}

/**********************************************************/

void EdgeInfoContainer<3>::modifyHanging()
{
  const HexLawAndOrder& HLAO = dynamic_cast<const HierarchicalMesh3d*>(_HMP)->HexLawOrder();
  fixarray<4,int>       vertexes;
  LocalVector           lugr,lulu,luru,lulo,luro;
  
  lugr.ReInit(_ncomp,4);
  lugr.zero();
  lulu.ReInit(_ncomp,4);
  lulu.zero();
  luru.ReInit(_ncomp,4);
  luru.zero();
  lulo.ReInit(_ncomp,4);
  lulo.zero();
  luro.ReInit(_ncomp,4);
  luro.zero();
  
  for (int i=0; i<size(); i++)
    {
      if ((*this)[i]!=NULL && (*this)[i]->getCount()==1 && (*this)[i]->getEdge().slave()!=-1)
	{
	  const Edge& quad = (*this)[i]->getEdge();
	  vertexes = (*this)[i]->vertex();
	  
	  int lu = HLAO.GlobalChildFace(vertexes,quad.master(),0);
	  int ru = HLAO.GlobalChildFace(vertexes,quad.master(),1);
	  int lo = HLAO.GlobalChildFace(vertexes,quad.master(),3);
	  int ro = HLAO.GlobalChildFace(vertexes,quad.master(),2);
	  
	  LocalVector help = (*this)[i]->getValue();
	  
	  for (int c=0; c<_ncomp; c++)
	    {
	      lugr(0,c) = ((*this)[lu]->getValue())(0,c);
	      lugr(1,c) = ((*this)[ru]->getValue())(1,c);
	      lugr(2,c) = ((*this)[ro]->getValue())(2,c);
	      lugr(3,c) = ((*this)[lo]->getValue())(3,c);
	      
	      lulu(0,c) = help(0,c);
	      lulu(1,c) = 0.5 * (help(0,c) + help(1,c));
	      lulu(2,c) = 0.25 * (help(0,c) + help(1,c) + help(2,c) + help(3,c));
	      lulu(3,c) = 0.5 * (help(0,c) + help(3,c));
	      
	      luru(0,c) = 0.5 * (help(0,c) + help(1,c));
	      luru(1,c) = help(1,c);
	      luru(2,c) = 0.5 * (help(1,c) + help(2,c));
	      luru(3,c) = 0.25 * (help(0,c) + help(1,c) + help(2,c) + help(3,c));
	      
	      lulo(0,c) = 0.5 * (help(0,c) + help(3,c));
	      lulo(1,c) = 0.25 * (help(0,c) + help(1,c) + help(2,c) + help(3,c));
	      lulo(2,c) = 0.5 * (help(2,c) + help(3,c));
	      lulo(3,c) = help(3,c);
	      
	      luro(0,c) = 0.25 * (help(0,c) + help(1,c) + help(2,c) + help(3,c));
	      luro(1,c) = 0.5 * (help(1,c) + help(2,c));
	      luro(2,c) = help(2,c);
	      luro(3,c) = 0.5 * (help(2,c) + help(3,c));
	    }
	  (*this)[i]->addNodes(lugr);
	  (*this)[lu]->addNodes(lulu);
	  (*this)[ru]->addNodes(luru);
	  (*this)[lo]->addNodes(lulo);
	  (*this)[ro]->addNodes(luro);
	}
    }
}

/**********************************************************/

template class EdgeInfoContainer<2>;
template class EdgeInfoContainer<3>;
