#include "backup.h"
#include "meshinterpolator.h"
#include "q12d.h"
#include "q22d.h"
#include "stdtimesolver.h"
#include "domainrighthandside.h"

using namespace std;

namespace Gascoigne
{
/**********************************************************/

MeshInterpolator::MeshInterpolator() : _MA(NULL), _DI(NULL)
{
}

/**********************************************************/

MeshInterpolator::~MeshInterpolator()
{
  if (_DI)
  {
    delete GetDiscretization();
  }
  if (_MA)
  {
    delete GetMeshAgent();
  }
}

/**********************************************************/

void MeshInterpolator::CheckCell(int oldNumber, int newNumber)
{
  if (_Old->sleep(oldNumber) && _New->sleep(newNumber))
  {
    for (int i=0; i<_Old->nchilds(oldNumber); i++)
    {
      int oi = _Old->child(oldNumber,i);
      int ni = _New->child(newNumber,i);
      CheckCell(oi,ni);
    }
  }
  else
  {
    for (int i=0; i<_Old->nodes_per_cell(oldNumber); i++)
    {
      int oi = _Old->vertex_of_cell(oldNumber,i);
      _NewNodeNumber[oi] = _New->vertex_of_cell(newNumber,i);
    }
    if (!_VecNew.empty() && _Old->sleep(oldNumber) && !_New->sleep(newNumber))
    {
      _ToBeRefNew.push_back(newNumber);
    }
    else if (!_VecOld.empty() && !_Old->sleep(oldNumber) && _New->sleep(newNumber))
    {
      _ToBeRef.push_back(oldNumber);
    }
  }
}

/**********************************************************/

void MeshInterpolator::Coarsen(int newNumber)
{
  for (int i=0; i<_New->nchilds(newNumber); i++)
  {
    if (_New->sleep(_New->child(newNumber,i)))
    {
      Coarsen(_New->child(newNumber,i));
    }
  }
  int npc = _New->nodes_per_cell(newNumber);
  for (int i=0; i<npc; i++)
  {
    int k = _New->vertex_of_cell(newNumber,i);
    for (int j=0; j<npc; j++)
    {
      int l = _New->vertex_of_cell(_New->child(newNumber,i),j);
      for (int s=0; s<_VecInt.size(); s++)
      {
	double w = _weights(i,j);
	_VecInt[s].add_node(k,w,l);
      }
    }
  }
  for (int i=0; i<npc; i++)
  {
    int k = _New->vertex_of_cell(newNumber,i);
    for (int j=0; j<npc; j++)
    {
      if (j!=i)
      {
	int ni = _New->child(newNumber,i);
        for (int s=0; s<_VecInt.size(); s++)
        {
	  int l = _New->vertex_of_cell(ni,j);
	  _VecInt[s].zero_node(l);
        }
      }
    }
  }
}

/**********************************************************/

void MeshInterpolator::Distribute(int oldNumber, int newNumber)
{
  if (_Old->sleep(oldNumber) && _New->sleep(newNumber))
  {
    for (int i=0; i<_Old->nchilds(oldNumber); i++)
    {
      Distribute(_Old->child(oldNumber,i),_New->child(newNumber,i));
    }
  }
  else if (!_Old->sleep(oldNumber) && _New->sleep(newNumber))
  {
    Coarsen(newNumber);
  }
  else if (_Old->sleep(oldNumber) || _New->sleep(newNumber))
  {
    cerr << "Das darf gar nicht passieren!!!" << endl;
    abort();
  }
}

/**********************************************************/

void MeshInterpolator::RefineAndInterpolate(HierarchicalMesh* Mesh, vector<GlobalVector>& u, const IntVector& refine, vector<bool>& done)
{
  IntVector coarse(0);
  int oldcells = Mesh->ncells();
  Mesh->refine(refine,coarse);

  int nn = Mesh->nnodes();
  done.resize(nn,false);
  for (int s=0; s<u.size(); s++)
  {
    u[s].resize(nn,0.);
  }
  set<int> fathers;
  for (int cell=oldcells; cell<Mesh->ncells(); cell++)
  {
    fathers.insert(Mesh->Vater(cell));
  }
  for (set<int>::const_iterator p = fathers.begin(); p!=fathers.end(); p++)
  {
    int npc = Mesh->nodes_per_cell(*p);
    for (int i=0; i<npc; i++)
    {
      int ci = Mesh->child(*p,i);
      for (int j=0; j<npc; j++)
      {
	int l = Mesh->vertex_of_cell(ci,i);
        if (j!=i && !done[l])
        {
	  int k = Mesh->vertex_of_cell(ci,j);
          for (int s=0; s<u.size(); s++)
          {
	    double w = _weights(i,j);
	    u[s].add_node(k,w,l);
          }
        }
      }
    }
    for (int i=0; i<npc; i++)
    {
      int ci = Mesh->child(*p,i);
      for (int j=0; j<npc; j++)
      {
        if (j!=i)
        {
	  int k = Mesh->vertex_of_cell(ci,j);
          done[k] = true;
        }
      }
    }
  }
}

/**********************************************************/

void MeshInterpolator::AddVectorIntermediate(const GlobalVector& u)
{
  _VecInt.push_back(u);
}

/**********************************************************/

void MeshInterpolator::AddVectorOld(const GlobalVector& u)
{
  _VecOld.push_back(u);
  //GetDiscretization()->HNAverage(u);
}

/**********************************************************/

void MeshInterpolator::AddVectorNew(const GlobalVector& u)
{
  _VecNew.push_back(u);
  //GetDiscretization()->HNAverage(u);
}

/**********************************************************/

void MeshInterpolator::BasicInit(DiscretizationInterface* DI, MeshAgentInterface* MA, const string& name)
{
  // Klassenvariablen initialisieren
  _BaseCells.clear();
  _ToBeRef.clear();
  _ToBeRefNew.clear();
  _NewNodeNumber.clear();
  _VecInt.clear();
  _VecOld.clear();
  _VecNew.clear();

  // Original-Solver und -MeshAgent speichern
  GetOriginalMeshAgentPointer() = MA;

  int dim = GetOriginalMeshAgent()->GetMesh(0)->dimension();
  _ODI = DI;

  // neuen MeshAgent anlegen
  _MA = new MeshAgent;

  GetMeshAgent()->GetShapes2d() = GetOriginalMeshAgent()->GetShapes2d();
  GetMeshAgent()->GetShapes3d() = GetOriginalMeshAgent()->GetShapes3d();
  GetMeshAgent()->SetDefaultValues(dim,name+".gup",0);
  GetMeshAgent()->BasicInit(NULL);

  // neue Discretization anlegen
  _DI = new Q12d;
  GetDiscretization()->BasicInit(NULL);
  GetDiscretization()->ReInit(GetMeshAgent()->GetMesh(0));

  const Q1* Q1DP = dynamic_cast<const Q1*>(GetDiscretization());
  if (Q1DP)
  {
    int size = static_cast<int>(pow(2.,static_cast<double>(dim)));
    _weights.resize(size,size);
    _weights = Q1DP->GetLocalInterpolationWeights();
  }
  else
  {
    cerr << "nur Q1-Elemente unterstuetzt!" << endl;
    abort();
  }

  // Indizes aufgrund unterschiedlicher Nummerierung im HierarchicalMesh tauschen!
  for (int i=0; i<_weights.m(); i++)
  {
    swap(_weights(i,2),_weights(i,3));
    if (dim==3)
    {
      swap(_weights(i,6),_weights(i,7));
    }
  }
  for (int j=0; j<_weights.n(); j++)
  {
    swap(_weights(2,j),_weights(3,j));
    if (dim==3)
    {
      swap(_weights(6,j),_weights(7,j));
    }
  }

  // Pointer auf die HierarchicalMeshs holen
  _Old = GetOriginalMeshAgent()->GetHierarchicalMesh();
  _New = GetMeshAgent()->GetHierarchicalMesh();
  assert(_Old);
  assert(_New);
}

/**********************************************************/

       class ProjectionRightHandSide : public Gascoigne::DomainRightHandSide
         {
           protected:
	   mutable Gascoigne::FemFunction _U;
	   
	 public:
	   
	   ProjectionRightHandSide() : Gascoigne::DomainRightHandSide()  {}
	   ~ProjectionRightHandSide() { }
	   
	   int GetNcomp() const { return 1;}//_U.size(); }
	   std::string GetName() const { return "ProjectionRightHandSide"; }
	   
	   void SetFemData(Gascoigne::FemData& q) const
	   {
	     assert(q.count("U")==1);
	     _U = q["U"];
	   }
	   
	   void operator()(Gascoigne::VectorIterator b, const Gascoigne::TestFunction& N, 
			   const Gascoigne::Vertex2d& v) const 
	   {
	     for (int i=0; i<_U.size(); i++)
	       b[i] += _U[i].m() * N.m();
	   }
       };

/**********************************************************/

void MeshInterpolator::RhsForProjection(GlobalVector& f, const GlobalVector& u)
{
  AddVectorNew(u);

  vector<bool> doneOld(_Old->nnodes(),true),doneNew(_New->nnodes(),true);
  HierarchicalMesh* Mesh;
  if (_Old->ncells()<_New->ncells())
  {
    Mesh = _Old;
  }
  else
  {
    Mesh = _New;
  }
  for (int c=0; c<Mesh->ncells(); c++)
  {
    if (Mesh->level(c)==0)
    {
      _BaseCells.insert(c);
    }
  }

  do
  {
    _NewNodeNumber.resize(_Old->nnodes(),-1);
    _ToBeRef.clear();
    _ToBeRefNew.clear();
    for (set<int>::const_iterator pbc = _BaseCells.begin(); pbc!=_BaseCells.end(); pbc++)
    {
      CheckCell(*pbc,*pbc);
    }
    if (!_ToBeRef.empty())
    {
      RefineAndInterpolate(_Old,_VecOld,_ToBeRef,doneOld);
    }
    if (!_ToBeRefNew.empty())
    {
      RefineAndInterpolate(_New,_VecNew,_ToBeRefNew,doneNew);
    }
  }
  while (!_ToBeRef.empty() || !_ToBeRefNew.empty());

  GetMeshAgent()->global_refine(0);
  GetDiscretization()->ReInit(GetMeshAgent()->GetMesh(0));

  ProjectionRightHandSide  PD;
  
  GlobalVector _help;
  _help.ReInit(u.ncomp(),u.n());
  GetDiscretization()->AddNodeVector("U",&u);
  GetDiscretization()->Rhs(_help,PD,1.);
  GetDiscretization()->DeleteNodeVector("U");

  AddVectorIntermediate(_help);
  for (set<int>::const_iterator pbc = _BaseCells.begin(); pbc!=_BaseCells.end(); pbc++)
  {
    Distribute(*pbc,*pbc);
  }

  f.zero();

  assert(_VecInt.size()==1);
  GlobalVector help = _VecInt[0];

  for (int i=0; i<f.n(); i++)
  {
    for (int c=0; c<f.ncomp(); c++)
    {
      f(i,c) = help(_NewNodeNumber[i],c);
    }
  }
  GetOriginalDiscretization()->HNDistribute(f);

  _Old = NULL;
  _New = NULL;
}

/**********************************************************/
}

