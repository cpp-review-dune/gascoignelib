#include "backup.h"
#include "meshinterpolator.h"
#include "q1.h"
#include "q2.h"
#include "stdtimesolver.h"

using namespace std;

namespace Gascoigne
{
/**********************************************************/

MeshInterpolator::MeshInterpolator() : _MA(NULL), _SI(NULL)
{
}

/**********************************************************/

MeshInterpolator::~MeshInterpolator()
{
  if (GetSolverPointer())
  {
    delete GetSolver();
  }
  if (GetMeshAgentPointer())
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
	int l = Mesh->vertex_of_cell(ci,j);
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

void MeshInterpolator::AddVectorIntermediate(GlobalVector u)
{
  _VecInt.push_back(u);
}

/**********************************************************/

void MeshInterpolator::AddVectorOld(GlobalVector u)
{
  assert(GetOriginalSolverPointer());
  _VecOld.push_back(u);
  GetOriginalSolver()->HNAverage(u);
}

/**********************************************************/

void MeshInterpolator::AddVectorNew(GlobalVector u)
{
  assert(GetSolverPointer());
  _VecNew.push_back(u);
  GetSolver()->HNAverage(u);
}

/**********************************************************/

void MeshInterpolator::BasicInit(SolverInterface* SI, MeshAgentInterface* MA, const string& name, const ProblemDescriptorInterface* PD)
{
  // Klassenvariablen initialisieren
  _name = name;
  _BaseCells.clear();
  _ToBeRef.clear();
  _ToBeRefNew.clear();
  _NewNodeNumber.clear();
  _VecInt.clear();
  _VecOld.clear();
  _VecNew.clear();

  // Original-Solver und -MeshAgent speichern
  GetOriginalMeshAgentPointer() = MA;
  _dim = GetOriginalMeshAgent()->GetMesh(0)->dimension();
  GetOriginalSolverPointer() = SI;

  // neuen MeshAgent anlegen
  GetMeshAgentPointer() = new MeshAgent;

  GetMeshAgent()->GetShapes2d() = GetOriginalMeshAgent()->GetShapes2d();
  GetMeshAgent()->GetShapes3d() = GetOriginalMeshAgent()->GetShapes3d();
  GetMeshAgent()->SetDefaultValues(_dim,_name+".gup",0);
  GetMeshAgent()->BasicInit(NULL);

  // ProblemDescriptor speichern
  GetProblemPointer() = PD;

  // neuen Solver anlegen
  GetSolverPointer() = new StdSolver;
  GetSolver()->BasicInit(GetMeshAgent()->nlevels()-1,GetOriginalSolver()->GetParamfile(),GetMeshAgent()->GetMesh(0));
  GetSolver()->RegisterVector(_help);
  GetSolver()->NewMesh(GetMeshAgent()->nlevels()-1,GetMeshAgent()->GetMesh(0));
  GetSolver()->SetProblem(*GetProblem());
  GetSolver()->ReInitVector();
  GetSolver()->SetDistribute(false);

  const Q1* Q1DP = dynamic_cast<const Q1*>(GetSolver()->GetDiscretization());
  if (Q1DP)
  {
    int size = static_cast<int>(pow(2.,static_cast<double>(_dim)));
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
    if (_dim==3)
    {
      swap(_weights(i,6),_weights(i,7));
    }
  }
  for (int j=0; j<_weights.n(); j++)
  {
    swap(_weights(2,j),_weights(3,j));
    if (_dim==3)
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

void MeshInterpolator::RhsForProjection(BasicGhostVector& gf)
{
  GlobalVector u;
  ReadBackUpResize(u,_name+".bup");
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
  GetSolver()->NewMesh(GetMeshAgent()->nlevels()-1,GetMeshAgent()->GetMesh(0));
  GetSolver()->ReInitVector();

  GetSolver()->AddNodeVector("U",&u);
  GetSolver()->Rhs(_help);
  GetSolver()->DeleteNodeVector("U");

  AddVectorIntermediate(GetSolver()->GetGV(_help));
  for (set<int>::const_iterator pbc = _BaseCells.begin(); pbc!=_BaseCells.end(); pbc++)
  {
    Distribute(*pbc,*pbc);
  }

  GlobalVector& f = GetOriginalSolver()->GetGV(gf);
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
  GetOriginalSolver()->HNDistribute(gf);

  _Old = NULL;
  _New = NULL;
}

/**********************************************************/
}

