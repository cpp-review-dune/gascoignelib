#include  "meshagent.h"
#include  "visualization.h"
#include  "filescanner.h"
#include  "gascoignemeshconstructor.h"
#include  "stringutil.h"


using namespace std;
using namespace Gascoigne;

/*-----------------------------------------*/

MeshAgent::MeshAgent() : MeshAgentInterface(), HMP(NULL), GMG(NULL), _dimension(-1)
{
}

/*-----------------------------------------*/

MeshAgent::~MeshAgent()
{
  if (HMP!=NULL) { delete HMP; HMP=NULL;}
  if (GMG!=NULL) { delete GMG; GMG=NULL;}
}

/*-----------------------------------------*/

void MeshAgent::ReInit()
{
  GMG->ReInit(_dimension,HMP->nlevels());

  GascoigneMeshConstructor MGM(HMP,GMG);
  MGM.BasicInit();
}

/*-----------------------------------------*/

void MeshAgent::ReadMesh(int dim, string meshname, int prerefine)
{
  vector<string> s = StringSplit(meshname.c_str(),'.');
  string suff = s[s.size()-1];

  if(suff=="inp")
    {
      HMP->read_inp(meshname);
    }
  else if(suff=="gup")
    {
      HMP->read_gup(meshname);
    }
  else
    {
      cerr << "MeshAgent::read():\tunknown suffix " << suff << endl;
      assert(0);
    }
  HMP->global_refine(prerefine);
}

/*-----------------------------------------*/

void MeshAgent::BasicInit(int dim, string meshname, int prerefine)
{
  _dimension = dim;
  if (_dimension==2)
    {
      HMP = new HierarchicalMesh2d;
    }
  else if (_dimension==3)
    {
      HMP = new HierarchicalMesh3d;
    }
  else
    {
      cout << "dimension of Mesh ? " << _dimension << endl;
    }
  assert(HMP);

  ReadMesh(dim, meshname, prerefine);

  GMG = NewMultiGridMesh();

  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::BasicInit(string meshname, int prerefine, std::map<int,BoundaryFunction<2>* >& curved)
{
  _dimension = 2;
  HMP = new HierarchicalMesh2d;

  std::map<int,BoundaryFunction<2>* >::iterator p=curved.begin();
  for(p=curved.begin();p!=curved.end();p++)
    {
      HMP->AddShape(p->first,p->second);
    }

  ReadMesh(_dimension, meshname, prerefine);

  GMG = NewMultiGridMesh();

  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::BasicInit(string meshname, int prerefine, std::map<int,BoundaryFunction<3>* >& curved)
{
  _dimension = 3;
  HMP = new HierarchicalMesh3d;

  std::map<int,BoundaryFunction<3>* >::iterator p=curved.begin();
  for(p=curved.begin();p!=curved.end();p++)
    {
//       HMP->AddShape(p->first,p->second);
    }

  ReadMesh(_dimension, meshname, prerefine);

  GMG = NewMultiGridMesh();

  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::BasicInit(const ParamFile* paramfile)
{
  DataFormatHandler DFH;
  DFH.insert("dimension",&_dimension,2);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(paramfile,"Mesh");

  if (_dimension==2)
    {
      HMP = new HierarchicalMesh2d(paramfile);
    }
  else if (_dimension==3)
    {
      HMP = new HierarchicalMesh3d(paramfile);
    }
  else
    {
      cout << "dimension of Mesh ? " << _dimension << endl;
    }
  assert(HMP);

  GMG = NewMultiGridMesh();

  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::read_gup(const string& fname)
{
  assert(HMP);
  HMP->read_gup(fname);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::write_gup(const string& fname) const
{
  assert(HMP);
  HMP->write_gup(fname);
}

/*-----------------------------------------*/

void MeshAgent::global_refine(int n)
{
  assert(HMP);
  HMP->global_refine(n);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::random_patch_refine(double p, int n)
{
  assert(HMP);
  HMP->random_patch_refine(p,n);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::refine_nodes(nvector<int>& refnodes)
{
  nvector<int> coarsenodes(0);
  refine_nodes(refnodes,coarsenodes);
}

/*-----------------------------------------*/

void MeshAgent::refine_nodes(nvector<int>& refnodes, nvector<int>& coarsenodes)
{
  assert(HMP);
  HMP->vertex_patch_refine(refnodes,coarsenodes);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::refine_cells(nvector<int>& ref)
{
  nvector<int> refnodes;
  
  for (int i=0; i<ref.size(); i++)
    {
      int cell = ref[i];
      for (int j=0; j<HMP->nodes_per_cell(cell); j++)
	{
	  refnodes.push_back(HMP->vertex_of_cell(cell,j));
	}
    }
  refine_nodes(refnodes);
}
