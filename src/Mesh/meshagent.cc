#include  "meshagent.h"
#include  "visualization.h"
#include  "filescanner.h"
#include  "gascoignemeshconstructor.h"
#include  "stringutil.h"


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
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
  cerr << "*********************************************************" << endl;
  cerr << "Der Aufruf BasicInit(dim,meshname,prerefine) sollte durch" << endl;
  cerr << "    SetDefaultValues(dim,meshname,prerefine)" << endl;
  cerr << "ersetzt werden!!!" << endl;
  cerr << "*********************************************************" << endl;
  assert(0);

  _dimension = dim;
  if (_dimension==2)
    {
      HMP = new HierarchicalMesh2d;
      for(map<int,BoundaryFunction<2>* >::const_iterator p=_curved2d.begin();p!=_curved2d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
    }
  else if (_dimension==3)
    {
      HMP = new HierarchicalMesh3d;
      for(map<int,BoundaryFunction<3>* >::const_iterator p=_curved3d.begin();p!=_curved3d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
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

void MeshAgent::BasicInit(const ParamFile* paramfile)
{
  DataFormatHandler DFH;
  DFH.insert("dimension",&_dimension);
  DFH.insert("gridname" ,&_gridname);      // inp oder gup Format 
  DFH.insert("prerefine",&_prerefine);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(paramfile,"Mesh");

  if (_dimension==2)
    {
      HMP = new HierarchicalMesh2d;
      for(map<int,BoundaryFunction<2>* >::const_iterator p=_curved2d.begin();p!=_curved2d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
    }
  else if (_dimension==3)
    {
      HMP = new HierarchicalMesh3d;
      for(map<int,BoundaryFunction<3>* >::const_iterator p=_curved3d.begin();p!=_curved3d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
    }
  else
    {
      cout << "dimension of Mesh ? " << _dimension << endl;
    }
  assert(HMP);
 
  ReadMesh(_dimension, _gridname, _prerefine);

  GMG = NewMultiGridMesh();

  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::SetDefaultValues(int dim, string gridname, int prerefine)
{
  _dimension = dim;
  _gridname  = gridname;
  _prerefine = prerefine;
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

void MeshAgent::global_patch_coarsen(int n)
{
  assert(HMP);
  HMP->global_patch_coarsen(n);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::global_refine(int n)
{
  assert(HMP);
  HMP->global_refine(n);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::random_patch_coarsen(double p, int n)
{
  assert(HMP);
  HMP->random_patch_coarsen(p,n);
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

void MeshAgent::refine_nodes(IntVector& refnodes)
{
  IntVector coarsenodes(0);
  refine_nodes(refnodes,coarsenodes);
}

/*-----------------------------------------*/

void MeshAgent::refine_nodes(IntVector& refnodes, IntVector& coarsenodes)
{
  assert(HMP);
  HMP->vertex_patch_refine(refnodes,coarsenodes);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::refine_cells(IntVector& ref)
{
  IntVector refnodes;
  
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
}
