#include  "meshagent.h"
#include  "visualization.h"
#include  "filescanner.h"
#include  "gascoignemeshconstructor.h"

/*-----------------------------------------*/

MeshAgent::MeshAgent() : MeshAgentInterface(), HMP(NULL), GMG(NULL), dimension(-1)
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
  GMG->ReInit(dimension,HMP->nlevels());

  GascoigneMeshConstructor MGM(HMP,GMG);
  MGM.BasicInit();
}

/*-----------------------------------------*/

void MeshAgent::BasicInit(const ParamFile* paramfile)
{
  DataFormatHandler DFH;
  DFH.insert("dimension",&dimension,2);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(paramfile,"Mesh");

  if (dimension==2)
    {
      HMP = new HierarchicalMesh2d(paramfile);
    }
  else if (dimension==3)
    {
      HMP = new HierarchicalMesh3d(paramfile);
    }
  else
    {
      std::cout << "dimension of Mesh ? " << dimension << std::endl;
    }
  assert(HMP);

  GMG = NewMultiGridMesh();

  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::read_gup(const std::string& fname)
{
  assert(HMP);
  HMP->read_gup(fname);
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::write_gup(const std::string& fname) const
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
