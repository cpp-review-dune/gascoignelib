#define __STL_NO_DRAND48
#include  "hierarchicalmesh.h"
#undef  __STL_NO_DRAND48

#include  "filescanner.h"
#include  "giota.h"
#include  "grandom.h"
#include  "stringutil.h"

using namespace std;

/*------------------------------------------------------*/

HierarchicalMesh::HierarchicalMesh() 
  : mnlevels(1), pdepth(0), curvedshapes(0), withfaces(1), etapatcher(1)
 {}

/*------------------------------------------------------*/

HierarchicalMesh::HierarchicalMesh(const HierarchicalMesh& H)
  : mnlevels(1), pdepth(0), curvedshapes(0), withfaces(1), etapatcher(1)
{
  *this = H;
}

/*------------------------------------------------------*/

HierarchicalMesh::~HierarchicalMesh()
{
  if (curvedshapes!=NULL) delete curvedshapes; curvedshapes=NULL;
}

/*------------------------------------------------------*/

HierarchicalMesh& HierarchicalMesh::operator=(const HierarchicalMesh& H)
{
  // copy nearly all data
  mnlevels  = H.nlevels()-1;

  vo2n.memory(H.Vertexo2n().size());
  eo2n.memory(H.Edgeo2n().size());
  co2n.memory(H.Cello2n().size());

  vo2n      = H.Vertexo2n();
  eo2n      = H.Edgeo2n();
  co2n      = H.Cello2n();

  
  edges     = H.edge();
  pdepth    = H.patchdepth();
  withfaces = H.withfaces;
}

/*------------------------------------------------------*/

void HierarchicalMesh::NewCurvedShapes()
{
  if(!curvedshapes) curvedshapes = new CurvedShapes;
}

/*------------------------------------------------------*/

void HierarchicalMesh::InitShapes(const std::set<std::vector<std::string> >&  curved)
{
  NewCurvedShapes();

  std::set<std::vector<std::string> >::const_iterator  p=curved.begin();

  for (;p!=curved.end(); p++)
    {
      curvedshapes->BasicInit(*p);
    }
}

/*------------------------------------------------------*/

void HierarchicalMesh::clear_transfer_lists()
{
  vo2n.resize(0);
  eo2n.resize(0);
  co2n.resize(0);
}

/*------------------------------------------------------*/

void HierarchicalMesh::SetParameters(string gridname, int patchdepth, int epatcher)
{
  pdepth = patchdepth;
  etapatcher = epatcher;

  if(gridname=="none")
    {
      cerr << "no \"gridname\" " << endl;
      return;
      abort();
    }

  vector<string> s = StringSplit(gridname.c_str(),'.');
  string suff = s[s.size()-1];

  if(suff=="inp")
    {
      read_inp(gridname);
    }
  else if(suff=="gup")
    {
      read_gup(gridname);
    }
  else
    {
      cerr << "HierarchicalMesh::read():\tunknown suffix " << suff << endl;
      abort();
    }
}

/*------------------------------------------------------*/

void HierarchicalMesh::ReadParameters(const ParamFile* pf)
{
  int patchdepth,epatcher;
  int prerefine;
  string gridname;
  set<vector<string> >  curved;

  DataFormatHandler DFH;

  DFH.insert("gridname" ,&gridname,"none");      // inp oder gup Format 
  DFH.insert("prerefine",&prerefine,0);
  DFH.insert("curvedboundary",&curved);
  DFH.insert("patchdepth",&patchdepth,1);
  DFH.insert("etapatcher",&epatcher,1);

  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(pf, "Mesh");

  InitShapes(curved);
  SetParameters(gridname,patchdepth,epatcher);
  global_refine(prerefine);
}

/*---------------------------------------------------*/

void HierarchicalMesh::global_refine(int k)
{
  IntVector cell_coarse(0);

  for (int i=0; i<k; i++)
    {
      IntVector v(ncells());
      iota(v.begin(),v.end(),0);
      refine(v,cell_coarse);
    }
}

/*---------------------------------------------------*/

void HierarchicalMesh::random_refine(double p, int c)
{
  int nq = ncells();
  int nc = 1+(int) (p*nq);      
  
  if (p<0) nc=0;

  IntVector cell_ref(nc), cell_coarse;

  IntVector v(ncells());

  iota(v.begin(),v.end(),0);

  random_sample(v.begin(),v.end(),cell_ref.begin(),cell_ref.end());

  if (c)
    {
      int nq = ncells();
      cell_coarse.resize(nq);
      iota(cell_coarse.begin(),cell_coarse.end(),0);
    }
  refine(cell_ref,cell_coarse);
}

/*---------------------------------------------------*/

void HierarchicalMesh::random_patch_refine(double p, int c)
{
  int nq = ncells();
  int nc = 1+(int) (p*nq);      

  nc = GascoigneMath::min_int(nc,nq);
  if (p<0) nc=0;

  IntVector cell_ref(nc), cell_coarse;
  IntVector v(ncells());
  iota(v.begin(),v.end(),0);
  random_sample(v.begin(),v.end(),cell_ref.begin(),cell_ref.end());

  if (c)
    {
      int nq = ncells();
      cell_coarse.resize(nq);
      iota(cell_coarse.begin(),cell_coarse.end(),0);
    }
  patch_refine(cell_ref,cell_coarse);
}


/*---------------------------------------------------*/

void HierarchicalMesh::vertex_patch_refine(IntVector& refnodes, IntVector& coarsenodes)
{
  IntVector ref, coarse, vertexlevel;
  IntSet    refcoarsenodes, coarsecoarsenodes;

  FillVertexLevels  (vertexlevel);

  // RefineCoarseNodes verteilt auf die Patche
  // can be switched off for more localized refinement
  // (ausgefranzt)
  if (etapatcher)
    {
      RefineCoarseNodes (refcoarsenodes,refnodes,vertexlevel);
    }
  else
    {
      refcoarsenodes.insert(refnodes.begin(),refnodes.end());
    }
  coarsecoarsenodes.insert(coarsenodes.begin(),coarsenodes.end());

  VertexToCells          (ref   ,refcoarsenodes   ,vertexlevel);
  VertexToCellsCoarsening(coarse,coarsecoarsenodes,vertexlevel);

  patch_refine(ref,coarse);
}

/*---------------------------------------------------*/

void HierarchicalMesh::vertex_patch_refine(IntVector& refnodes)
{
  IntVector ref, coarse, vertexlevel;
  IntSet    refcoarsenodes;

  FillVertexLevels  (vertexlevel);

  // RefineCoarseNodes verteilt auf die Patche
  // can be switched off for more localized refinement
  // (ausgefranzt)
  if (etapatcher)
    {
      RefineCoarseNodes (refcoarsenodes,refnodes,vertexlevel);
    }
  else
    {
      refcoarsenodes.insert(refnodes.begin(),refnodes.end());
    }
  VertexToCells     (ref,refcoarsenodes,vertexlevel);

  patch_refine(ref,coarse);
}

/*---------------------------------------------------*/

void HierarchicalMesh::update_edges(IntVector& SwappedEdge)
{
  for (int i=0; i<edges.size(); i++)
    {
      int m = edges[i].master();
      int s = edges[i].slave();

      assert(m>=0);

      int nm = co2n[m];
      if (nm>=0)
	{
	  edges[i].master() = nm;
	  if (s>=0)
	    {
	      edges[i].slave() = co2n[s];
	    }
	}
      else
	{
	  if (s<0)
	    {
	      edges[i].master() = -1;
	      SwappedEdge.push_back(i);
	      continue;
	    }
	  edges[i].swapping(co2n[s]);
	}
    }
}
