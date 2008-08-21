#include  "meshagent.h"
#include  "visualization.h"
#include  "filescanner.h"
#include  "gascoignemeshconstructor.h"
#include  "stringutil.h"
#ifdef __NEWER_THAN_GCC_4_2__
#include <tr1/unordered_map>
#define HASHMAP std::tr1::unordered_map
#else
#include  <ext/hash_map>
#define HASHMAP __gnu_cxx::hash_map
#endif

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  MeshAgent::MeshAgent() : MeshAgentInterface(), _goc2nc(false), HMP(NULL), GMG(NULL)
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
  GMG->ReInit(GetDimension(),HMP->nlevels()-HMP->patchdepth());

  GascoigneMeshConstructor MGM(HMP,GMG);
  MGM.BasicInit();
  _celll2g = MGM.Celll2g();
  _cellg2l = MGM.Cellg2l();
  if(HMP->patchdepth()>=2)
  {
    BuildQ4PatchList(MGM.Patchl2g());
    assert(_q4patch.size()==_q4toq2.size());
    PatchIndexHandler &PIH = GMesh(0)->GetPatchIndexHandler();
    nvector<IntVector> &q4patch2cell = PIH.GetAllQ4Patch2Cell();
    q4patch2cell.clear();
    for(int i=0; i<_q4patch.size(); i++)
    {
      IntVector q2p = _q4toq2[i];
      IntVector q4p2c(0);
      for(int j=0; j<q2p.size(); j++)
      {
        int p = q2p[j];
        IntVector cells = PIH.GetPatch2Cell(p);
        q4p2c.insert(q4p2c.end(),cells.begin(),cells.end());
      }
      q4patch2cell.push_back(q4p2c);
    }
    PIH.GetIndexQ4() = _q4patch;
  }

  if(_goc2nc)
  {
      _co2n.clear();
      for(int i =0 ; i < _cl2g.size(); i++)
      {
	  int cn_i = _cl2g[i];//HM Num nach alter nummerierung
	  cn_i = HMP->Cello2n(cn_i);//HM Num nach neuer nummerierung
	  if(cn_i <0)
	  {
	      //Hier wurde vergroebert
	  }
	  else
	  {
	      set<int> kinder;
	      // Zuordnung alte GM Nummern zu neuen
	      if(HMP->sleep(cn_i))
	      {
		  //Zelle verfeinert
		  for(int j = 0; j <HMP->nchilds(cn_i); j++)
		  {
		      kinder.insert(_cellg2l[HMP->child(cn_i,j)]);
		  }
	      }
	      else
	      {
		  kinder.insert(_cellg2l[cn_i]);
	      }
	      _co2n[i] = kinder;
	  }
      }
      

      // Die Var fuer die alten Werteumschreiben
      _cl2g = _celll2g;
      _cg2l = _cellg2l;
      
      //den fathers Vektor neu fuellen
      _fathers.clear();
      _fathers.resize(_cl2g.size());
      for(int i = 0; i < _fathers.size(); i++)
      {
	  _fathers[i] = HMP->Vater(_cl2g[i]);
      }
  }
}

/*-----------------------------------------*/

void MeshAgent::BuildQ4PatchList(const IntVector &patchl2g)
{
  _q4patch.resize(0);
  _q4toq2.resize(0);
  IntSet q2patch,q4patchset;
  HMP->GetAwakePatchs(q2patch);
  for(IntSet::const_iterator p=q2patch.begin(); p!=q2patch.end(); p++)
  {
    assert(*p!=-1);
    int vater = HMP->Vater(*p);
    assert(vater!=-1);
    q4patchset.insert(vater);
  }
  for(IntSet::const_iterator p=q4patchset.begin(); p!=q4patchset.end(); p++)
  {
    _q4patch.push_back(HMP->ConstructQ4Patch(*p));
  }
  if(patchl2g.size()!=q2patch.size())
  {
    cerr << "MeshAgent::BuildQ4PatchList: patchl2g must be same size as q2patch!!!" << endl;
    abort();
  }
  HASHMAP<int,int> patchg2l;
  for(int i=0; i<patchl2g.size(); i++)
  {
    //Edit patchg2l[patchl2g[i]] = i;
    patchg2l.insert(make_pair<int,int>(patchl2g[i],i));
  }
  IntVector perm(4);
  perm[0]=0;
  perm[1]=1;
  perm[2]=3;
  perm[3]=2;
  if(HMP->dimension()==3)
  {
    perm.resize(8);
    perm[4]=4;
    perm[5]=5;
    perm[6]=7;
    perm[7]=6;
  }
  _q4toq2.resize(q4patchset.size());
  int q4_l = 0;
  for(IntSet::const_iterator p=q4patchset.begin(); p!=q4patchset.end(); p++,q4_l++)
  {
    const IntVector &K = HMP->Kinder(*p);
    assert(K.size()==perm.size());
    for(int i=0; i<K.size(); i++)
    {
      assert(patchg2l.find(K[perm[i]])!=patchg2l.end());
      _q4toq2[q4_l].push_back(patchg2l[K[perm[i]]]);
    }
  }
}

/*-----------------------------------------*/

void MeshAgent::BasicInit(const ParamFile* paramfile)
{
  assert(HMP==NULL);
  int dim = 0;

  DataFormatHandler DFH;
  DFH.insert("dimension",&dim);
  //um die zuordnung alte GMNr. -> GMNr. an/abzuschalten
  DFH.insert("cellnumtrans",&_goc2nc,false);
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(paramfile,"Mesh");

  if (dim==2)
    {
      HMP = new HierarchicalMesh2d;
      for(map<int,BoundaryFunction<2>* >::const_iterator p=_curved2d.begin();p!=_curved2d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
    }
  else if (dim==3)
    {
      HMP = new HierarchicalMesh3d;
      for(map<int,BoundaryFunction<3>* >::const_iterator p=_curved3d.begin();p!=_curved3d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
    }
  else
    {
      cout << "dimension of Mesh ? " << dim << endl;
    }
  assert(HMP);
  HMP->BasicInit(paramfile);
  
  GMG = NewMultiGridMesh();

  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::BasicInit(const string& gridname, int dim, int patchdepth, int epatcher, bool goc2nc)
{
  assert(HMP==NULL);
  _goc2nc = goc2nc;
  if (dim==2)
    {
      HMP = new HierarchicalMesh2d;
      for(map<int,BoundaryFunction<2>* >::const_iterator p=_curved2d.begin();p!=_curved2d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
    }
  else if (dim==3)
    {
      HMP = new HierarchicalMesh3d;
      for(map<int,BoundaryFunction<3>* >::const_iterator p=_curved3d.begin();p!=_curved3d.end();p++)
        {
          HMP->AddShape(p->first,p->second);
        }
    }
  else
    {
      cout << "dimension of Mesh ? " << dim << endl;
    }
  assert(HMP);
  HMP->SetParameters(gridname,patchdepth,epatcher);
  
  GMG = NewMultiGridMesh();

  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::read_gup(const string& fname)
{
  assert(HMP);
  HMP->read_gup(fname);
  ClearCl2g();
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::read_gip(const string& fname)
{
  assert(HMP);
  HMP->read_gip(fname);
  ClearCl2g();
  ReInit();
}

/*-----------------------------------------*/

void MeshAgent::write_gup(const string& fname) const
{
  assert(HMP);
  HMP->write_gup(fname);
}

/*-----------------------------------------*/

void MeshAgent::write_gip(const string& fname) const
{
  assert(HMP);
  HMP->write_gip(fname);
}

/*-----------------------------------------*/

void MeshAgent::write_inp(const string& fname) const
{
  assert(HMP);
  HMP->write_inp(fname);
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

/*----------------------------------------*/

inline const set<int> MeshAgent::Cello2n(int i)const
{
    map<int,set<int> >::const_iterator p = _co2n.find(i);
    if(p == _co2n.end())
    {
	return set<int>();
    }
    else
    {
	return p->second;
    }
}

/*----------------------------------------*/

inline const int MeshAgent::Cello2nFather(int i)const
{
    assert(_co2n.find(i)==_co2n.end());
    //Umrechnung alte HM nummer in neue GM nummer
    return _cg2l.find(HMP->Cello2n(_fathers[i]))->second;
}


}

#undef HASHMAP
