#include  "levelmesh3d.h"
#include  "nmatrix.h"
#include  "levelsorter3d.h"
#include  "leveljumper.h"
#include  "set2vec.h"


using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

LevelMesh3d::LevelMesh3d(const HierarchicalMesh* hmp) : 
  Index()   
{
  HMP = dynamic_cast<const HierarchicalMesh3d*>(hmp);
}

/*---------------------------------------------------*/

LevelMesh3d::~LevelMesh3d() 
{}

/*---------------------------------------------------*/

void LevelMesh3d::BasicInit(const set<int>& newh, const set<int>& oldh)
{
  int n = newh.size()+oldh.size();
  Index::Hexl2g().memory(n);
  IntVector::const_iterator p = set_union(newh.begin(),newh.end(),oldh.begin(),oldh.end(),
			   Index::Hexl2g().begin());
  n = p-Index::Hexl2g().begin();

  InitCells(n);
  InitNodes(n);
  InitEdges(n);
}

/*-----------------------------------------*/

void LevelMesh3d::InitCells(int n)
{
  Index::Hexl2g().memory(n);

  sort(Index::Hexl2g().begin(), Index::Hexl2g().end(), LevelSorter3d(*HMP));

  Index::InitHexs();
}

/*-----------------------------------------*/

void LevelMesh3d::InitNodes(int n)
{
  IntSet nodes;
  for(int i=0;i<n;i++)
    {
      int ind = Index::Hexl2g()[i];
      for(int ii=0;ii<8;ii++)  
	{
	  nodes.insert(HMP->vertex_of_cell(ind,ii));
	}
    }
  Index::InitNodes(nodes);
}

/*-----------------------------------------*/

void LevelMesh3d::InitEdges(int n)
{
  IntSet edges;
  for(int i=0;i<n;i++)
    {
      const Hex& Q = HMP->hex(Index::Hexl2g()[i]);
      for(int ii=0;ii<6;ii++)  
	{
	  edges.insert(Q.edge(ii));
	}
    }
  Index::InitEdges(edges);

  //stable_sort(Edgel2g().begin(),Edgel2g().end(),HangEdgeSort5(*this));

  Edgeg2l().clear();
  for(int i=0;i<Edgel2g().size();i++)
    {
      Edgeg2l().insert(make_pair(Edgel2g()[i],i));
    }
}

/*-----------------------------------------*/

bool LevelMesh3d::BuildFathers(set<int>&  Vaeter) const
{
  for(int i=0; i<ncells();i++)
    {
      const Hex& h = hex(i);
      int findex = h.father();
      if(findex==-1) 
	{
	  return 0;
	}

      const Hex& hf = HMP->hex(findex);
      for(int ii=0;ii<hf.nchilds();ii++) {
	int cindex = hf.child(ii);
	if(Hexg2lCheck(cindex)==-2) {
	  return 0;
	}
      }
      Vaeter.insert(findex);
    }
  return 1;
}

/*---------------------------------------------------*/

void LevelMesh3d::ConstructIndOfPatch(nvector<IntVector>& dst) const
{
  set<int>  Vaeter;
  BuildFathers(Vaeter);

  int nh = ncells()/8;
  dst.reserve(nh);
  dst.resize (nh,IntVector(27));

  int count=0;
  set<int>::const_iterator pf = Vaeter.begin();
  nmatrix<int> A(27,2);

  A(0,0) = 0; A(0,1) = 0;
  A(1,0) = 0; A(1,1) = 1;
  A(2,0) = 1; A(2,1) = 1;
  A(3,0) = 0; A(3,1) = 3;
  A(4,0) = 0; A(4,1) = 2;
  A(5,0) = 1; A(5,1) = 2;
  A(6,0) = 3; A(6,1) = 3;
  A(7,0) = 3; A(7,1) = 2;
  A(8,0) = 2; A(8,1) = 2;
  
  A(9 ,0) = 4; A(9 ,1) = 0;
  A(10,0) = 4; A(10,1) = 1;
  A(11,0) = 5; A(11,1) = 1;
  A(12,0) = 4; A(12,1) = 3;
  A(13,0) = 4; A(13,1) = 2;
  A(14,0) = 5; A(14,1) = 2;
  A(15,0) = 7; A(15,1) = 3;
  A(16,0) = 7; A(16,1) = 2;
  A(17,0) = 6; A(17,1) = 2;

  A(18,0) = 4; A(18,1) = 4;
  A(19,0) = 4; A(19,1) = 5;
  A(20,0) = 5; A(20,1) = 5;
  A(21,0) = 4; A(21,1) = 7;
  A(22,0) = 4; A(22,1) = 6;
  A(23,0) = 5; A(23,1) = 6;
  A(24,0) = 7; A(24,1) = 7;
  A(25,0) = 7; A(25,1) = 6;
  A(26,0) = 6; A(26,1) = 6;

  while(pf!=Vaeter.end())
    {
      int findex = *pf;
      const Hex& hf = HMP->hex(findex);
      fixarray<8,int> FineHexs;
      for(int ii=0;ii<hf.nchilds();ii++) 
	{
	  FineHexs[ii] = hf.child(ii);
	}
      for (int i=0; i<27; i++)
	{
	  int fh = FineHexs[A(i,0)];
	  int gi = Vertexg2l( HMP->hex(fh).vertex(A(i,1)) );
	  dst[count][i] = gi;
	}
      count++;
      pf++;
    }
}

/*---------------------------------------------------*/

bool LevelMesh3d::ConstructCellIndOfPatch(nvector<int>& dst) const
{
  set<int>  Vaeter;
  BuildFathers(Vaeter);

  int nh = ncells()/8;
  dst.reservesize(nh);

  int count=0;
  set<int>::const_iterator pf = Vaeter.begin();

  while(pf!=Vaeter.end())
    {
      int findex = *pf;
      dst[count] = findex;

      count++;
      pf++;
    }
  return 1;
}

/*---------------------------------------------------*/

void LevelMesh3d::ConstructHangingStructureQuadratic(QuadraticHNStructure3& hnq2,
						     QuadraticHNStructure9& hnq2face) const
{
  hnq2.clear();
  hnq2face.clear();
  set<int> habschon; 

  const HexLawAndOrder& HLaO = HMP->HexLawOrder();

  for (int i=0;i<ncells();i++)
    {
      const Hex& q = HMP->hex(Hexl2g(i));
      int father = q.father();
      if (father<0) continue;
      if (habschon.find(father)!=habschon.end()) continue;
      habschon.insert(father);
      for(int in=0; in<6; in++)
	{
	  int neighbour = HMP->neighbour(father,in);
	  if (neighbour<0) continue;
	  const Hex& qfn = HMP->hex(neighbour);
	  if (qfn.nchilds()==0) continue;

	  int ne = in;
	  {
	    int start = 0;
	    int neighbourneighbour = HMP->neighbour(neighbour,ne);
	    while ((neighbourneighbour!=father) && (start<6))
	      {
		start++;
		ne = (ne+1)%6;
		neighbourneighbour = HMP->neighbour(neighbour,ne);
	      }
	    assert(neighbourneighbour==father);
	  }

	  fixarray<4,int> childs;
	  HLaO.childs_of_face(childs,qfn,ne);
	  
	  {
	    int child = childs[0];
	    if (Hexg2lCheck(child)>=0) continue;

	    const Hex& qfc = HMP->hex(child);
	    
	    if (qfc.nchilds()==0) continue;
	    
	    if (Hexg2lCheck(qfc.child(0))<0) continue;

	  }

	  // jetzt haengt er

	  fixarray<9,int> F = HLaO.PatchVerticesOfFace(neighbour,ne);
	  int nec = HLaO.ChildFace(ne);

	  // ordne F;
	  for (int j=0; j<4; j++)
	    {
	      const Hex& qfcc = HMP->hex(childs[j]);

	      int hnl = Vertexg2l(HLaO.face_vertex(qfcc,nec));
	      
	      // permutiere F
	      fixarray<9,int> G = HLaO.GiveOrdering(F,qfcc);
	  
	      fixarray<4,int> face;
	      face[0] = F[G[0]];
	      face[1] = F[G[1]];
	      face[2] = F[G[3]];
	      face[3] = F[G[4]];

	      for (int k=0; k<9; k++)  G[k] = Vertexg2l(F[G[k]]);
	      
	      hnq2face.insert(make_pair(hnl,G));

	      fixarray<4,int> childface;
	      HLaO.GetFace(childface,childs[j],nec);
	      // nun die hanging lines
	      for (int e=0; e<4; e++)
		{
		  int hne = Vertexg2l(HLaO.EdgeVertexOfFace(qfcc,face,e));
		  if (hnq2.find(hne)==hnq2.end())
		    {
		      fixarray<3,int> line;
		      int e0 = childface[e];
		      int e1 = childface[(e+1)%4];

		      if (e0==F[4]) swap(e1,e0);
		      if (e1!=F[4])
			{
			  if ((e0==F[1]) || (e0==F[3]) || (e0==F[5]) || (e0==F[7]))
			    {
			      swap(e1,e0);
			    }
			}

		      line[0] = Vertexg2l(e0);
		      line[1] = Vertexg2l(e1);

		      int last;
		      if      ( (e0==F[0]) && (e1==F[1])) last = F[2];
		      else if ( (e0==F[2]) && (e1==F[1])) last = F[0];
		      else if ( (e0==F[3]) && (e1==F[4])) last = F[5];
		      else if ( (e0==F[5]) && (e1==F[4])) last = F[3];
		      else if ( (e0==F[6]) && (e1==F[7])) last = F[8];
		      else if ( (e0==F[8]) && (e1==F[7])) last = F[6];
		      else if ( (e0==F[0]) && (e1==F[3])) last = F[6];
		      else if ( (e0==F[6]) && (e1==F[3])) last = F[0];
		      else if ( (e0==F[1]) && (e1==F[4])) last = F[7];
		      else if ( (e0==F[7]) && (e1==F[4])) last = F[1];
		      else if ( (e0==F[2]) && (e1==F[5])) last = F[8];
		      else if ( (e0==F[8]) && (e1==F[5])) last = F[2];
		      else  assert(0);

		      line[2] = Vertexg2l(last);
		      hnq2.insert(make_pair(hne,line));
		    }
		}
	    }
	}
    }
  //  cout << "****** LevelM 3d " << hnq2face.size() << " " << hnq2 .size() << endl;
}

/*---------------------------------------------------*/

bool LevelMesh3d::EnkelUniform(const Hex& Q) const
{
  for (int ii=0; ii<Q.nchilds(); ii++)
    {
      int qccindex = Q.child(ii);
      const Hex& qcc = HMP->hex(qccindex);
      for(int iii=0;iii<qcc.nchilds();iii++)
	{
	  int qcindex = qcc.child(iii);
	  if(Hexg2lCheck(qcindex)==-2)
	    {
	      return 0;
	    }
	}
    }
  return 1;
}

/*---------------------------------------------------*/

void LevelMesh3d::fill_opis(IntSet& dst, IntSet& oldhexs) const
{
  for(int i=0; i<ncells(); i++)
    {
      const Hex& Q = hex(i);

      int f = Q.father();
      assert(f>=0);

      int opi = HMP->hex(f).father();

      if (opi<0) 
	{
	  int j = Hexl2g(i);
	  oldhexs.insert(j);
	}
      else
	{
	  dst.insert(opi);
	}
    }
}

/*---------------------------------------------------*/

void LevelMesh3d::fill_childs(IntSet& dst, const Hex& Q) const
{
  for (int i=0; i<Q.nchilds(); i++)
    {
      dst.insert(Q.child(i));
    }
}

/*---------------------------------------------------*/

void LevelMesh3d::fill_enkel(IntSet& dst, const Hex& Q) const
{
  for (int i=0; i<Q.nchilds(); i++)
    {
      const Hex& C = HMP->hex(Q.child(i));
      for (int j=0; j<C.nchilds(); j++)
	{
	  int cc = C.child(j);
	  if (Hexg2lCheck(cc)>=0)
	    {
	      dst.insert(cc);
	    }
	}
    }
}

/*---------------------------------------------------*/

void LevelMesh3d::check_leveljump() const
{
  LevelJumper  Phi;
  for(int c=0;c<ncells();c++)
    {
      Phi.update(hex(c));
    }
  assert(! Phi.check());
}

/*---------------------------------------------------*/

void LevelMesh3d::construct_lists(IntSet& newhexs, IntSet& oldhexs) const
{
  newhexs.clear();  
  oldhexs.clear();

  check_leveljump();

  set<int>  Opis;
  fill_opis(Opis,oldhexs);
  for(set<int>::const_iterator p=Opis.begin();p!=Opis.end();p++)
    {
      const Hex& Q = HMP->hex(*p);
      
      if (EnkelUniform(Q))
	{
	  fill_childs(newhexs,Q);
	}
      else
	{
	  fill_enkel(oldhexs,Q);
	}
    }

  // Iteration zum Regulaer machen (LevelJump)
  while(1)
    {
      LevelJumper  Phi;
      for(set<int>::const_iterator p=newhexs.begin(); p!=newhexs.end(); p++)
	{
	  Phi.update(HMP->hex(*p));
	}
      for(set<int>::const_iterator p=oldhexs.begin(); p!=oldhexs.end(); p++)
	{
	  Phi.update(HMP->hex(*p));
	}
      if(!Phi.check()) break;

      int rep=0;
      IntSet help(newhexs);
      for(set<int>::const_iterator p=newhexs.begin(); p!=newhexs.end(); p++)
	{
	  const Hex& q = HMP->hex(*p);
	  if (!Phi.VertexOK(q))
	    {
	      rep++;
	      const Hex& qf = HMP->hex(q.father());
	      for(int ii=0;ii<8;ii++)
		{
		  help.erase(qf.child(ii));
		}
	      fill_enkel(oldhexs,qf);
	    }
	}
      newhexs = help;
//       cerr << "\t Regular Iteration\t" << count++ << " " << rep << endl;
    }
}

/*---------------------------------------------------*/

void LevelMesh3d::InitBoundaryHandler(BoundaryIndexHandler& BI) const
{
    // bquads
  // probably slowly (could be in multigridmesh !)
  // since we cannnot go from quad -> bline
  IntSet bquads;
  for(int i=0;i<HMP->nbquads();i++)
    {
      int q = HMP->bquad(i).of_quad();
      if(Hexg2l().find(q)!=Hexg2l().end())
	{
	  bquads.insert(i);
	}
    }

  // which colors are there ?
  BI.clear();
  for(IntSet::const_iterator p=bquads.begin();
      p!=bquads.end(); p++)
    {
      const BoundaryQuad& bl = HMP->bquad(*p);
      int col = bl.material();
      BI.GetColors().insert(col);
    }
  IntVector colorvec;
  Set2Vec(colorvec,BI.GetColors());

  // compute inverse positions
  map<int,int> inv;

  for (int i=0; i<colorvec.size(); i++)
    {
      inv.insert(make_pair(colorvec[i],i));
    }
  
  int nc = colorvec.size(); 
  vector<set<int>  > H1(nc);  // for verteces
  // for cells and local indices
  vector<set<fixarray<2,int> > >  H2(nc); 

  for(IntSet::const_iterator q=bquads.begin();
      q!=bquads.end(); q++)
    {
      const BoundaryQuad& bl = HMP->bquad(*q);
      int col = bl.material();

      map<int,int>::const_iterator p = inv.find(col);
      if (p==inv.end())
	{
	  cout << "BoundaryIndexHandler::init3d"<< endl;
	  abort();
	}
      int pos = p->second;

      for(int ii=0;ii<4;ii++)
	{
	  int vindex = Vertexg2l(bl.vertex(ii));
	  H1[pos].insert(vindex);
	}
      fixarray<2,int> ind;
      ind[0] = Hexg2l(bl.of_quad());
      ind[1] = bl.edge_in_quad();
      H2[pos].insert(ind);
    }
  BI.CopySetToVector(H1,colorvec,BI.GetVertex());
  for (int i=0; i<H2.size(); i++)
    {
      IntVector v1(H2[i].size());
      IntVector v2(H2[i].size());
      int j = 0;
      
      set<fixarray<2,int> >::const_iterator p;
      for (p=H2[i].begin(); p!=H2[i].end(); p++)
	{
	  v1[j] = (*p)[0];
	  v2[j] = (*p)[1];
	  j++;
	}
      int color = colorvec[i];

      BI.GetCell ().insert(make_pair(color,v1));
      BI.GetLocal().insert(make_pair(color,v2));
    }
}
