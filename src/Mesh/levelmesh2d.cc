#include  "levelmesh2d.h"
#include  "levelsorter.h"
#include  "leveljumper.h"
#include  "set2vec.h"

using namespace std;
 
/*---------------------------------------------------*/

LevelMesh2d::LevelMesh2d(const HierarchicalMesh* hmp) 
  : Index()
{
  HMP = dynamic_cast<const HierarchicalMesh2d*>(hmp);
}

/*---------------------------------------------------*/

LevelMesh2d::~LevelMesh2d() 
{}

/*---------------------------------------------------*/

bool LevelMesh2d::EdgeIsHangingGlobalIndex(int i) const
{
  int igm = HMP->edge(i).master();
  int igs = HMP->edge(i).slave();

  // rand oder kleien hedges
  if(igs<0) return 0;

  bool m_in_lmesh = (Quadg2l().find(igm)!=Quadg2l().end());
  bool s_in_lmesh = (Quadg2l().find(igs)!=Quadg2l().end());

  if(m_in_lmesh && s_in_lmesh) return 0;

  int ivg = HMP->NodeOnEdge(i);
  if(Vertexg2l().find(ivg)==Vertexg2l().end()) return 0;

  return 1;
}

/*---------------------------------------------------*/

void LevelMesh2d::BasicInit(const std::set<int>& newq, const std::set<int>& oldq) 
{
  // doch sortiert

  int n = newq.size()+oldq.size();
  Index::Quadl2g().memory(n);
  IntVector::const_iterator p = set_union(newq.begin(),newq.end(),oldq.begin(),oldq.end(),
					  Index::Quadl2g().begin());
  n = p-Index::Quadl2g().begin();

  InitCells(n);
  InitNodes(n);
  //  InitEdges(n);
}

/*-----------------------------------------*/

void LevelMesh2d::InitCells(int n)
{
  Index::Quadl2g().memory(n);

  std::sort(Index::Quadl2g().begin(), Index::Quadl2g().end(), LevelSorter2d(*HMP));

  Index::InitQuads();
}

/*-----------------------------------------*/

void LevelMesh2d::InitNodes(int n)
{
  IntSet nodes;
  for(int i=0;i<n;i++)
    {
      int ind = Index::Quadl2g()[i];
      for(int ii=0;ii<4;ii++)  
	{
	  nodes.insert(HMP->vertex_of_cell(ind,ii));
	}
    }
  Index::InitNodes(nodes);
}

/*-----------------------------------------*/

void LevelMesh2d::InitEdges(int n)
{
  // edges
  IntSet edges;
  for(int i=0;i<n;i++)
    {
      const Quad& Q = HMP->quad(Index::Quadl2g()[i]);
      for(int ii=0;ii<4;ii++)  
	{
	  edges.insert(Q.edge(ii));
	}
    }
  Index::InitEdges(edges);

  // sortiere haengende edges nach hinten

  stable_sort(Edgel2g().begin(),Edgel2g().end(),HangEdgeSort3(*this));

  Edgeg2l().clear();
  for(int i=0;i<Edgel2g().size();i++)
    {
      Edgeg2l().insert(std::make_pair(Edgel2g()[i],i));
    }
}

/*-----------------------------------------*/

bool LevelMesh2d::BuildFathers(std::set<int>&  Vaeter) const
{
  for(int i=0; i<ncells();i++)
    {
      const Quad& q = quad(i);
      int findex = q.father();
      if(findex==-1) 
	{
	  return 0;
	}

      const Quad& qf = HMP->quad(findex);
      for(int ii=0;ii<qf.nchilds();ii++) {
	int cindex = qf.child(ii);
	if(Quadg2lCheck(cindex)==-2) {
	  return 0;
	}
      }
      Vaeter.insert(findex);
    }
  return 1;
}

/*-----------------------------------------*/

bool LevelMesh2d::ConstructCellIndOfPatch(nvector<int>& dst) const
{
  std::set<int>  Vaeter;
  BuildFathers(Vaeter);

  int nq = ncells()/4;
  dst.reservesize(nq);

  int count=0;
  std::set<int>::const_iterator pf = Vaeter.begin();
  while(pf!=Vaeter.end())
    {
      int findex = *pf;
      const Quad& qf = HMP->quad(findex);

      dst[count] = findex;
     
      count++;
      pf++;
    }
  return 1;
}


/*-----------------------------------------*/

void LevelMesh2d::ConstructIndOfPatch(nvector<IntVector>& dst) const
{
  std::set<int>  Vaeter;
  BuildFathers(Vaeter);

  int nq = ncells()/4;
  dst.reserve(nq);
  dst.resize (nq,IntVector(9));

  int count=0;
  std::set<int>::const_iterator pf = Vaeter.begin();
  while(pf!=Vaeter.end())
    {
      int findex = *pf;
      const Quad& qf = HMP->quad(findex);

      fixarray<4,int> FineQuads;
      for(int ii=0;ii<qf.nchilds();ii++) 
	{
	  FineQuads[ii] = qf.child(ii);
	}
      dst[count][0] = Vertexg2l( HMP->quad(FineQuads[0]).vertex(0) );
      dst[count][1] = Vertexg2l( HMP->quad(FineQuads[0]).vertex(1) );
      dst[count][2] = Vertexg2l( HMP->quad(FineQuads[1]).vertex(1) );
      dst[count][3] = Vertexg2l( HMP->quad(FineQuads[0]).vertex(3) );
      dst[count][4] = Vertexg2l( HMP->quad(FineQuads[0]).vertex(2) );
      dst[count][5] = Vertexg2l( HMP->quad(FineQuads[1]).vertex(2) );
      dst[count][6] = Vertexg2l( HMP->quad(FineQuads[3]).vertex(3) );
      dst[count][7] = Vertexg2l( HMP->quad(FineQuads[3]).vertex(2) );
      dst[count][8] = Vertexg2l( HMP->quad(FineQuads[2]).vertex(2) );
     
      count++;
      pf++;
    }
}

/*---------------------------------------------------*/

void LevelMesh2d::ConstructHangingStructureQuadratic(QuadraticHNStructure3& hnq2) const
{
  hnq2.clear();
  int count=0;
  std::set<int> habschon; 
  for (int i=0;i<ncells();i++)
    {
      const Quad& q = HMP->quad(Quadl2g(i));
      int father = q.father();
      if (father<0) continue;
      if (habschon.find(father)!=habschon.end()) continue;
      habschon.insert(father);
      for(int in=0; in<4; in++)
	{
	  int neighbour = HMP->neighbour(father,in);
	  if (neighbour<0) continue;
	  const Quad& qfn = HMP->quad(neighbour);
	  if (qfn.nchilds()==0) continue;

	  fixarray<2,int> childs;
	  int ne = in;
	  
	  {
	    int start = 0;
	    int neighbourneighbour = HMP->neighbour(neighbour,ne);
	    while ((neighbourneighbour!=father) && (start<4))
	      {
		start++;
		ne = (ne+1)%4;
		neighbourneighbour = HMP->neighbour(neighbour,ne);
	      }
	    
	    assert(neighbourneighbour==father);
	  }	  
	  HMP->QuadLawOrder().childs_of_edge(childs,qfn,ne);
	  int child = childs[0];
	  if (Quadg2lCheck(child)>=0) continue;
	  const Quad& qfc = HMP->quad(child);
	  
	  if (qfc.nchilds()==0) continue;

	  int enkel = qfc.child(0);
	  if (Quadg2lCheck(enkel)<0) continue;
	  
	  // jetzt haengt er
	  int hn = Vertexg2l( HMP->QuadLawOrder().edge_vertex(qfc,ne) );
	  fixarray<3,int> F;
	  F[0] = qfn[ne];
	  F[1] = HMP->QuadLawOrder().edge_vertex(qfn,ne);
	  F[2] = qfn[(ne+1)%4];
	  
	  assert((qfc[0]==F[0]) || (qfc[1]==F[0]) ||
		 (qfc[2]==F[0]) || (qfc[3]==F[0]));
	  
	  for (int k=0; k<3; k++)  F[k] = Vertexg2l(F[k]);
	  
	  hnq2.insert(std::make_pair(hn,F));
	  
	  const Quad& qfc2 = HMP->quad(childs[1]);
	  hn = Vertexg2l(HMP->QuadLawOrder().edge_vertex(qfc2,ne));
	  
	  std::swap(F[0],F[2]);
	  
	  hnq2.insert(std::make_pair(hn,F));
	}
    }
}

/*---------------------------------------------------*/

void LevelMesh2d::check_leveljump() const
{
  LevelJumper  Phi;
  for(int c=0;c<ncells();c++)
    {
      Phi.update(quad(c));
    }
  assert(! Phi.check());
  //if(Phi.check()) std::cerr << "LevelMesh2d::check_leveljump() aenderb\n";
}

/*---------------------------------------------------*/

void LevelMesh2d::fill_opis(IntSet& dst, IntSet& oldquads) const
{
  for(int i=0; i<ncells(); i++)
    {
      const Quad& Q = quad(i);

      int f = Q.father();
      assert(f>=0);

      int opi = HMP->quad(f).father();

      if (opi<0) 
	{
	  int j = Quadl2g(i);
	  oldquads.insert(j);
	}
      else
	{
	  dst.insert(opi);
	}
    }
}

/*---------------------------------------------------*/

void LevelMesh2d::fill_childs(IntSet& dst, const Quad& Q) const
{
  for (int ii=0; ii<Q.nchilds(); ii++)
    {
      int qccindex = Q.child(ii);
      dst.insert(qccindex);
    }
}

/*---------------------------------------------------*/

void LevelMesh2d::fill_enkel(IntSet& oldquads, const Quad& Q) const
{
  for (int ii=0; ii<Q.nchilds(); ii++)
    {
      int qccindex = Q.child(ii);
      const Quad& qcc = HMP->quad(qccindex);
      for(int iii=0;iii<qcc.nchilds();iii++)
	{
	  int qcindex = qcc.child(iii);
	  if(Quadg2lCheck(qcindex)>=0)
	    {
	      oldquads.insert(qcindex);
	    }
	}
    }
}

/*---------------------------------------------------*/

bool LevelMesh2d::EnkelUniform(const Quad& Q) const
{
  bool regular=1;
  
  for (int ii=0; ii<Q.nchilds(); ii++)
    {
      int qccindex = Q.child(ii);
      const Quad& qcc = HMP->quad(qccindex);
      for(int iii=0;iii<qcc.nchilds();iii++)
	{
	  int qcindex = qcc.child(iii);
	  if(Quadg2lCheck(qcindex)==-2)
	    {
	      return 0;
	    }
	}
    }
  return 1;
}

/*---------------------------------------------------*/

void LevelMesh2d::construct_lists(IntSet& newquads, IntSet& oldquads) const
{
  newquads.clear();  
  oldquads.clear();

  check_leveljump();

  std::set<int>     Opis;
  fill_opis(Opis,oldquads);
  for(std::set<int>::const_iterator p=Opis.begin();p!=Opis.end();p++)
    {
      const Quad& Q = HMP->quad(*p);
      
      if ( EnkelUniform(Q) )
	{
	  fill_childs(newquads,Q);
	}
      else
	{
	  fill_enkel(oldquads,Q);
	}
    }

  // Iteration zum Regulaer machen (LevelJump)
  int count = 0;
  while(1)
    {
      LevelJumper  Phi;
      std::set<int>::const_iterator p;
      for(p=newquads.begin(); p!=newquads.end(); p++)
	{
	  Phi.update(HMP->quad(*p));
	}
      for(p=oldquads.begin(); p!=oldquads.end(); p++)
	{
	  Phi.update(HMP->quad(*p));
	}
      if(!Phi.check()) break;

      //int rep=0;
      IntSet help(newquads);
      for(p=newquads.begin(); p!=newquads.end(); p++)
	{
	  const Quad& q = HMP->quad(*p);
	  if (!Phi.VertexOK(q))
	    {
	      //rep++;
	      const Quad& qf = HMP->quad(q.father());
	      for(int ii=0;ii<4;ii++)
		{
		  help.erase(qf.child(ii));
		}
	      fill_enkel(oldquads,qf);
	    }
	}
      newquads = help;
      //std::cerr << "\t Regular Iteration\t" << count++ << " " << rep << std::endl;
    }
}

/*---------------------------------------------------*/

void LevelMesh2d::InitBoundaryHandler(BoundaryIndexHandler& BI) const
{
  IntSet blines;
  for(int i=0;i<HMP->nblines();i++)
    {
      int q = HMP->bline(i).of_quad();
      if(Quadg2l().find(q)!=Quadg2l().end())
	{
	  blines.insert(i);
	}
    }

  // which colors are there ?
  BI.clear();
  for(IntSet::const_iterator p=blines.begin();
      p!=blines.end(); p++)
    {
      const BoundaryLine& bl = HMP->bline(*p);
      int col = bl.material();

      BI.GetColors().insert(col);
    }
  IntVector colorvec;
  Set2Vec(colorvec,BI.GetColors());

  // compute inverse positions
  std::map<int,int> inv;

  for (int i=0; i<colorvec.size(); i++)
    {
      inv.insert(std::make_pair(colorvec[i],i));
    }
  
  int nc = colorvec.size(); 
  std::vector<std::set<int>  > H1(nc);  // for verteces
  // for cells and local indices
  std::vector<std::set<fixarray<2,int> > >  H2(nc); 

  for(IntSet::const_iterator q=blines.begin();
      q!=blines.end(); q++)
    {
      const BoundaryLine& bl = HMP->bline(*q);
      int col = bl.material();

      std::map<int,int>::const_iterator p = inv.find(col);
      if (p==inv.end())
	{
	  std::cout << "LevelMesh2d::BuildVertexOnBoundary()"<< std::endl;
	  abort();
	}
      int pos = p->second;

      for(int ii=0;ii<2;ii++)
	{
	  int vindex = Vertexg2l(bl.vertex(ii));
	  H1[pos].insert(vindex);
	}
      fixarray<2,int> ind;
      ind[0] = Quadg2l(bl.of_quad());
      ind[1] = bl.edge_in_quad();
      H2[pos].insert(ind);
    }
  BI.CopySetToVector(H1,colorvec,BI.GetVertex());

  for (int i=0; i<H2.size(); i++)
    {
      IntVector v1(H2[i].size());
      IntVector v2(H2[i].size());
      int j = 0;
      
      std::set<fixarray<2,int> >::const_iterator p;
      for (p=H2[i].begin(); p!=H2[i].end(); p++)
	{
	  v1[j] = (*p)[0];
	  v2[j] = (*p)[1];
	  j++;
	}
      int color = colorvec[i];

      BI.GetCell().insert(std::make_pair(color,v1));
      BI.GetLocal().insert(std::make_pair(color,v2));
    }
}

/*--------------------------------------------------------------*/

