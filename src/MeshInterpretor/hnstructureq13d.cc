#include  "hnstructureq13d.h"
#include  "gascoignemesh.h"

/*-----------------------------------------*/

HNStructureQ13d::HNStructureQ13d() : HNStructureQ12d(), faces(NULL)
{}

/*-----------------------------------------*/

HNStructureQ13d::~HNStructureQ13d()
{}

/*--------------------------------------------------------*/

void HNStructureQ13d::ReInit(const MeshInterface* M)
{
  const GascoigneMesh* GM = dynamic_cast<const GascoigneMesh*>(M);
  edges = GM->GetHangingIndexHandler().GetStructure();
  faces = GM->GetHangingIndexHandler().GetStructureFace();
  assert(edges);
  assert(faces);
}

/*-----------------------------------------*/

bool HNStructureQ13d::ZeroCheck(const GlobalVector& u) const
{  
  bool r = HNStructureQ12d::ZeroCheck(u);
  if (r) return r;

  for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
    {
      int i = p->first;
      for(int c=0;c<u.ncomp();c++) 
	{
	  if(u(i,c)!=0.)  return 1;
	}
    }
  return 0;
}

/*-----------------------------------------*/

void HNStructureQ13d::Zero(GlobalVector& u) const
{
  HNStructureQ12d::Zero(u);

  for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
    {
      int i = p->first;
      for(int c=0;c<u.ncomp();c++) u.zero_node(i);
    }
}

/*-----------------------------------------*/

void HNStructureQ13d::Average(GlobalVector& u) const
{
  HNStructureQ12d::Average(u);

  for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
    {
      int i = p->first;
      const fixarray<9,int>& f = p->second;
      u.equ_node(i, 0.25,f[0], 0.25,f[1], 0.25,f[3], 0.25,f[4]);
    }
}

/*-----------------------------------------*/

void HNStructureQ13d::Distribute(GlobalVector& u) const
{
  HNStructureQ12d::Distribute(u);

  for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
    {
      int i = p->first;
      const fixarray<9,int>& f = p->second;

      u.add_node(f[0],0.25,i);
      u.add_node(f[1],0.25,i);
      u.add_node(f[3],0.25,i);
      u.add_node(f[4],0.25,i);
      u.zero_node(i);
    }
}

/*-----------------------------------------*/

int HNStructureQ13d::hanging(int i) const 
{ 
  int r = HNStructureQ12d::hanging(i);

  if (r>0) return 2;
  if (faces->find(i)!=faces->end()) return 4;
  return 0;
}

/*-----------------------------------------*/

void HNStructureQ13d::MatrixDiag(int ncomp, MatrixInterface& A) const
{
  HNStructureQ12d::MatrixDiag(ncomp,A);

  nmatrix<double> M(ncomp);
  M.identity();
  for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
    {
      A.entry_diag(p->first,M);
    }
}

/*-----------------------------------------*/

void HNStructureQ13d::SparseStructureDiag(SparseStructure* S) const
{
  HNStructureQ12d::SparseStructureDiag(S);

  for(const_fiterator p=faces->begin(); p!=faces->end(); p++)
    {
      int i = p->first;
      S->build_add(i,i);
    }
}

/*----------------------------------------------*/

fixarray<4,int> HNStructureQ13d::GetHangingFace(int i) const
{
  const_fiterator p = faces->find(i);
  assert(p!=faces->end());
	  
  fixarray<4,int> Face;
  Face[0] = p->second[0];
  Face[1] = p->second[1];
  Face[2] = p->second[3];
  Face[3] = p->second[4];
  
  return Face;
}

/*----------------------------------------------*/

fixarray<2,int> HNStructureQ13d::GetHangingEdge(int i) const
{
  std::map<int,EdgeVector>::const_iterator p = edges->find(i);
  assert(p!=edges->end());
	  
  fixarray<2,int> Edge;
  Edge[0] = p->second[0];
  Edge[1] = p->second[1];
  
  return Edge;
}

/*----------------------------------------------*/

void HNStructureQ13d::CondenseHanging(nvector<int>& indices) const
{
  CondenseHanging2er(indices);
  CondenseHanging4er(indices);
}

/*----------------------------------------------*/

void HNStructureQ13d::CondenseHanging(EntryMatrix& E, nvector<int>& indices) const
{
  CondenseHanging2er(E,indices);
  CondenseHanging4er(E,indices);
}

/*----------------------------------------------*/

void HNStructureQ13d::CondenseHanging4er(EntryMatrix& E, nvector<int>& indices) const
{
  nvector<int> x(0), y(0);

  for(int ii=0; ii<8; ii++)
    {
      int j = indices[ii];

      if (hanging(j)==4) // 4er haengender Knoten
	{
	  fixarray<4,int> Face = GetHangingFace(j);

	  x.push_back(ii);
	  for(int i=0; i<4; i++)
	    {
	      int FaceIndex = Face[i];
	      //
	      // suche ob FaceIndex schon in indices sind
	      //
	      int jj = 0;
	      bool found = 0;
	      while ( (jj<8) && !found)
		{
		  found = (indices[jj] == FaceIndex);
		  jj++;
		}
	      jj--;
	      if (found) y.push_back(jj); // merke Kopplung in Nachbar vertex
	      else       indices[ii] = FaceIndex; // ersetze Kopplung
	    }
	}
    }
  int n1 = y.size();
  int n2 = x.size();

  assert(n1==3*n2);

  int counter = 0;
  for(int i=0;i<x.size();i++)
    {
      int i1 = x[i];           // new node !
      
      E.multiply_column_row(i1,0.25);
      int last = counter+3;
      
      assert (last<=y.size());

      for ( ; counter<last; counter++)
	{
	  int i2 = y[counter];   // already there
	  E.add_column_row(i2,i1);	  
	}
    }
}

/*----------------------------------------------*/

void HNStructureQ13d::Couplings(nvector<int>& indices) const
{
  // fuer Structure (nicht bigstencil)
  //
  // erst alle haengenden lines ersetzen
  //
  int linecount = 0;
  int quadcount = 0;

  nvector<int>::const_iterator p0 = indices.begin();
  nvector<int>::const_iterator p1 = indices.end();

  for(int i=0; i<8; i++)
    {
      int& ind = indices[i];

      if (hanging(ind)!=2) continue;
      
      linecount++;

      //const IntVector2& line = hang(ind);
      fixarray<2,int> line = GetHangingEdge(ind);
      
      for(int k=0; k<2; k++) 
	{
	  // entweder gibt es newindex schon oder muss hinzugefuegt werden
	  //
	  int newindex = line[k];
	  nvector<int>::const_iterator p = find(p0,p1,newindex);
	  if (p==p1)
	    {
	      ind = newindex;
	      break;
	    }
	}
    }
  // jetzt duerften zu haengenden quads
  // nur noch ein vertex nicht eingetragen sein
  for(int i=0; i<8; i++)
    {
      int& ind = indices[i];

      if (hanging(ind)!=4) continue;

      quadcount++;

      fixarray<4,int> face = GetHangingFace(ind);

      for(int k=0; k<4; k++) 
	{
	  int newindex = face[k];
	  nvector<int>::const_iterator p = find(p0,p1,newindex);
	  if (p==p1)
	    {
	      ind = newindex;
	      break;
	    }
	} 
    }
  assert (quadcount<=3);
  if (quadcount==1) assert (linecount>=2);
  if (quadcount>=2) assert (linecount==3);
}

/*----------------------------------------------*/

void HNStructureQ13d::CondenseHanging2er(EntryMatrix& E, nvector<int>& indices) const
{
  nvector<int> x(0), y(0);

  for(int ii=0; ii<8; ii++)
    {
      int i = indices[ii];

      if (hanging(i)==2) // 2er haengender Knoten
	{
	  fixarray<2,int> Edge = GetHangingEdge(i);

	  x.push_back(ii);

	  for(int iii=0; iii<2; iii++)
	    {
	      int ir = Edge[iii];
	      bool found = 0;
	      for(int iiii=0; (iiii<8) && !found; iiii++)
		{
		  if(ir==indices[iiii])
		    {
		      found = 1;
		      y.push_back(iiii);
		    }
		}
	      if(!found) indices[ii] = ir;
	    }
	}
    }
  assert(x.size()==y.size());
  assert(x.size()<=3);

  for(int i=0;i<x.size();i++)
    {
      int i1 = x[i];        // new node !
      int i2 = y[i];        // already there
      
      E.multiply_column_row(i1,0.5);
      E.add_column_row (i2,i1);
    }
}

/*----------------------------------------------*/

void HNStructureQ13d::CondenseHanging2er(nvector<int>& indices) const
{
  for(int ii=0; ii<8; ii++)
    {
      int i = indices[ii];

      if (hanging(i)==2) // 2er haengender Knoten
	{
	  fixarray<2,int> Edge = GetHangingEdge(i);

	  for(int iii=0; iii<2; iii++)
	    {
	      int ir = Edge[iii];
	      bool found = 0;
	      for(int iiii=0; (iiii<8) && !found; iiii++)
		{
		  if(ir==indices[iiii])
		    {
		      found = 1;
		    }
		}
	      if(!found) indices[ii] = ir;
	    }
	}
    }
}

/*----------------------------------------------*/

void HNStructureQ13d::CondenseHanging4er(nvector<int>& indices) const
{
  for(int ii=0; ii<8; ii++)
    {
      int j = indices[ii];

      if (hanging(j)==4) // 4er haengender Knoten
	{
	  fixarray<4,int> Face = GetHangingFace(j);

	  for(int i=0; i<4; i++)
	    {
	      int FaceIndex = Face[i];
	      //
	      // suche ob FaceIndex schon in indices sind
	      //
	      int jj = 0;
	      bool found = 0;
	      while ( (jj<8) && !found)
		{
		  found = (indices[jj] == FaceIndex);
		  jj++;
		}
	      jj--;
	      if (!found) indices[ii] = FaceIndex; // ersetze Kopplung
	    }
	}
    }
}

/*----------------------------------------------*/


