#include  "mginterpolatornested.h"
#include  "gascoignemeshtransfer.h"

/*--------------------------------------------------------*/

void MgInterpolatorNested::init(const MeshTransferInterface* MT)
{
  const GascoigneMeshTransfer* GT = dynamic_cast<const GascoigneMeshTransfer*>(MT);

  c2f.reservesize(GT->GetC2f().size());

  zweier = GT->GetZweier();
  vierer = GT->GetVierer();
  achter = GT->GetAchter();
  c2f    = GT->GetC2f();
}

/*-----------------------------------------*/
  
void MgInterpolatorNested::restrict_zero(Vector& uL, const Vector& ul) const
{
  for(int i=0;i<c2f.size();i++)  uL.equ_node(i,1.,c2f[i],ul);
  for(std::map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) {
    int il = p->first;
    fixarray<2,int> n2 = p->second;
    uL.add_node(n2[0],0.5,il,ul);
    uL.add_node(n2[1],0.5,il,ul);
  }
  for(std::map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    uL.add_node(n4[0],0.25,il,ul);
    uL.add_node(n4[1],0.25,il,ul);
    uL.add_node(n4[2],0.25,il,ul);
    uL.add_node(n4[3],0.25,il,ul);
  }
  for(std::map<int,fixarray<8,int> >::const_iterator p=achter.begin();
      p!=achter.end();p++) {
    int il = p->first;
    fixarray<8,int> n8 = p->second;
    for (int i=0; i<8; i++)
      {
	uL.add_node(n8[i],0.125,il,ul);
      }
  }
}

/*-----------------------------------------*/

void MgInterpolatorNested::prolongate_add(Vector& ul, const Vector& uL) const
{
  for(int i=0;i<c2f.size();i++)  ul.add_node(c2f[i],1.,i,uL);
  for(std::map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) 
    {
      int il = p->first;
      fixarray<2,int> n2 = p->second;
      ul.add_node(il,0.5,n2[0],uL);
      ul.add_node(il,0.5,n2[1],uL);
    }
  for(std::map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) 
    {
      int il = p->first;
      fixarray<4,int> n4 = p->second;
      ul.add_node(il,0.25,n4[0],uL);
      ul.add_node(il,0.25,n4[1],uL);
      ul.add_node(il,0.25,n4[2],uL);
      ul.add_node(il,0.25,n4[3],uL);
    }
  for(std::map<int,fixarray<8,int> >::const_iterator p=achter.begin();
      p!=achter.end();p++) 
    {
      int il = p->first;
      fixarray<8,int> n8 = p->second;
      for (int i=0; i<8; i++)
	{
	  ul.add_node(il,0.125,n8[i],uL);
	}
    }
}

/*-----------------------------------------*/

void MgInterpolatorNested::SolutionTransfer(Vector& uL, const Vector& ul) const
{
  for(int i=0;i<c2f.size();i++)  uL.equ_node(i,1.,c2f[i],ul);
}

/*-----------------------------------------*/

void MgInterpolatorNested::Pi(Vector& u) const
{
  for(std::map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
      p!=zweier.end();p++) {
    int il = p->first;
    fixarray<2,int> n2 = p->second;
    u.add_node(il,-0.5,c2f[n2[0]],u);
    u.add_node(il,-0.5,c2f[n2[1]],u);
  }
  for(std::map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
      p!=vierer.end();p++) {
    int il = p->first;
    fixarray<4,int> n4 = p->second;
    u.add_node(il,-0.25,c2f[n4[0]],u);
    u.add_node(il,-0.25,c2f[n4[1]],u);
    u.add_node(il,-0.25,c2f[n4[2]],u);
    u.add_node(il,-0.25,c2f[n4[3]],u);
  }
  for(int i=0;i<c2f.size();i++)  u.node_zero(c2f[i]);
}
