#include  "aleinterpolator.h"
#include  "gascoignemeshtransfer.h"


using namespace std;

#define __DIM__ 2

/*--------------------------------------------------------*/

namespace Gascoigne
{

  /*-----------------------------------------*/
  
  void AleInterpolator::restrict_zero(GlobalVector& uL, const GlobalVector& ul) const
  {

    // we need to check, that no values are copied onto the interface 
    // from the side, where test-functions are deleted.

    assert((ul.ncomp()==5)||(ul.ncomp()==7));
    // int DIM = 0;
    // if      (ul.ncomp() == 5) DIM = 2;
    // else if (ul.ncomp() == 7) DIM = 3;
    // else abort();
    

    for(int i=0;i<c2f.size();i++)  uL.equ_node(i,1.,c2f[i],ul);
    for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
	p!=zweier.end();p++) 
      {
	int il = p->first;
	fixarray<2,int> n2 = p->second;
	for (int i=0;i<2;++i)
	  {
	    if (__inodes->find(n2[i])==__inodes->end())
	      {
		uL(n2[i],0) += 0.5 * ul(il,0);
		for (int c=1;c<uL.ncomp();++c)
		  uL(n2[i],c) += 0.5 * ul(il,c);
	      }
	    else
	      {
		if (__snodes->find(il)!=__snodes->end())
		  uL(n2[i],0) += 0.5 * ul(il,0);

		if (__fnodes->find(il)!=__fnodes->end())
		  for (int c=1;c<uL.ncomp();++c)
		    uL(n2[i],c) += 0.5 * ul(il,c);
	      }
	    
	  }
	
      }
    for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
	p!=vierer.end();p++) 
      {
	int il = p->first;
	fixarray<4,int> n4 = p->second;
	for (int i=0;i<4;++i)
	  {
	    if (__inodes->find(n4[i])==__inodes->end())
	      {
		uL(n4[i],0) += 0.25 * ul(il,0);
		for (int c=1;c<uL.ncomp();++c)
		  uL(n4[i],c) += 0.25 * ul(il,c);
	      }
	    else
	      {
		if (__snodes->find(il)!=__snodes->end())
		  uL(n4[i],0) += 0.25 * ul(il,0);
		
		if (__fnodes->find(il)!=__fnodes->end())
		  for (int c=1;c<uL.ncomp();++c)
		    uL(n4[i],c) += 0.25 * ul(il,c);
	      }
	  }
      }
    for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
	p!=achter.end();p++) 
      {
	int il = p->first;
	fixarray<8,int> n8 = p->second;
	for (int i=0; i<8; i++)
	  {
	    if (__inodes->find(n8[i])==__inodes->end())
	      {
		uL(n8[i],0) += 0.125 * ul(il,0);
		for (int c=1;c<uL.ncomp();++c)
		  uL(n8[i],c) += 0.125 * ul(il,c);
	      }
	    else
	      {
		if (__snodes->find(il)!=__snodes->end())
		  uL(n8[i],0) += 0.125 * ul(il,0);
		
		if (__fnodes->find(il)!=__fnodes->end())
		  for (int c=1;c<uL.ncomp();++c)
		    uL(n8[i],c) += 0.125 * ul(il,c);
	      }	    
	  }
      }
  }

  /*-----------------------------------------*/

  void AleInterpolator::prolongate_add(GlobalVector& ul, const GlobalVector& uL) const
  {
    for(int i=0;i<c2f.size();i++)  ul.add_node(c2f[i],1.,i,uL);
    for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
	p!=zweier.end();p++) 
      {
	int il = p->first;
	fixarray<2,int> n2 = p->second;

	if (__inodes->find(il)==__inodes->end())
	  for (int i=0;i<2;++i)
	    {
	      ul(il,0) += 0.5 * uL(n2[i],0);
	      for (int c=1;c<ul.ncomp();++c)
		ul(il,c) += 0.5 * uL(n2[i],c);
	    }
	else
	  for (int i=0;i<2;++i)
	    {
	      if (__snodes->find(n2[i])!=__snodes->end())
		ul(il,0) += 0.5 * uL(n2[i],0);
	      if (__fnodes->find(n2[i])!=__fnodes->end())
		for (int c=1;c<ul.ncomp();++c)
		  ul(il,c) += 0.5 * uL(n2[i],c);
	    }
      }

    for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
	p!=vierer.end();p++) 
      {
	int il = p->first;
	fixarray<4,int> n4 = p->second;
	if (__inodes->find(il)==__inodes->end())
	  for (int i=0;i<4;++i)
	    {
	      ul(il,0) += 0.25 * uL(n4[i],0);
	      for (int c=1;c<ul.ncomp();++c)
		ul(il,c) += 0.25 * uL(n4[i],c);
	    }
	else
	  for (int i=0;i<4;++i)
	    {
	      if (__snodes->find(n4[i])!=__snodes->end())
		ul(il,0) += 0.25 * uL(n4[i],0);
	      if (__fnodes->find(n4[i])!=__fnodes->end())
		for (int c=1;c<ul.ncomp();++c)
		  ul(il,c) += 0.25 * uL(n4[i],c);
	    }
      }
    for(map<int,fixarray<8,int> >::const_iterator p=achter.begin();
	p!=achter.end();p++) 
      {
	int il = p->first;
	fixarray<8,int> n8 = p->second;
	if (__inodes->find(il)==__inodes->end())
	  for (int i=0; i<8; i++)
	    {
	      ul(il,0) += 0.125 * uL(n8[i],0);
	      for (int c=1;c<ul.ncomp();++c)
		ul(il,c) += 0.125 * uL(n8[i],c);
	    }
	else
	  for (int i=0; i<8; i++)
	    {
	      if (__snodes->find(n8[i])!=__snodes->end())
		ul(il,0) += 0.125 * uL(n8[i],0);
	      if (__fnodes->find(n8[i])!=__fnodes->end())
		for (int c=1;c<ul.ncomp();++c)
		  ul(il,c) += 0.125 * uL(n8[i],c);
	    }
      }
  }

  // /*-----------------------------------------*/

  // void AleInterpolator::SolutionTransfer(GlobalVector& uL, const GlobalVector& ul) const
  // {
  //   for(int i=0;i<c2f.size();i++) uL.equ_node(i,1.,c2f[i],ul);
  // }

  // /*-----------------------------------------*/

  // void AleInterpolator::Pi(GlobalVector& u) const
  // {
  //   for(map<int,fixarray<2,int> >::const_iterator p=zweier.begin();
  // 	p!=zweier.end();p++) {
  //     int il = p->first;
  //     fixarray<2,int> n2 = p->second;
  //     u.add_node(il,-0.5,c2f[n2[0]],u);
  //     u.add_node(il,-0.5,c2f[n2[1]],u);
  //   }
  //   for(map<int,fixarray<4,int> >::const_iterator p=vierer.begin();
  // 	p!=vierer.end();p++) {
  //     int il = p->first;
  //     fixarray<4,int> n4 = p->second;
  //     u.add_node(il,-0.25,c2f[n4[0]],u);
  //     u.add_node(il,-0.25,c2f[n4[1]],u);
  //     u.add_node(il,-0.25,c2f[n4[2]],u);
  //     u.add_node(il,-0.25,c2f[n4[3]],u);
  //   }
  //   for(int i=0;i<c2f.size();i++)  u.node_zero(c2f[i]);
  // }
}
