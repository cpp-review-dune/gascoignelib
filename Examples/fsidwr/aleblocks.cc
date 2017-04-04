
#include "aleblocks.h"

#include "sparseblockmatrix.xx"
#include "sparseblockilu.xx"


namespace Gascoigne
{
  
  
  template class FluidBlock<3>;

  template class SparseBlockMatrix<FluidBlock<3> >;
  template class SparseBlockIlu<FluidBlock<3> >;
  



  

  template<class FD, class FS>
  void CopySubMatrixSolid(SparseBlockMatrix<FD>* D,
			  const SparseBlockMatrix<FS>* S,
			  const std::vector<int>& nodes_l2g,
			  const HASHMAP<int,int>& nodes_g2l,
			  const HASHSET<int>& interface_nodes)
  {
    if (nodes_l2g.size()==0) return;
    
    assert(nodes_l2g.size()>0);
    assert(nodes_l2g.size()==nodes_g2l.size());
    int n = nodes_l2g.size();

    // Copy Stencil
    const ColumnDiagStencil& ST_S = dynamic_cast<const ColumnDiagStencil&> (*(S->GetStencil())); 

    int neglected = 0;
    int included  = 0;
    SparseStructure SS;
    SS.build_begin(n);
    for (int i=0;i<n;++i)
      {
	int rS = nodes_l2g[i];
	for (int pS = ST_S.start(rS);pS<ST_S.stop(rS);++pS)
	  {
	    int cS = ST_S.col(pS);
	    HASHMAP<int,int>::const_iterator hit = nodes_g2l.find(cS);
	    if (hit==nodes_g2l.end()) 
	      {
		++neglected;
		continue;
	      }
	    ++included;
	    SS.build_add(i,hit->second);
	  }
      }
    
    SS.build_end();
    D->ReInit(&SS);
    
    // cout << "New Structure: " << n << " rows " 
    // 	 << included << " couplings, "
    // 	 << neglected << " neglected." << endl;
    

    // Copy Entries
    D->zero();
    const ColumnDiagStencil* ST = dynamic_cast<const ColumnDiagStencil*> (D->GetStencil());
    assert(n==ST->n());
    typename vector<FD>::iterator IT = D->mat().begin();

    for(int r=0;r<n;++r)
      {
	assert(r<nodes_l2g.size());
	int rS = nodes_l2g[r];
	int pS = ST_S.start(rS);

	for (int p=ST->start(r);p!=ST->stop(r);++p,++IT)
	  {
	    int c = ST->col(p);
	    int cS = -1;
	    for (;pS<ST_S.stop(rS); ++pS)
	      {
		int cc = ST_S.col(pS);
		HASHMAP<int,int>::const_iterator it_g2l = nodes_g2l.find(cc);
		if (it_g2l!=nodes_g2l.end())
		  cS = it_g2l->second;
		if (c==cS) break;
	      }
	    assert(c==cS);

	    *IT = *S->mat(pS);

	    
	    // SET DIRICHLET-VALUES in PRESSURE_BLOCK
	    int col = ST_S.col(pS);

	    // Spalte nicht auf 0!
	    // if (interface_nodes.find(col)!=interface_nodes.end())
	    //   {
	    // 	for (int j=0;j<IT->ncomp();++j)
	    // 	  (*IT)(0,j)=0;
	    //   }
	    if (interface_nodes.find(rS)!=interface_nodes.end())
	      {
		for (int j=0;j<IT->ncomp();++j)
		  (*IT)(j,0) =0;

		
		// set 1 to diagonal
		if (col==rS)
		  {
		    (*IT)(0,0)=1.0;
		  }
	      }
	  }
      }
  }


  template<class FD, class FS>  
  void CopySubMatrixFluid(SparseBlockMatrix<FD >* D,
  			  const SparseBlockMatrix<FS >* S,
  			  const std::vector<int>& nodes_l2g,
  			  const HASHMAP<int,int>& nodes_g2l,
  			  const HASHSET<int>& interface_nodes)
  {
    cerr << " template function not written!" << endl;
    abort();
  }
  

  template<>
  void CopySubMatrixFluid(SparseBlockMatrix<FMatrixBlock<7> >* D,
  			  const SparseBlockMatrix<FMatrixBlock<7> >* S,
  			  const std::vector<int>& nodes_l2g,
  			  const HASHMAP<int,int>& nodes_g2l,
  			  const HASHSET<int>& interface_nodes)
  {
    assert(nodes_l2g.size()>0);
    assert(nodes_l2g.size()==nodes_g2l.size());
    int n = nodes_l2g.size();

    // Copy Stencil
    const ColumnDiagStencil& ST_S = dynamic_cast<const ColumnDiagStencil&> (*(S->GetStencil())); 

    int neglected = 0;
    int included  = 0;
    SparseStructure SS;
    SS.build_begin(n);
    for (int i=0;i<n;++i)
      {
  	int rS = nodes_l2g[i];
  	for (int pS = ST_S.start(rS);pS<ST_S.stop(rS);++pS)
  	  {
  	    int cS = ST_S.col(pS);
  	    HASHMAP<int,int>::const_iterator hit = nodes_g2l.find(cS);
  	    if (hit==nodes_g2l.end()) 
  	      {
  		++neglected;
  		continue;
  	      }
  	    ++included;
  	    SS.build_add(i,hit->second);
  	  }
      }
    
    SS.build_end();
    D->ReInit(&SS);
    
    // Copy Entries
    D->zero();
    const ColumnDiagStencil* ST = dynamic_cast<const ColumnDiagStencil*> (D->GetStencil());
    assert(n==ST->n());
    vector<FMatrixBlock<7> >::iterator IT = D->mat().begin();

    for(int r=0;r<n;++r)
      {
  	assert(r<nodes_l2g.size());
  	int rS = nodes_l2g[r];
  	int pS = ST_S.start(rS);

  	for (int p=ST->start(r);p!=ST->stop(r);++p,++IT)
  	  {
  	    int c = ST->col(p);
  	    int cS = -1;
  	    for (;pS<ST_S.stop(rS); ++pS)
  	      {
  		int cc = ST_S.col(pS);
  		HASHMAP<int,int>::const_iterator it_g2l = nodes_g2l.find(cc);
  		if (it_g2l!=nodes_g2l.end())
  		  cS = it_g2l->second;
  		if (c==cS) break;
  	      }
  	    assert(c==cS);

  	    *IT = *S->mat(pS);
	    

  	    // SET BOUNDARY VALUES!!!
  	    // on the interface set boundary values v=0 and u=0
  	    // this is also set in the smoothing-function for the residual

  	    int col = ST_S.col(pS);

	    if (0)
  	    if (interface_nodes.find(col)!=interface_nodes.end())
  	      {
  		for (int i=0;i<3;++i)
  		  for (int j=0;j<7;++j)
  		    {
  		      (*IT)(j,i+1)  =0;
  		      (*IT)(j,i+1+3)=0;
  		    }
  	      }
  	    if (interface_nodes.find(rS)!=interface_nodes.end())
  	      {
		
  		assert(IT->ncomp()==7);
  		for (int i=0;i<3;++i)
  		  for (int j=0;j<7;++j)
  		    {
  		      (*IT)(i+1,j)  =0;
  		      (*IT)(i+1+3,j)=0;
  		    }
  		if (col==rS)
  		  {
  		    for (int i=0;i<3;++i)
  		      {
  			(*IT)(i+1,i+1) = 1.0;
  			(*IT)(i+4,i+4) = 1.0;
  		      }
  		  }
  	      }

  	    for (int i=0;i<3;++i)
  	      for (int j=0;j<4;++j)
  		assert((*IT)(i+4,j)==0);
  	    for (int i=0;i<3;++i)
  	      for (int j=0;j<3;++j)
  		if (i!=j)
  		  assert((*IT)(i+4,j+4)==0);
  	  }
      }
  }






  template<>
  void CopySubMatrixFluid(SparseBlockMatrix<FMatrixBlock<5> >* D,
			  const SparseBlockMatrix<FMatrixBlock<5> >* S,
			  const std::vector<int>& nodes_l2g,
			  const HASHMAP<int,int>& nodes_g2l,
			  const HASHSET<int>& interface_nodes)
  {
    assert(nodes_l2g.size()>0);
    assert(nodes_l2g.size()==nodes_g2l.size());
    int n = nodes_l2g.size();

    // Copy Stencil
    const ColumnDiagStencil& ST_S = dynamic_cast<const ColumnDiagStencil&> (*(S->GetStencil())); 

    int neglected = 0;
    int included  = 0;
    SparseStructure SS;
    SS.build_begin(n);
    for (int i=0;i<n;++i)
      {
	int rS = nodes_l2g[i];
	for (int pS = ST_S.start(rS);pS<ST_S.stop(rS);++pS)
	  {
	    int cS = ST_S.col(pS);
	    HASHMAP<int,int>::const_iterator hit = nodes_g2l.find(cS);
	    if (hit==nodes_g2l.end()) 
	      {
		++neglected;
		continue;
	      }
	    ++included;
	    SS.build_add(i,hit->second);
	  }
      }
    
    SS.build_end();
    D->ReInit(&SS);
    
    // Copy Entries
    D->zero();
    const ColumnDiagStencil* ST = dynamic_cast<const ColumnDiagStencil*> (D->GetStencil());
    assert(n==ST->n());
    vector<FMatrixBlock<5> >::iterator IT = D->mat().begin();

    for(int r=0;r<n;++r)
      {
	assert(r<nodes_l2g.size());
	int rS = nodes_l2g[r];
	int pS = ST_S.start(rS);

	for (int p=ST->start(r);p!=ST->stop(r);++p,++IT)
	  {
	    int c = ST->col(p);
	    int cS = -1;
	    for (;pS<ST_S.stop(rS); ++pS)
	      {
		int cc = ST_S.col(pS);
		HASHMAP<int,int>::const_iterator it_g2l = nodes_g2l.find(cc);
		if (it_g2l!=nodes_g2l.end())
		  cS = it_g2l->second;
		if (c==cS) break;
	      }
	    assert(c==cS);

	    *IT = *S->mat(pS);
	    

	    // SET BOUNDARY VALUES!!!
	    // on the interface set boundary values v=0 and u=0
	    // this is also set in the smoothing-function for the residual

	    int col = ST_S.col(pS);

	    if (interface_nodes.find(col)!=interface_nodes.end())
	      {
		for (int i=0;i<2;++i)
		  for (int j=0;j<5;++j)
		    {
		      (*IT)(j,i+1)  =0;
		      (*IT)(j,i+1+2)=0;
		    }
	      }
	    if (interface_nodes.find(rS)!=interface_nodes.end())
	      {
		assert(IT->ncomp()==5);
		for (int i=0;i<2;++i)
		  for (int j=0;j<5;++j)
		    {
		      (*IT)(i+1,j)  =0;
		      (*IT)(i+1+2,j)=0;
		    }

		if (col==rS)
		  {
		    for (int i=0;i<2;++i)
		      {
			(*IT)(i+1,i+1) = 1.0;
			(*IT)(i+1+2,i+1+2) = 1.0;
		      }
		  }
	      }
	  }
      }
  }

  
  template<>
  void CopySubMatrixFluid(SparseBlockMatrix<FluidBlock<3> >* D,
			  const SparseBlockMatrix<FMatrixBlock<7> >* S,
			  const std::vector<int>& nodes_l2g,
			  const HASHMAP<int,int>& nodes_g2l,
			  const HASHSET<int>& interface_nodes)
  {
    assert(nodes_l2g.size()>0);
    assert(nodes_l2g.size()==nodes_g2l.size());
    int n = nodes_l2g.size();

    // Copy Stencil
    const ColumnDiagStencil& ST_S = dynamic_cast<const ColumnDiagStencil&> (*(S->GetStencil())); 

    int neglected = 0;
    int included  = 0;
    SparseStructure SS;
    SS.build_begin(n);
    for (int i=0;i<n;++i)
      {
	int rS = nodes_l2g[i];
	for (int pS = ST_S.start(rS);pS<ST_S.stop(rS);++pS)
	  {
	    int cS = ST_S.col(pS);
	    HASHMAP<int,int>::const_iterator hit = nodes_g2l.find(cS);
	    if (hit==nodes_g2l.end()) 
	      {
		++neglected;
		continue;
	      }
	    ++included;
	    SS.build_add(i,hit->second);
	  }
      }
    
    SS.build_end();
    D->ReInit(&SS);
    

    // Copy Entries
    D->zero();
    const ColumnDiagStencil* ST = dynamic_cast<const ColumnDiagStencil*> (D->GetStencil());
    assert(n==ST->n());
    vector<FluidBlock<3> >::iterator IT = D->mat().begin();

    for(int r=0;r<n;++r)
      {
	assert(r<nodes_l2g.size());
	int rS = nodes_l2g[r];
	int pS = ST_S.start(rS);

	for (int p=ST->start(r);p!=ST->stop(r);++p,++IT)
	  {
	    int c = ST->col(p);
	    int cS = -1;
	    for (;pS<ST_S.stop(rS); ++pS)
	      {
		int cc = ST_S.col(pS);
		HASHMAP<int,int>::const_iterator it_g2l = nodes_g2l.find(cc);
		if (it_g2l!=nodes_g2l.end())
		  cS = it_g2l->second;
		if (c==cS) break;
	      }
	    assert(c==cS);

	    *IT = *S->mat(pS);
	    // XXX
	    int col = ST_S.col(pS);
	    assert(IT->ncomp()==7);
	    if (interface_nodes.find(col)!=interface_nodes.end())
	      {
		// delete column entries
		IT->__S  = 0;
		for (int i=0;i<4;++i)
		  for (int j=0;j<6;++j)
		    (*IT)(i,j+1) = 0;
	      }
	    if (interface_nodes.find(rS)!=interface_nodes.end())
	      {
		// delete row
		IT->__S = 0;
		
		for (int i=0;i<3;++i)
		  for (int j=0;j<7;++j)
		    (*IT)(i+1,j)  =0;


		// set 1 to diagonal
		if (col==rS)
		  {
		    for (int i=0;i<3;++i)
		      (*IT)(i+1,i+1) = 1.0;
		    IT->__S = 1.0;
		  }
	      }

	  }
	
      }
  }
  
    
  template 
  void CopySubMatrixSolid(SparseBlockMatrix<FMatrixBlock<7> >* D,
			  const SparseBlockMatrix<FMatrixBlock<7> >* S,
			  const std::vector<int>& nodes_l2g,
			  const HASHMAP<int,int>& nodes_g2l,		     
			  const HASHSET<int>& interface_nodes);
  
  template
  void CopySubMatrixFluid(SparseBlockMatrix<FMatrixBlock<5> >* D,
			  const SparseBlockMatrix<FMatrixBlock<5> >* S,
			  const std::vector<int>& nodes_l2g,
			  const HASHMAP<int,int>& nodes_g2l,
			  const HASHSET<int>& interface_nodes);
  // template
  // void CopySubMatrixFluid(SparseBlockMatrix<FMatrixBlock<7> >* D,
  // 			  const SparseBlockMatrix<FMatrixBlock<7> >* S,
  // 			  const std::vector<int>& nodes_l2g,
  // 			  const HASHMAP<int,int>& nodes_g2l,
  // 			  const HASHSET<int>& interface_nodes);


  template
  void CopySubMatrixSolid(SparseBlockMatrix<FMatrixBlock<5> >* D,
			  const SparseBlockMatrix<FMatrixBlock<5> >* S,
			  const std::vector<int>& nodes_l2g,
			  const HASHMAP<int,int>& nodes_g2l,		     
			  const HASHSET<int>& interface_nodes);
  
    
  

}
