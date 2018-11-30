#include "fsimatrix.h"

namespace Gascoigne
{
  
  /*
   * Construct two small matrices: 
   *    __F only fluid-nodes, p, v, u, 
   *    __S only solid-nodes, v,u
   */
  template<int DIM>
  void FSIMatrix<DIM>::ReInit(const SparseStructureInterface* _S)
  {
    const SparseStructure* S = dynamic_cast<const SparseStructure*> (_S);
    assert(S);
    
    SparseStructure SF,SS;
    SF.build_begin(n_fluid_nodes());
    SS.build_begin(n_solid_nodes());

    for (int r=0;r<S->n();++r)
      {
	int fr = Fg2l(r);
	int sr = Sg2l(r);
	
	const IntSet& row = S->row(r);
	for (IntSet::const_iterator it = row.begin();it!=row.end();++it)
	  {
	    int c = *it;
	    int fc = Fg2l(c);
	    int sc = Sg2l(c);
	    if ((fr!=-1)&&(fc!=-1)) SF.build_add(fr,fc);
	    if ((sr!=-1)&&(sc!=-1)) SS.build_add(sr,sc);

	    if (  ((fr==-1)||(fc==-1)) && ((sr==-1)||(sc==-1)) )
	      {
		cerr << "row/col " << r << " " << c << " not in solid nor fluid!" << endl;
		abort();
	      }
	  }
      }
    SF.build_end();
    SS.build_end();
    cout << "Stencil F: " << SF.n() << "\t" << SF.ntotal() << endl;
    cout << "Stencil S: " << SS.n() << "\t" << SS.ntotal() << endl;

    _xs1.ncomp()=DIM; _xs1.resize(SS.n());
    _ys1.ncomp()=DIM; _ys1.resize(SS.n());
    _xs2.ncomp()=DIM; _xs2.resize(SS.n());
    _ys2.ncomp()=DIM; _ys2.resize(SS.n());

    _xf_ns.ncomp()=DIM+1; _xf_ns.resize(SF.n());
    _yf_ns.ncomp()=DIM+1; _yf_ns.resize(SF.n());
    _xf_ext.ncomp()=DIM; _xf_ext.resize(SF.n());
    _yf_ext.ncomp()=DIM; _yf_ext.resize(SF.n());
    
    __S11.ReInit(&SS);
    __S12.ReInit(&SS);
    __S21.ReInit(&SS);
    __S22.ReInit(&SS);
    __F_NS.ReInit(&SF);
    __F_EXT.ReInit(&SF);
    __F_ALE.ReInit(&SF);
  }

  template<int DIM>
  void FSIMatrix<DIM>::entry(nvector<int>::const_iterator start, 
			     nvector<int>::const_iterator stop, 
			     const EntryMatrix& M, double s)
  {
    int nfluid=0;
    int nsolid=0;
    
       
    nvector<int>::const_iterator p, q;
    int nnodes=0;
    
    for (p=start;p!=stop;++p)
      {
	if (fluid_node(*p)) ++nfluid;
	if (solid_node(*p)) ++nsolid;
	++nnodes;
      }
    
    assert((nfluid==nnodes)||(nsolid==nnodes));
    
    if (nfluid==nnodes) 
      {
	nvector<int> fi;
	for (nvector<int>::const_iterator it = start;it!=stop;++it)
	  fi.push_back(Fg2l(*it));
	for (int i=0;i<fi.size();++i) assert(fi[i]>=0);

	__F_NS.entry_some (0,    0,     fi.begin(),fi.end(),M,s);
	__F_ALE.entry_some(0,    DIM+1, fi.begin(),fi.end(),M,s);
	__F_EXT.entry_some(DIM+1,DIM+1, fi.begin(),fi.end(),M,s);

      }
    else if (nsolid==nnodes) 
      {
	nvector<int> si;
	for (nvector<int>::const_iterator it = start;it!=stop;++it)
	  si.push_back(Sg2l(*it));
	for (int i=0;i<si.size();++i) assert(si[i]>=0);
	__S11.entry_some(1,1,         si.begin(),si.end(),M,s);
	__S12.entry_some(1,1+DIM,     si.begin(),si.end(),M,s);
	__S21.entry_some(1+DIM,1,     si.begin(),si.end(),M,s);
	__S22.entry_some(1+DIM,1+DIM, si.begin(),si.end(),M,s);
      }
    else
      {
	cerr << "FSIMatrix entry: " << nfluid << " " << nsolid 
	     << " " << nnodes << endl;
	assert(0);
      }
  }

  template<int DIM>
  void FSIMatrix<DIM>::dirichlet (const nvector<int>& bv, const std::vector<int>& cv) 
  {
    nvector<int> bf,bs;
    for (int i=0;i<bv.size();++i)
      {
	if (fluid_node(bv[i])) bf.push_back(Fg2l(bv[i]));
	if (solid_node(bv[i])) bs.push_back(Sg2l(bv[i]));
      }
    vector<int> cns,cext,cs1,cs2;
    for (int i=0;i<cv.size();++i)
      {
	int c = cv[i];
	if (c<=DIM) cns.push_back(c);
	else cext.push_back(c-DIM-1);
	
	if (c<=DIM) cs1.push_back(c-1);
	else cs2.push_back(c-1-DIM);	
      }
    
    // fluid
    __F_NS.dirichlet(bf,cns);
    __F_EXT.dirichlet(bf,cext);
    __F_ALE.dirichlet_only_row_no_diag(bf,cns);
    __F_ALE.dirichlet_only_column(bf,cext);
    
    // solid
    __S11.dirichlet(bs,cs1);
    __S22.dirichlet(bs,cs2);
    __S12.dirichlet_only_row_no_diag(bs,cs1);
    __S21.dirichlet_only_row_no_diag(bs,cs2);
    
    __S12.dirichlet_only_column(bs,cs2);
    __S21.dirichlet_only_column(bs,cs1);
  }









  template<int DIM>
  void FSIMatrix<DIM>::dirichlet (int ii, const std::vector<int>& cv)
  {
    if (fluid_node(ii)) 
      {
	//	__F.dirichlet(Fg2l(ii),cv);




	vector<int> cv_ns;
	vector<int> cv_ext;
	for (int i=0;i<cv.size();++i)
	  {
	    if (cv[i]<=DIM)  cv_ns.push_back(cv[i]);
	    else cv_ext.push_back(cv[i]-DIM-1);
	  }
	for (int i=0;i<cv_ns.size();++i) assert(cv_ns[i]<=DIM);
	for (int i=0;i<cv_ext.size();++i) assert(cv_ext[i]<DIM);
	for (int i=0;i<cv_ext.size();++i) assert(cv_ext[i]>=0);

	__F_NS .dirichlet(Fg2l(ii),cv_ns);
	__F_EXT.dirichlet(Fg2l(ii),cv_ext);
	
	__F_ALE.dirichlet_only_row_no_diag(Fg2l(ii),cv_ns);
	__F_ALE.dirichlet_only_column(Fg2l(ii)     ,cv_ext);
      }
    
    if (solid_node(ii)) 
      {
	nvector<int> cv1,cv2;
	for (int i=0;i<cv.size();++i) 
	  {
	    if (cv[i]<=DIM) cv1.push_back(cv[i]-1);
	    else cv2.push_back(cv[i]-1-DIM);
	  }
	
	for (int i=0;i<cv1.size();++i) assert((cv1[i]>=0)&&(cv1[i]<DIM));
	for (int i=0;i<cv2.size();++i) assert((cv2[i]>=0)&&(cv2[i]<DIM));
	
	__S11.dirichlet(Sg2l(ii),cv1);
	__S22.dirichlet(Sg2l(ii),cv2);

	__S12.dirichlet_only_row_no_diag(Sg2l(ii),cv1);
	__S21.dirichlet_only_row_no_diag(Sg2l(ii),cv2);

	__S12.dirichlet_only_column(Sg2l(ii),cv2);
	__S21.dirichlet_only_column(Sg2l(ii),cv1);
      }
  }


  template<int DIM>
  void FSIMatrix<DIM>::Fluid_g2l_set(GlobalVector& xf, const GlobalVector& x) const
  {
    assert(xf.n()     == __fluid_l2g->size());
    assert(xf.ncomp() == FCOMP);
    int l=0;
    for (vector<int>::const_iterator it = __fluid_l2g->begin();it!=__fluid_l2g->end();++it,++l)
      xf.equ_node(l,*it,x);
  }
  template<int DIM>
  void FSIMatrix<DIM>::Solid1_g2l_set(GlobalVector& xs, const GlobalVector& x) const
  {
    assert(xs.n()     == __solid_l2g->size());
    assert(xs.ncomp() == DIM);
    int l=0;
    for (vector<int>::const_iterator it = __solid_l2g->begin();it!=__solid_l2g->end();++it,++l)
      for (int c=0;c<DIM;++c)
	xs(l,c)=x(*it,c+1);
  }
  template<int DIM>
  void FSIMatrix<DIM>::Solid2_g2l_set(GlobalVector& xs, const GlobalVector& x) const
  {
    assert(xs.n()     == __solid_l2g->size());
    assert(xs.ncomp() == DIM);
    int l=0;
    for (vector<int>::const_iterator it = __solid_l2g->begin();it!=__solid_l2g->end();++it,++l)
      for (int c=0;c<DIM;++c)
	xs(l,c)=x(*it,c+1+DIM);
  }



  template<int DIM>
  void FSIMatrix<DIM>::Fluid_l2g_add(GlobalVector& x, const GlobalVector& xf, double s) const
  {
    assert(xf.n()     == __fluid_l2g->size());
    assert(xf.ncomp() == FCOMP);
    int l=0;
    for (vector<int>::const_iterator it = __fluid_l2g->begin();it!=__fluid_l2g->end();++it,++l)
      x.add_node(*it,s,l,xf);
  }

  template<int DIM>
  void FSIMatrix<DIM>::Solid1_l2g_add(GlobalVector& x, const GlobalVector& xs,double s) const
  {
    assert(xs.n()     == __solid_l2g->size());
    assert(xs.ncomp() == DIM);
    int l=0;
    for (vector<int>::const_iterator it = __solid_l2g->begin();it!=__solid_l2g->end();++it,++l)
      for (int c=0;c<DIM;++c)
	x(*it,c+1) += s*xs(l,c);
  }
  template<int DIM>
  void FSIMatrix<DIM>::Solid2_l2g_add(GlobalVector& x, const GlobalVector& xs,double s) const
  {
    assert(xs.n()     == __solid_l2g->size());
    assert(xs.ncomp() == DIM);
    int l=0;
    for (vector<int>::const_iterator it = __solid_l2g->begin();it!=__solid_l2g->end();++it,++l)
      for (int c=0;c<DIM;++c)
	x(*it,c+1+DIM) += s*xs(l,c);
  }




  template<int DIM>
  void FSIMatrix<DIM>::Fluid_NS_g2l_set(GlobalVector& xf, const GlobalVector& x) const
  {
    assert(xf.n()     == __fluid_l2g->size());
    assert(xf.ncomp() == DIM+1);
    int l=0;
    for (vector<int>::const_iterator it = __fluid_l2g->begin();it!=__fluid_l2g->end();++it,++l)
      for (int c=0;c<DIM+1;++c)
	xf(l,c) = x(*it,c);
  }
  template<int DIM>
  void FSIMatrix<DIM>::Fluid_NS_l2g_add(GlobalVector& x, const GlobalVector& xf, double s) const
  {
    assert(xf.n()     == __fluid_l2g->size());
    assert(xf.ncomp() == DIM+1);
    int l=0;
    for (vector<int>::const_iterator it = __fluid_l2g->begin();it!=__fluid_l2g->end();++it,++l)
      for (int c=0;c<DIM+1;++c)
	x(*it,c) += s*xf(l,c);
  }

  template<int DIM>
  void FSIMatrix<DIM>::Fluid_EXT_g2l_set(GlobalVector& xf, const GlobalVector& x) const
  {
    assert(xf.n()     == __fluid_l2g->size());
    assert(xf.ncomp() == DIM);
    assert(xf.ncomp()+DIM+1==x.ncomp());
    int l=0;
    for (vector<int>::const_iterator it = __fluid_l2g->begin();it!=__fluid_l2g->end();++it,++l)
      for (int c=0;c<DIM;++c)
	xf(l,c) = x(*it,c+DIM+1);
  }
  template<int DIM>
  void FSIMatrix<DIM>::Fluid_EXT_l2g_add(GlobalVector& x, const GlobalVector& xf, double s) const
  {
    assert(xf.n()     == __fluid_l2g->size());
    assert(xf.ncomp() == DIM);
    assert(xf.ncomp()+DIM+1==x.ncomp());
    int l=0;
    for (vector<int>::const_iterator it = __fluid_l2g->begin();it!=__fluid_l2g->end();++it,++l)
      for (int c=0;c<DIM;++c)
	x(*it,c+DIM+1) += s*xf(l,c);
  }

  
  

  
  template<int DIM>
  void FSIMatrix<DIM>::vmult(GlobalVector& y, const GlobalVector& x, double s) const
  {

    Fluid_NS_g2l_set(_xf_ns,x);
    Fluid_EXT_g2l_set(_xf_ext,x);
    
    // NS
    _yf_ns.zero();
    __F_NS.vmult(_yf_ns,_xf_ns,s);
    __F_ALE.vmult(_yf_ns,_xf_ext,s);
    Fluid_NS_l2g_add(y,_yf_ns,1.0);
    
    // EXT
    //    Fluid_EXT_g2l_set(_xf_ext,x);
    _yf_ext.zero();
    __F_EXT.vmult(_yf_ext,_xf_ext,s);
    Fluid_EXT_l2g_add(y,_yf_ext,1.0);
    
    // Solid
    Solid1_g2l_set(_xs1,x);
    Solid2_g2l_set(_xs2,x);
    _ys1.zero();
    __S11.vmult(_ys1,_xs1,s);
    __S12.vmult(_ys1,_xs2,s);
    _ys2.zero();
    __S21.vmult(_ys2,_xs1,s);
    __S22.vmult(_ys2,_xs2,s);
    Solid1_l2g_add(y,_ys1,1.0);
    Solid2_l2g_add(y,_ys2,1.0);
  }
  
    
  


  template class FSIMatrix<2>;
  template class FSIMatrix<3>;
  

}



