#ifndef __dynamicblockilu_h
#define __dynamicblockilu_h

#include  "dynamicblockmatrix.h"
#include  "iluinterface.h"

/*-------------------------------------------------------------*/

namespace Gascoigne
{
    template<class B>
	class DynamicBlockIlu: public virtual IluInterface, public DynamicBlockMatrix<B>
    {
    protected:

	typedef typename list<int>::const_iterator    const_citerator;
	typedef typename list<int>::iterator          citerator;
	typedef typename list<B>::const_iterator      const_viterator;
	typedef typename list<B>::iterator            viterator;


	nvector<int>          p,q;
	GlobalVector*         yp;
    
	void backward() const;
	void forward () const;
	virtual void hin(const GlobalVector& x) const;
	virtual void her(GlobalVector& x) const;

	int   n()          const { return DynamicBlockMatrix<B>::n();};

    public:

	DynamicBlockIlu<B>();
	DynamicBlockIlu<B>(const DynamicBlockIlu<B>& I);
	~DynamicBlockIlu();



	viterator       vdiag (int i)       { return DynamicBlockMatrix<B>::vdiag (i);}
	const_viterator vdiag (int i) const { return DynamicBlockMatrix<B>::vdiag (i);}
	citerator       cdiag (int i)       { return DynamicBlockMatrix<B>::cdiag (i);}
	const_citerator cdiag (int i) const { return DynamicBlockMatrix<B>::cdiag (i);}
	const_citerator cstart(int i) const { return DynamicBlockMatrix<B>::cstart(i);}
	const_citerator cstop (int i) const { return DynamicBlockMatrix<B>::cstop (i);}
	citerator       cstart(int i)       { return DynamicBlockMatrix<B>::cstart(i);}
	citerator       cstop (int i)       { return DynamicBlockMatrix<B>::cstop (i);}
	const_viterator vstart(int i) const { return DynamicBlockMatrix<B>::vstart(i);}
	const_viterator vstop (int i) const { return DynamicBlockMatrix<B>::vstop (i);}
	viterator       vstart(int i)       { return DynamicBlockMatrix<B>::vstart(i);}
	viterator       vstop (int i)       { return DynamicBlockMatrix<B>::vstop (i);}



	string GetName() const {return "DynamicBlockIlu";}
  
	nvector<int>&       GetP() {return p;}
	nvector<int>&       GetQ() {return q;}
	const nvector<int>& GetP() const {return p;}
	const nvector<int>& GetQ() const {return q;}

	void modify(int c, double s);
	void zero() { DynamicBlockMatrix<B>::zero(); }

	void compute_ilu ();
	void ReInit      (const SparseStructureInterface* SI);
	void ConstructStructure(const nvector<int>& perm, const MatrixInterface& A);
	void copy_entries(const MatrixInterface* A);
	void solve       (GlobalVector& x) const;
	void solvetrans  (GlobalVector& x) const { assert(0);};
    };
}

#endif
