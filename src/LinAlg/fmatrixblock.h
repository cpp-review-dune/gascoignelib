#ifndef __fmatrixblock_h
#define __fmatrixblock_h

#include  "entrymatrix.h"
#include  "nodematrix.h"

/*****************************************************/

namespace Gascoigne
{
template<int N>
class FMatrixBlock : public NodeMatrix<N,float>
{
  typedef typename NodeMatrix<N,float>::iterator        iterator;
  typedef typename NodeMatrix<N,float>::const_iterator  const_iterator;

  typedef nvector<double>::iterator        viterator;
  typedef nvector<double>::const_iterator  const_viterator;

public:

  int ncomp() const { return N;}

  void   operator *= (const FMatrixBlock<N>&);
  void   operator *= (double s);

  void   transpose();
  void   transpose(FMatrixBlock<N>& A);
  void   copy_transpose(const FMatrixBlock<N>& A);

  void   zero_row(int);
  void   uno_diag(int);
  float& diag(int i);

  void   DirichletRow (const std::vector<int>& cv);
  void   DirichletCol (const std::vector<int>& cv);
  void   DirichletDiag(const std::vector<int>& cv);
 


  void   entry     (const nmatrix<double>&);
  void   entry     (int i, int j, const EntryMatrix&, double s=1.);
  void   dual_entry(int i, int j, const EntryMatrix&, double s=1.);
  void   inverse ();
  void   vmult   (viterator) const;
  void   mult    (FMatrixBlock<N>&, const FMatrixBlock<N>&) const; 

  void   submult(const FMatrixBlock<N>& B, const FMatrixBlock<N>& C)
  {
    // this -= B*C
    nvector<float>::iterator p(begin());
    for (char i=0; i<N; i++)
      {
	for (char j=0; j<N; j++)
	  {
	    nvector<float>::const_iterator pC(C.begin()+j);
	    nvector<float>::const_iterator pB(B.begin()+i*N);
	    nvector<float>::const_iterator qB(pB+N);
	    //for (int k=0; k<N; k++)
	    for (; pB!=qB; pB++)
	      {
		//value(i,j) -= B(i,k) * C(k,j);
		*p -= *pB * *pC;
		pC += N;
	      }
	    p++;
        }      
      }
  }

  void add(double s, const FMatrixBlock<N>& A)
    {
      for (int i=0; i<N; i++)
	{
	  for (int j=0; j<N; j++)
	    {
	      value(i,j) += s*A(i,j);
	    }
	}
    }
  void adddiag(const nvector<double>& s, double l)
    {
      for (int i=0; i<N; i++)
	{
	  value(i,i) += s[i]*l;
	}
    }
  void add(double s, const TimePattern& TP);

  void cadd(double s, viterator p, const_viterator q0) const
    {
      const_iterator pm = begin();
      const_viterator pend = p+N;
      for ( ; p!=pend; p++)
	{
	  double sum = 0.;
	  const_viterator qend(q0+N);
	  for (const_viterator q=q0; q!=qend; q++)
	    {
	      sum +=  *pm++ * *q;
	    }
	  *p += s*sum;
	}
    }
  void caddtrans(double s, viterator p, const_viterator q0) const
    {
      const_iterator pm = begin();

      for (int k=0; k<N; k++)
	{
	  for (int h=0; h<N; h++)
	    {
	      //sum += M(k,h) * (*(q0+h));
	      *p++ += *pm++ * *q0;
	    }
	  p -= N;
	  q0++;
	}
    }
  void subtract(viterator p0, const_viterator q0) const
    {
      const_iterator pm = begin();

      for (viterator p(p0); p!=p0+N; p++)
	{
	  for (const_viterator q(q0); q!=q0+N; q++)
	    {
	      *p -= *pm++ * *q;
	    }
	}
    };
  std::ostream& print(std::ostream& s) const;

  // Zugriff auf Inhalt ueber ganzen Vektor, damits auch ohne
  // Struktur geht.
  void vector_get(nvector<float>& v) const
    {
      v.resize(size());
      for (int i=0;i<size();++i)
	v[i]=operator[](i);
    }
  void vector_set(nvector<float>& v)
    {
      assert(v.size()==size());
      for (int i=0;i<size();++i)
	operator[](i)=v[i];
    }
  void vector_add(double d, nvector<float>& v)
    {
      assert(v.size()==N*N);
      for (int i=0;i<size();++i)
	operator[](i)+=d*v[i];
    }
};
}

#endif


