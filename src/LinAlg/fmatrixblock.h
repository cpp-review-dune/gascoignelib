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
  void   submult (const FMatrixBlock<N>&, const FMatrixBlock<N>&); 
  void   vmult   (viterator) const;
  void   mult    (FMatrixBlock<N>&, const FMatrixBlock<N>&) const; 

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
      for (int k=0; k<N; k++)
	{
	  double sum = 0.;
	  for (int h=0; h<N; h++)
	    {
	      sum += *pm++ * (*(q0+h));
	    }
	  *p++ += s*sum;
	}
      p -= N;
/*       const_iterator ende = q0+N; */
/*       for (iterator pc = p; pc!=p+N; pc++) */
/* 	{ */
/* 	  double sum = 0.; */
/* 	  for (const_iterator q00 = q0; q00<ende; q00++) */
/* 	    { */
/* 	      sum += *pm++ * *q00; */
/* 	    } */
/* 	  *pc += s*sum; */
/* 	} */
/*       for (int k=0; k<N; k++) */
/* 	{ */
/* 	  double sum = 0.; */
/* 	  for (int h=0; h<N; h++) */
/* 	  { */
/* 	      sum += value(k,h) * (*(q0+h)); */
/* 	  } */
/* 	  *p++ += s*sum; */
/* 	} */
      /*      p -= N; */
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
      q0 -= N;
    }
  void subtract(viterator p, const_viterator q0) const
    {
      const_iterator pm = begin();

      for (int k=0; k<N; k++)
	{
	  for (int h=0; h<N; h++)
	    {
	      *p -= *pm++ * (*(q0+h));
	    }
	  p++;
	}
      p -= N;
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


