#include  "fmatrixblock.h"

namespace Gascoigne
{
static nvector<double>  vhelp(100);  // maximale Blockgroesse !!!

/**********************************************************/

template<int N>
std::ostream& FMatrixBlock<N>::print(std::ostream& s) const
{
  s << *this;
  return s;
}

/**********************************************************/

template<int N>
void  FMatrixBlock<N>::add(double s, const TimePattern& TP)
{
  for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
        {
          value(i,j) += s*TP(i,j);
        }
    }
}

/**********************************************************/

template<int N>
  void   FMatrixBlock<N>::DirichletRow (const std::vector<int>& cv)
{
  for (int i=0; i<cv.size(); i++) {
    int c = cv[i];
    for (int d=0; d<N; d++) value(c,d) = 0.;
  }
}

template<int N>
  void   FMatrixBlock<N>::DirichletCol (const std::vector<int>& cv)
{
  for (int i=0; i<cv.size(); i++) {
    int c = cv[i];
    for (int d=0; d<N; d++) value(d,c) = 0.;
  }
}

template<int N>
  void   FMatrixBlock<N>::DirichletDiag(const std::vector<int>& cv)
{
  for (int i=0; i<cv.size(); i++) {
    int c = cv[i];
    value(c,c) = 1.;
  }
}

template<int N>
void FMatrixBlock<N>::zero_row(int c)
{
  for (int i=0; i<N; i++) value(c,i) =0.;
}
 
template<int N>
void FMatrixBlock<N>::uno_diag(int c)
{
  value(c,c) =1.;
}

/**********************************************************/

template<int N>
float& FMatrixBlock<N>::diag(int i)
{
  return (*this)(i,i);
}

/**********************************************************/

template<int N>
void FMatrixBlock<N>::entry(const nmatrix<double>& E)
{
  for (int c=0; c<N; c++)
    {
      for (int d=0; d<N; d++)
        {
            value(c,d) += E(c,d); 
        }
    }
}

/**********************************************************/

template<int N>
void FMatrixBlock<N>::entry(int i, int j, const EntryMatrix& E, double s)
{
  for (int c=0; c<N; c++)
    {
      for (int d=0; d<N; d++)
        {
            value(c,d) += s * E(i,j,c,d); 
        }
    }
}

/**********************************************************/

template<int N>
void FMatrixBlock<N>::dual_entry(int i, int j, const EntryMatrix& E, double s)
{
  for (int c=0; c<N; c++)
    {
      for (int d=0; d<N; d++)
        {
            value(c,d) += s * E(j,i,d,c); 
        }
    }
}

/**********************************************************/

template<int N>
void FMatrixBlock<N>::transpose() 
{
  for(int i=0;i<N;i++)
    {
      for(int j=0;j<i;j++)
        {
          std::swap(value(i,j),value(j,i));
        }
    }
}

/**********************************************************/

template<int N>
void FMatrixBlock<N>::transpose(FMatrixBlock<N>& A) 
{
  for(int i=0;i<N;i++)
    {
      for(int j=0;j<N;j++)
        {
          std::swap(value(i,j),A(j,i));
        }
    }
}

/**********************************************************/

template<int N>
void FMatrixBlock<N>::copy_transpose(const FMatrixBlock<N>& A) 
{
  for(int i=0;i<N;i++)
    {
      for(int j=0;j<N;j++)
        {
          value(i,j)=A(j,i);
        }
    }
}

/**********************************************************/

template<int N>
void FMatrixBlock<N>::operator *= (const FMatrixBlock<N>& B)
{
  for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
        {
          vhelp[j] = value(i,j);
        }
      for (int j=0; j<N; j++)
        {
          value(i,j) = 0.;
          for (int k=0; k<N; k++)
            {
              value(i,j) += vhelp[k] * B(k,j);
            }
        }      
    }
}

/**********************************************************/

template<int N>
void FMatrixBlock<N>::operator *= (double s)
{
  nvector<float>::iterator       p(begin());
  nvector<float>::const_iterator q(p+N*N);
 
  for ( ; p!=q; p++)
    {
      *p *= s;
    }
}

/**********************************************************/

template<int N>
void FMatrixBlock<N>::inverse()
{ 
  gauss_jordan();
}

/**********************************************************/

template<int N>
void FMatrixBlock<N>::vmult(viterator p) const
{
  
  // copy old entries of vector in vhelp
  //
  nvector<double>::iterator       a = vhelp.begin();
  nvector<double>::const_iterator b = a+N;
  
  for ( ; a!=b; a++)
    {
      *a = *p++;
    }
  p -= N;
  a -= N;

  nvector<float>::const_iterator q(begin());
  for (viterator c=p; c!=p+N; c++)
    {
      *c = 0.;
      for (; a!=b; a++)
	{
	  *c += *q++ * *a;
	} 
      a -= N;
    }
  
}

/**********************************************************/

template FMatrixBlock<1>;
template FMatrixBlock<2>;
template FMatrixBlock<3>;
template FMatrixBlock<4>;
}

/**********************************************************/

