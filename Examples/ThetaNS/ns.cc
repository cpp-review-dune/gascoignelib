#include  "ns.h"
#include  "filescanner.h"

extern double __DT;
extern double __THETA;
extern double __TIME;



extern double __DT_OLD;
extern double __THETA_OLD;

extern bool __FIRSTSTEP, __LASTSTEP, __MIDDLESTEP;



using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  NS::NS(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("nu" ,    &__nu , 0.0);
    DFH.insert("lpsp" ,   &__lpsp0 , 0.0);
    DFH.insert("lpsv" ,   &__lpsv0 , 0.0);
    FileScanner FS(DFH, pf, "Equation");
    assert(__nu>0);
    assert(__lpsp0>0);
  }

  void NS::point(double h, const FemFunction& U, const Vertex2d& v) const
  {
    __h = h;
    __lpsp = __lpsp0 / (__nu/ __h / __h + 1.0/__h);
    __lpsv = __lpsv0 / (__nu/ __h / __h + 1.0/__h);
    
  }
  
  
  void NS::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
  {
    // stab
    b[0] += __lpsp * __DT * (U[0].x()*N.x()+U[0].y()*N.y());

    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
	for (int k=0;k<2;++k)
	  b[i+1] += __lpsv * __DT * U[j+1].m()*U[i+1][j+1] * U[k+1].m()*N[k+1];


    // // divergence
    b[0] += __DT * (U[1].x()+U[2].y())*N.m();

    
    
    b[1] -=  __DT * U[0].m() * N.x();
    b[2] -=  __DT * U[0].m() * N.y();
    
    // Time
    //    b[0] +=0.01* (U[0].m() - (*OLD)[0].m()) * N.m();
    for (int c=1;c<3;++c)
      b[c] += (U[c].m() - (*OLD)[c].m()) * N.m();

    
    // impl
    b[1] += __THETA * __DT * __nu * (U[1].x()*N.x()+U[1].y()*N.y());
    b[2] += __THETA * __DT * __nu * (U[2].x()*N.x()+U[2].y()*N.y());
    
    // expl
    b[1] += (1.0-__THETA)*__DT*__nu*((*OLD)[1].x()*N.x()+(*OLD)[1].y()*N.y());
    b[2] += (1.0-__THETA)*__DT*__nu*((*OLD)[2].x()*N.x()+(*OLD)[2].y()*N.y());

    // Convection
    b[1] += __THETA * __DT * (U[1].m() * U[1].x() + U[2].m() * U[1].y()) * N.m();
    b[2] += __THETA * __DT * (U[1].m() * U[2].x() + U[2].m() * U[2].y()) * N.m();

    b[1] += (1.0-__THETA)*__DT*((*OLD)[1].m()*(*OLD)[1].x()+(*OLD)[2].m()*(*OLD)[1].y())*N.m();
    b[2] += (1.0-__THETA)*__DT*((*OLD)[1].m()*(*OLD)[2].x()+(*OLD)[2].m()*(*OLD)[2].y())*N.m();


    
  }
  
  
  void NS::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const TestFunction& N) const
  {    
    double lap = M.x()*N.x()+M.y()*N.y();
    double mass = M.m() * N.m();

    // stab
    A(0,0) += __lpsp * __DT * lap;
    for (int i=0;i<2;++i)
      for (int j=0;j<2;++j)
	for (int k=0;k<2;++k)
	  {
	    A(i+1,j+1) += __lpsv * __DT * M.m()*U[i+1][j+1] * U[k+1].m()*N[k+1];
	    A(i+1,i+1) += __lpsv * __DT * U[j+1].m()*M[j+1] * U[k+1].m()*N[k+1];
	    A(i+1,k+1) += __lpsv * __DT * U[j+1].m()*U[i+1][j+1] * M.m()*N[k+1];
	  }
    




    // divergence
    A(0,1) += __DT * M.x()*N.m();
    A(0,2) += __DT * M.y()*N.m();


    A(1,0) -= __DT * M.m() * N.x();
    A(2,0) -= __DT * M.m() * N.y();
       
    //time
    //    A(0,0)+=0.01* M.m()*N.m();
    for (int c=1;c<3;++c)
      A(c,c) += mass;

    // impl
    A(1,1) += __THETA * __DT * __nu * lap;
    A(2,2) += __THETA * __DT * __nu * lap;


    // Convection
    A(1,1) += __THETA * __DT * (U[1].m() * M.x() + U[2].m() * M.y()) * N.m();
    A(2,2) += __THETA * __DT * (U[1].m() * M.x() + U[2].m() * M.y()) * N.m();

    A(1,1) += __THETA * __DT * M.m() * U[1].x() * N.m();
    A(2,2) += __THETA * __DT * M.m() * U[2].y() * N.m();
    A(1,2) += __THETA * __DT * M.m() * U[1].y() * N.m();
    A(2,1) += __THETA * __DT * M.m() * U[2].x() * N.m();
    
  }



 
 void NS::point_M(int j, const FemFunction& U, const TestFunction& M) const
  {
  }

  
  void NS::MatrixBlock(EntryMatrix& A, const FemFunction& U, const FemFunction& N) const
  {
    for (int j=0; j<N.size();++j)  // trial
      {
#define M N[j]
	point_M(j,U,M);

	for (int i=0; i<N.size();++i) // test
	  {
	    A.SetDofIndex(i,j);
	    Matrix(A,U,M , N[i]);
	  }
#undef M
      }
  }


  void NS::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
    {
      //      b[0] += __DT * __lps * (UP[0].x() * N.x() + UP[0].y() * N.y());
    }
    
  void NS::StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
    {
      //      A(0,0) += __DT * __lps * (Mp.x() * Np.x() + Mp.y() * Np.y());
    }


  ////////////////////////////////////////////////// ADJOINT





  
  AdjointNS::AdjointNS(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("nu" ,    &__nu , 0.0);
    DFH.insert("lpsp" ,   &__lpsp0 , 0.0);
    DFH.insert("lpsv" ,   &__lpsv0 , 0.0);
    FileScanner FS(DFH, pf, "Equation");
    assert(__nu>0);
    assert(__lpsp0>0);
  }



  void AdjointNS::point(double h, const FemFunction& U, const Vertex2d& v) const
  {
    __h = h;
    __lpsp = __lpsp0 / (__nu/ __h / __h + 1.0/__h);
    __lpsv = __lpsv0 / (__nu/ __h / __h + 1.0/__h);
  }
  
  
  void AdjointNS::Form(VectorIterator b, const FemFunction& Z, const TestFunction& N) const
  {
    // time new
    //    b[0] += 0.01*Z[0].m() * N.m();
    for (int c=1;c<3;++c)
      b[c] += Z[c].m() * N.m();

    if (__FIRSTSTEP||__MIDDLESTEP)
      {
	// stab
	b[0] += __lpsp * __DT * (Z[0].x()*N.x()+Z[0].y()*N.y());

	for (int i=0;i<2;++i)
	  for (int j=0;j<2;++j)
	    for (int k=0;k<2;++k)
	      {
		b[j+1] += __lpsv * __DT * N.m()*(*_U)[i+1][j+1] * (*_U)[k+1].m()*Z[i+1][k+1];
		b[i+1] += __lpsv * __DT * (*_U)[j+1].m()*N[j+1] * (*_U)[k+1].m()*Z[i+1][k+1];
		b[k+1] += __lpsv * __DT * (*_U)[j+1].m()*(*_U)[i+1][j+1] * N.m()*Z[i+1][k+1];	
	      }


	// divergence
	b[0] -=  __DT * N.m() * Z[1].x();
	b[0] -=  __DT * N.m() * Z[2].y();
	b[1] += __DT * N.x() * Z[0].m();
	b[2] += __DT * N.y() * Z[0].m();
	
	// new tensor
	b[1] += __THETA * __DT * __nu * (Z[1].x()*N.x()+Z[1].y()*N.y());
	b[2] += __THETA * __DT * __nu * (Z[2].x()*N.x()+Z[2].y()*N.y());

	// new convection
	b[1] += __THETA*__DT*((*_U)[1].m()*N.x()+(*_U)[2].m()*N.y())*Z[1].m();
	b[2] += __THETA*__DT*((*_U)[1].m()*N.x()+(*_U)[2].m()*N.y())*Z[2].m();
	
	b[1] += __THETA*__DT*((*_U)[1].x()*Z[1].m()+(*_U)[2].x()*Z[2].m())*N.m();
	b[2] += __THETA*__DT*((*_U)[1].y()*Z[1].m()+(*_U)[2].y()*Z[2].m())*N.m();
      }

    if (__LASTSTEP||__MIDDLESTEP)
      {
	// divergence
	// b[1] += (1.0-__THETA_OLD) * __DT_OLD * N.x() * (*ZOLD)[0].m();
	// b[2] += (1.0-__THETA_OLD) * __DT_OLD * N.y() * (*ZOLD)[0].m();

	// time old
	//	b[0] += -0.01*(*ZOLD)[0].m() * N.m();
	for (int c=1;c<3;++c)
	  b[c] += -(*ZOLD)[c].m() * N.m();

	// old tensor
	b[1] += (1.0-__THETA_OLD) * __DT_OLD * __nu * ((*ZOLD)[1].x()*N.x()+(*ZOLD)[1].y()*N.y());
	b[2] += (1.0-__THETA_OLD) * __DT_OLD * __nu * ((*ZOLD)[2].x()*N.x()+(*ZOLD)[2].y()*N.y());


	// old convection
	b[1] += (1.0-__THETA_OLD)*__DT_OLD*((*_U)[1].m()*N.x()+(*_U)[2].m()*N.y())*(*ZOLD)[1].m();
	b[2] += (1.0-__THETA_OLD)*__DT_OLD*((*_U)[1].m()*N.x()+(*_U)[2].m()*N.y())*(*ZOLD)[2].m();

	b[1] += (1.0-__THETA)*__DT_OLD*((*_U)[1].x()*(*ZOLD)[1].m()+(*_U)[2].x()*(*ZOLD)[2].m())*N.m();
	b[2] += (1.0-__THETA)*__DT_OLD*((*_U)[1].y()*(*ZOLD)[1].m()+(*_U)[2].y()*(*ZOLD)[2].m())*N.m();
      }

    if (__LASTSTEP)
      b[0] += Z[0].m()*N.m();
  }
  
  
  void AdjointNS::Matrix(EntryMatrix& A, const FemFunction& Z, const TestFunction& M, const TestFunction& N) const
  {    
    double lap = M.x()*N.x()+M.y()*N.y();
    double mas = M.m() * N.m();

    // time new
    //    A(0,0) += 0.01*mas;
    for (int c=1;c<3;++c)
      A(c,c) += mas;
	
    if (__FIRSTSTEP||__MIDDLESTEP)
      {
	// stab
	A(0,0) += __lpsp * __DT * (M.x()*N.x()+M.y()*N.y());
	for (int i=0;i<2;++i)
	  for (int j=0;j<2;++j)
	    for (int k=0;k<2;++k)
	      {
		A(j+1,i+1) += __lpsv * __DT * N.m()*(*_U)[i+1][j+1] * (*_U)[k+1].m()*M[k+1];
		A(i+1,i+1) += __lpsv * __DT * (*_U)[j+1].m()*N[j+1] * (*_U)[k+1].m()*M[k+1];
		A(k+1,i+1) += __lpsv * __DT * (*_U)[j+1].m()*(*_U)[i+1][j+1] * N.m()*M[k+1];
	      }

	// divergence
	A(0,1) -= __DT * N.m() * M.x();
	A(0,2) -= __DT * N.m() * M.y();
	A(1,0) += __DT * N.x() * M.m();
	A(2,0) += __DT * N.y() * M.m();

	// new tensor
	A(1,1) += __THETA * __DT * __nu * lap;
	A(2,2) += __THETA * __DT * __nu * lap;


	// new convection
	A(1,1) += __THETA*__DT*((*_U)[1].m()*N.x()+(*_U)[2].m()*N.y())*M.m();
	A(2,2) += __THETA*__DT*((*_U)[1].m()*N.x()+(*_U)[2].m()*N.y())*M.m();
	
	A(1,1) += __THETA*__DT*(*_U)[1].x()*M.m()*N.m();
	A(2,1) += __THETA*__DT*(*_U)[1].y()*M.m()*N.m();
	A(1,2) += __THETA*__DT*(*_U)[2].x()*M.m()*N.m();
	A(2,2) += __THETA*__DT*(*_U)[2].y()*M.m()*N.m();
      }

    if (__LASTSTEP)
      A(0,0)+= M.m()*N.m();

  }



 
 void AdjointNS::point_M(int j, const FemFunction& U, const TestFunction& M) const
  {
  }

  
  void AdjointNS::MatrixBlock(EntryMatrix& A, const FemFunction& U, const FemFunction& N) const
  {
    for (int j=0; j<N.size();++j)  // trial
      {
#define M N[j]
	point_M(j,U,M);

	for (int i=0; i<N.size();++i) // test
	  {
	    A.SetDofIndex(i,j);
	    Matrix(A,U,M , N[i]);
	  }
#undef M
      }
  }

  
  void AdjointNS::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP, const TestFunction& N) const
  {
    //b[0] += __DT * __lps * (UP[0].x() * N.x() + UP[0].y() * N.y());
  } 
  void AdjointNS::StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np, const TestFunction& Mp) const
  {
    //A(0,0) += __DT * __lps * (Mp.x() * Np.x() + Mp.y() * Np.y());
  }
 
  
}

/*-----------------------------------------*/
