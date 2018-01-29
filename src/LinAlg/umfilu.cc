/**
*
* Copyright (C) 2004, 2007 by the Gascoigne 3D authors
*
* This file is part of Gascoigne 3D
*
* Gascoigne 3D is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version.
*
* Gascoigne 3D is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* Please refer to the file LICENSE.TXT for further information
* on this license.
*
**/


#include  "umfilu.h"
#include  <fstream>
#include <sstream>

#ifdef __WITH_UMFPACK__

extern double DDD;

using namespace std;
 
#define UMFPACK_OK       0
#define UMFPACK_INFO    90
#define UMFPACK_CONTROL 20

#define UMFPACK_A	(0)	/* Ax=b		*/
#define UMFPACK_At	(1)	/* A'x=b	*/

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/

namespace Gascoigne
{
extern "C" int umfpack_di_symbolic
(
    int n,
    int m,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;
extern "C" int umfpack_di_numeric
(
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    void *Symbolic,
    void **Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;
extern "C" int umfpack_di_solve
(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    double X [ ],
    const double B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;
extern "C" void umfpack_di_free_symbolic
(
    void **Symbolic
) ;
extern "C" void umfpack_di_free_numeric
(
    void **Numeric
) ;

extern "C" int umfpack_triplet_to_col
(
    int n,
    int nz,
    const int Ti [ ],
    const int Tj [ ],
    const double Tx [ ],
    int Bp [ ],
    int Bi [ ],
    double Bx [ ]
) ;

extern "C" void umfpack_di_report_status
(
    const double Control [UMFPACK_CONTROL],
    int status
) ;
extern "C" void umfpack_di_report_info
(
    const double Control [UMFPACK_CONTROL],
    const double Info [UMFPACK_INFO]
) ;
extern "C" void umfpack_di_report_control
(
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" int umfpack_di_report_symbolic
(
    const char name [ ],
    void *Symbolic,
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" int umfpack_di_report_numeric
(
    const char name [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" void umfpack_di_report_control
(
    const double Control [UMFPACK_CONTROL]
) ;
extern "C" void umfpack_di_defaults
(
    const double Control [UMFPACK_CONTROL]
) ;

/* ----------------------------------------- */

UmfIlu::UmfIlu(const MatrixInterface* A) 
  : SimpleMatrix(), Control(NULL), Info(NULL), Symbolic(NULL), Numeric(NULL)
{
  AP = dynamic_cast<const SimpleMatrix*>(A);
  assert(AP);

  Control = new double[UMFPACK_CONTROL];
  umfpack_di_defaults(Control);
  Control[0] = 2;
}

/* ----------------------------------------- */

UmfIlu::~UmfIlu()
{
  umfpack_di_free_symbolic (&Symbolic) ;
  umfpack_di_free_numeric (&Numeric) ;

  if(Control) delete[] Control;
  Control=NULL;
  if(Info) delete[] Info;
  Info=NULL;
}

/*-------------------------------------------------------------*/

  double UmfIlu::exec(std::string cmd) const
  {
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return 0;
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
      if(fgets(buffer, 128, pipe) != NULL)
	result += buffer;
    }
    pclose(pipe);
    return atof(result.c_str());
  }
  
  double UmfIlu::condition_number()  const
  {
    
    const ColumnStencil* SA = dynamic_cast<const ColumnStencil*>(AP->GetStencil());
    assert(SA);
    
    const int* sb = &(*SA->start().begin());
    const int* cb = &(*SA->col().begin());
    const double* mb = &AP->GetValue(0);
    
    
    int N = SA->n();
    assert(N>0);

    //    if (N<100000) return 0.0;
    
    ostringstream ss;
    ss << "mat" << N;
    string datnamebase = ss.str();
    ss << ".dat";
    string datname = ss.str();
    ofstream OUT(datname.c_str());

    
    for (int r=0;r<SA->n();++r)
      {
	for (int p=SA->start(r);p<SA->stop(r);++p)
	  {
	    int c = SA->col(p);

	    double v = AP->GetValue(p);
	    if ((v==0)&&(r!=c)) continue;
	    OUT <<  r+1 <<  "\t" << c+1 << "\t" << v << std::endl;
	  }
      }
    OUT.close();
    ss.clear();
    ss.str("");
    ss << datnamebase << ".m";
    string mname = ss.str();
    OUT.open(mname.c_str());
    
    OUT << "load " << datnamebase << ".dat" << endl;
    OUT << "A = spconvert(" << datnamebase << ")" << endl;
    //    OUT << "condest(A)" << endl;
    //        OUT << "eigs(A,1)" << endl;
    OUT << "eigs(A,1)" << endl;
    OUT << "exit" << endl;
    OUT.close();

    cout << "DDD="<< DDD << endl;
    ss.clear();
    ss.str("");
    ss << "/Applications/MATLAB_R2017a.app/bin/matlab  -nojvm -nodesktop -nodisplay -nosplash < " << mname
       << " test.m | tail -3 | head -1 " ;
    double cond =  exec(ss.str().c_str());

    ss.clear();
    ss.str("");
    ss << "rm " << mname << " " << datname << endl;
    //    assert(system(ss.str().c_str())==0);

    return cond;


  }




  

void UmfIlu::ReInit(const SparseStructureInterface* SS)
{
  const ColumnStencil* SA = dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);
  SimpleMatrix::ReInit(SA->n(),SA->nentries());

  umfpack_di_free_symbolic (&Symbolic) ;

  int n = SA->n();
  const int* sb = &(*SA->start().begin());
  const int* cb = &(*SA->col().begin());

  int status = umfpack_di_symbolic(n, n, sb, cb, NULL, &Symbolic, Control, Info);

  if(status != UMFPACK_OK)
    {
      umfpack_di_report_info(Control,Info);
      umfpack_di_report_status(Control,status);
      ofstream file("MATRIX");
//       AP->Write(file);
      cerr << "umfpack_symbolic failed\n"; exit(1);
    }
}

/*-------------------------------------------------*/

void UmfIlu::copy_entries(const MatrixInterface&  A)
{
}

/*-------------------------------------------------------------*/

void UmfIlu::ConstructStructure(const IntVector& perm, const MatrixInterface& A)
{
}

/*-----------------------------------------*/

void UmfIlu::Factorize()
{
  //      cout << "Konditionszahl   " << condition_number() << endl;
				

  //
  // baue LU auf
  //

  umfpack_di_free_numeric (&Numeric) ;

  const ColumnStencil* SA = dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);

  const int* sb = &(*SA->start().begin());
  const int* cb = &(*SA->col().begin());
  const double* mb = &AP->GetValue(0);
  int status = umfpack_di_numeric(sb, cb, mb, Symbolic, &Numeric, Control, Info) ;
  if(status != UMFPACK_OK)
    {
      umfpack_di_report_info(Control,Info);
      umfpack_di_report_status(Control,status);
      cerr << "umfpack_numeric failed\n"; exit(1);
    }
  //   umfpack_report_numeric("LU von A\n",Numeric,Control);
}

/*-----------------------------------------*/

void UmfIlu::Solve(DoubleVector& x, const DoubleVector& b)
{
  const ColumnStencil* SA = dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);

  const int* sb = &(*SA->start().begin());
  const int* cb = &(*SA->col().begin());
  const double* mb = &AP->GetValue(0);
  double* xb = &(*x.begin());
  const double* bb = &(*b.begin());
  int status = umfpack_di_solve (UMFPACK_At, sb, cb, mb, xb, bb, Numeric, Control, Info) ;

    if(status != UMFPACK_OK)
      {
	umfpack_di_report_info(Control,Info);
	umfpack_di_report_status(Control,status);
	cerr << "umfpack_di_solve failed\n"; exit(1);
      }
}

/*-----------------------------------------*/

void UmfIlu::SolveTranspose(DoubleVector& x, const DoubleVector& b)
{
  const ColumnStencil* SA = dynamic_cast<const ColumnStencil*>(AP->GetStencil());
  assert(SA);

  const int* sb = &(*SA->start().begin());
  const int* cb = &(*SA->col().begin());
  const double* mb = &AP->GetValue(0);
  double* xb = &(*x.begin());
  const double* bb = &(*b.begin());
  int status = umfpack_di_solve (UMFPACK_A, sb, cb, mb, xb, bb, Numeric, Control, Info) ;

    if(status != UMFPACK_OK)
      {
	umfpack_di_report_info(Control,Info);
	umfpack_di_report_status(Control,status);
	cerr << "umfpack_di_solve failed\n"; exit(1);
      }
}
}
 
#undef UMFPACK_OK     
#undef UMFPACK_INFO   
#undef UMFPACK_CONTROL

#undef UMFPACK_A	
#undef UMFPACK_At	

#endif
