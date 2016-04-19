#include "loop.h"
#include "gascoignemesh2d.h"
#include <time.h>
#include  "backup.h"
#include  "stdmultilevelsolver.h"
#include  "sparseblockmatrix.h"
#include  "fmatrixblock.h"
using namespace Gascoigne;
using namespace std;




double TIME, DT, DELTAMIN;

string Loop::PrecondCGMass(GlobalVector& u, GlobalVector& f, const TimePattern& TP, double s)
{
  double global_tol = 1.e-12;
  double tol        = 1.e-8;
  int    maxiter    = 50;
    
  bool reached;
  int iter = 0;
    
  GlobalVector g(u.ncomp(),u.n());
  GlobalVector r(u.ncomp(),u.n());
  GlobalVector d(u.ncomp(),u.n());
    
  assert(u.ncomp()==g.ncomp());
  assert(u.n()==g.n());
    
  MM.PrepareJacobi(s);
    
  MM.vmult_time(f,u,TP,-s);

  MM.JacobiVector(u);
  MM.Jacobi(f);
    
  r.equ(1,f);
  d.equ(1,f);
  double Res = r*r;
  double FirstRes = Res;

    
  if (sqrt(Res)<global_tol)
    {
      reached = true;
    }
  else
    {
      reached = false;
    }
    
  while(!reached && iter<maxiter)
    {
      iter++;
      g.zero();
      MM.vmult_time_Jacobi(g,d,TP,s);
      double lambda = Res/(g*d);
	
      u.add(lambda,d);
      r.add(-lambda,g);
	
      Res = r*r;
      if (Res < tol*tol * FirstRes || sqrt(Res)<global_tol)
	{
	  reached = true;
	}
      double betacg = -(r*g)/(d*g);
      d.sequ(betacg,1.,r);
    }
    
  MM.Jacobi(u);

  cout << "\t PCG: " << iter << "\t" << FirstRes << "\t" << sqrt(Res) << endl << endl;

  if(iter==maxiter)
    {
      return "too many iterations";
    }
  else
    {
      return "converged";
    }
}
  
void Loop::AssembleMassMatrix()
{
   GetMultiLevelSolver()->SetProblem("other");
  SparseStructure SA;
  GetMultiLevelSolver()->GetSolver()->GetDiscretization()->Structure(&SA);
  MM.ReInit(&SA);
  MM.zero();

  K.ReInit(&SA);
  K.zero();
  
  GetMultiLevelSolver()->GetSolver()->GetDiscretization()->MassMatrix(MM);

    const ColumnDiagStencil* ST = dynamic_cast<const ColumnDiagStencil* >(MM.GetStencil());
  assert(ST);
  
  LMM.resize(ST->n());
  LMM.zero();
  
  for (int row = 0; row<ST->n(); ++row)
    {
      for (int pos = ST->start(row); pos<ST->stop(row); ++pos)
	LMM[row] += MM.GetValue(pos);
     
	
      assert(LMM[row]>0);
    }
  
  // C-Matrix aufbauen.
  GetMultiLevelSolver()->SetProblem("test");
 
  VectorInterface f("f"); // nur, weil Gascoigne einen Vektor braucht. Wird aber nicht verwendet.
  GetMultiLevelSolver()->AssembleMatrix(f);



  
  // D-Matrix initialisieren
  D.ReInit(&SA);

  
  

}

void Loop::SolveTransport(double dt, VectorInterface& h, VectorInterface& uh, VectorInterface& ul, VectorInterface& oldh, VectorInterface& oldu, VectorInterface& f) 
  

{
  
  GetMultiLevelSolver()->SetProblem("test");
  
  
  // Loesen von UL
  GlobalVector& UL = GetMultiLevelSolver()->GetSolver()->GetGV(ul);
  GlobalVector& UH = GetMultiLevelSolver()->GetSolver()->GetGV(uh);
  GlobalVector& B  = GetMultiLevelSolver()->GetSolver()->GetGV(f);
  GlobalVector& H  = GetMultiLevelSolver()->GetSolver()->GetGV(h);
  GlobalVector& OLDH  = GetMultiLevelSolver()->GetSolver()->GetGV(oldh);
  GlobalVector& OLDU = GetMultiLevelSolver()->GetSolver()->GetGV(oldu);
  B.zero();
  //Loesung UL
    
  // C-Matrix
  const StdSolver* S = dynamic_cast<const StdSolver*> (GetMultiLevelSolver()->GetSolver());
  assert(S);
  const SparseBlockMatrix<FMatrixBlock<2> >* C =
    dynamic_cast<const SparseBlockMatrix<FMatrixBlock<2> >* > (S->GetMatrix());
  assert(C);
  
  
  
  // D-Matrix aufbauen
  D.zero();
  // B = operator K
  const ColumnDiagStencil* ST = dynamic_cast<const ColumnDiagStencil*> (C->GetStencil());
  assert(ST);
  Dii.resize(ST->n());
  Dii.zero();
  for (int row=0;row<ST->n();++row)
    {
      for (int p=ST->start(row);p<ST->stop(row);++p)
	{
	  int col = ST->col(p);
	  
	  const FMatrixBlock<2>& Cij = *(C->mat(p));
	  double CijVj = Cij(0,0)*OLDU(col,0) + Cij(1,1)*OLDU(col,1);
	  
	  	  // suche index (j,i)
	  int pt=ST->start(col);
	  for (;pt<ST->stop(col);++pt)
	    if (ST->col(pt)==row) break;
	  assert(pt<ST->stop(col));
	 	 
 
	  const FMatrixBlock<2>& Cji = *(C->mat(pt));
	  double CjiVi = Cji(0,0)*OLDU(row,0) + Cji(1,1)*OLDU(row,1);
	  double Cji0=Cji(0,0);
	  double Cji1=Cji(1,1);
	  double eps= 1.e-8;
	  
	  if(row!=col)
	  {
	    
	    
	      D.GetValue(p) = max(fabs(CijVj),fabs(CjiVi));
	       Dii[row]+=D.GetValue(p);
	       
	       /*  if( D.GetValue(p)<= eps)
		{
		  D.GetValue(p)=(D.GetValue(p)*D.GetValue(p)+eps*eps)/(2.0*eps);
	  
		Dii[row]+=D.GetValue(p);
	         }
	       */
               }

    }
    } 
    
  
  // D Diagonale
  
  for (int row=0;row<ST->n();++row)
    {      
      int p = ST->diag(row);
      D.GetValue(p) = -Dii[row];
    }
  
  
  
   
  
	
	GlobalVector dd(UL.ncomp(),UL.n());
	assert(UL.ncomp()==dd.ncomp());
	dd.zero();

	
	GlobalVector kk(UL.ncomp(),UL.n());
	assert(UL.ncomp()==kk.ncomp());
	kk.zero();


	
	
	
	TimePattern TP(2); TP(0,0)=1.0; TP(1,1) = 1.0;
	
  // (Ku)i= (-cij vj) ui
	
	for (int row=0;row<ST->n();++row)
	  {
	    for (int p=ST->start(row);p<ST->stop(row);++p)
	      {
		int col = ST->col(p);
		const FMatrixBlock<2>& Cij = *(C->mat(p));
		//  rr(row,0) +=(Cij(0,0)+Cij(1,1));
		// rr(row,1) += OLDU(col,0)+OLDU(col,1);
		 kk(row,0) += (-Cij(0,0)*OLDU(col,0)-Cij(1,1)*OLDU(col,1)) *  OLDH(col,0);
		 kk(row,1) += (-Cij(0,0)*OLDU(col,0)-Cij(1,1)*OLDU(col,1)) *  OLDH(col,1);
		
		 //	B(row,0) += (-Cij(0,0)*OLDU(col,0)-Cij(1,1)*OLDU(col,1)) *  OLDH(col,0);
		 //	B(row,1) += (-Cij(0,0)*OLDU(col,0)-Cij(1,1)*OLDU(col,1)) *  OLDH(col,1);
	      }
	  }
	
	
	// + Du
	// D.vmult_time(B,OLDH,TP,1.0);
	
	
	//  das waere Matrix-Vector Komponentenweise
	for (int row=0;row<ST->n();++row)
	  {
	    for (int p=ST->start(row);p<ST->stop(row);++p)
	      {
		int col = ST->col(p);
		kk(row,0) += D.GetValue(p) * OLDH(col,0);
		kk(row,1) += D.GetValue(p) * OLDH(col,1);
		
		B(row,0) += D.GetValue(p) * OLDH(col,0);
		B(row,1) += D.GetValue(p) * OLDH(col,1);

		
		//	dd(row,0) += D.GetValue(p) * OLDH(col,0);
		//	dd(row,1) += D.GetValue(p) * OLDH(col,1);
		
	      }
	  }
	
	
	// UL mit exp Euler, zwischenloesung fuer Heun
	assert(UL.n() == LMM.size());
	for (int i=0;i<UL.n();++i)
	  {
	    for (int c=0;c<UL.ncomp();++c)
	      UL(i,c)+=(kk(i,c)*dt)/LMM[i];
	  }
	

	//heun verfahren
	for (int row=0;row<ST->n();++row)
	  {
	    for (int p=ST->start(row);p<ST->stop(row);++p)
	      {
		int col = ST->col(p);
		const FMatrixBlock<2>& Cij = *(C->mat(p));
		//  rr(row,0) +=(Cij(0,0)+Cij(1,1));
		// rr(row,1) += OLDU(col,0)+OLDU(col,1);
		
		kk(row,0) += ((-Cij(0,0)*OLDU(col,0)-Cij(1,1)*OLDU(col,1))+D.GetValue(p)) *  UL(col,0);
		kk(row,1) += ((-Cij(0,0)*OLDU(col,0)-Cij(1,1)*OLDU(col,1))+D.GetValue(p)) *  UL(col,1);
	      }
	  }
	//	GetMultiLevelSolver()->GetSolver()->Visu("Results/b",f,_iter);

       
	
	/*
       /// rechte seite von Heun verfahren
       // B = B + (K+D)UL
       // Gascoigne braucht nun statt oldh UL
       // erster Schritt: B = B + K UL
       GetMultiLevelSolver()->GetSolver()->AddNodeVector("oldh",ul);
       GetMultiLevelSolver()->GetSolver()->AddNodeVector("oldu",oldu);
       GetMultiLevelSolver()->GetSolver()->Rhs(f);
       GetMultiLevelSolver()->GetSolver()->DeleteNodeVector("oldh");
       GetMultiLevelSolver()->GetSolver()->DeleteNodeVector("oldu");
       // zweiter SChritt B = B + D UL
       D.vmult_time(B,UL,TP,1.0);
       // B = k/2 B
       */
	kk *= 1.0/2.0;
	//	rr*=1.0/2.0;
	// M UL = B
	for (int i=0;i<UL.n();++i)
	  for (int c=0;c<UL.ncomp();++c)
	    UL(i,c)=(kk(i,c)*dt)/LMM[i];
	// UL = UL + B, jetzt ist UL die Heun-Loesung
	UL+=OLDH;
	
	
	
	
	
	
	
	
  // Flux
	
	GlobalVector FLUX(UL.ncomp(),UL.n());
	assert(UL.ncomp()==FLUX.ncomp());
	FLUX.zero();
	
	for (int row = 0; row<ST->n(); ++row)
	  {
	    for (int c=0;c<UH.ncomp();++c)
	      {
		for (int pos = ST->start(row); pos<ST->stop(row); ++pos)
		  {
		    
		    int col = ST->col(pos);
		    
		    double fij= MM.GetValue(pos)*(UH(row,c)-UH(col,c))+D.GetValue(pos)*(UL(row,c)-UL(col,c));
		    
		    FLUX(row,c) += fij;								    
		  }
	      }
	  }
	for (int i=0;i<UL.n();++i)
	  for (int c=0;c<UL.ncomp();++c)
	    UL(i,c)=UL(i,c)+dt*FLUX(i,c)/LMM[i];
	
	
   //U^n+1
   for (int i=0;i<UL.n();++i)
     for (int c=0;c<UL.ncomp();++c)
       H(i,c)=UL(i,c);
 

   // GetMultiLevelSolver()->GetSolver()->Visu("Results/ul",ul,_iter);
}



void Loop::run(const std::string& problemlabel)
{
  TIME=0.0;
  double tref;
  double dtmax, endtime;
  int prerefine;
  if (1)
    {
      DataFormatHandler DFH;
      DFH.insert("dt",   &DT, 0.);
      DFH.insert("dtmax",   &dtmax, 0.);
      DFH.insert("time", &TIME, 0.);
      DFH.insert("endtime", &endtime, 0.);
      FileScanner FS(DFH);
      FS.NoComplain();
      FS.readfile(_paramfile, "Loop");
      assert(DT>0.0);
    }
  if (1)
    {
      DataFormatHandler DFH;
      DFH.insert("prerefine",   &prerefine, 0);
      FileScanner FS(DFH);
      FS.NoComplain();
      FS.readfile(_paramfile, "Mesh");
    }
  if (1)
    {
      DataFormatHandler DFH;
      DFH.insert("deltamin",   &DELTAMIN, 0.);
      DFH.insert("Tref", &tref, 0.);
      FileScanner FS(DFH);
      FS.NoComplain();
      FS.readfile(_paramfile, "Equation");
      assert(DELTAMIN>0.0);
      assert(tref>0.0);
    }

 
 
  
  // vectors for solution and right hand side
  VectorInterface u("u"), h("h"), f("f"), oldu("oldu"),oldoldu("oldoldu"), oldh("oldh"), uh("uh"), ul("ul"),
    extu("extu"), other("other"),w("w");



  

  

  PrintMeshInformation();
  // initialize problem, solver and vectors 
  GetMultiLevelSolver()->ReInit("test");
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(h);
  GetMultiLevelSolver()->ReInitVector(oldu);
  GetMultiLevelSolver()->ReInitVector(oldoldu);
  GetMultiLevelSolver()->ReInitVector(oldh);
  GetMultiLevelSolver()->ReInitVector(extu);
  GetMultiLevelSolver()->ReInitVector(f);
  GetMultiLevelSolver()->ReInitVector(w);
  GetMultiLevelSolver()->ReInitVector(uh);
  GetMultiLevelSolver()->ReInitVector(ul);
    
  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();

  GetMultiLevelSolver()->ReInitVector(other);

  
  cout << "------------------------------" << endl;
  cout << "See-Eis" << endl;
  cout << "------------------------------" << endl << endl;
  
  // See-Eis
  GetMultiLevelSolver()->SetProblem("seaice");
  InitSolution(u);
  InitSolution(oldu);
  
  // Transport
  GetMultiLevelSolver()->SetProblem("transport");
  
  InitSolution(h);
  InitSolution(oldh);


  
  GetMultiLevelSolver()->SetProblem("seaice");
  // dimitri
  AssembleMassMatrix();


  
 


  // set up system matrix

  double writeeveryhour=24.0; // 1.0;

  double writenext= writeeveryhour;
  int writeiter = 0;
  string res;
  
  
  nvector<double> functionals;
  stringstream str;
  str << "func_" << prerefine << "_" << dtmax << ".txt";
  ofstream OUTF(str.str().c_str());

  
    
  int timeinc=0;
  
  clock_t start, end;
  double cpu_time_used;

      

      


  GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,writeiter);
  GetMultiLevelSolver()->GetSolver()->Visu("Results/h",h,writeiter);
  // GetMultiLevelSolver()->GetSolver()->Visu("Results/o",other,writeiter);
  
 
  
  
  for (_iter=1; _iter<=_niter; _iter++)
    {
      if ((TIME>10)&&(timeinc==0)) { DT*=2.0; ++timeinc; }
      if ((TIME>20)&&(timeinc==1)) { DT*=2.0; ++timeinc; }
      //  if ((TIME>30)&&(timeinc==2)) { DT*=2.0; ++timeinc; }
	//   if ((TIME>40)&&(timeinc==2)) { DT*=2.0; ++timeinc; }
     
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

      cout << "\n Zeitschritt " << _iter << " " << TIME << " -> " 
	   << TIME+DT << "  [" << DT << "]\t"  
	   << "\t" << TIME*tref/60/60/24 << " Tage" << endl;
      TIME += DT;
      if (TIME>endtime*24*60*60/tref)
	{
	  DT = endtime*24*60*60/tref - (TIME-DT);
	  TIME = endtime*24*60*60/tref;
		 
 	  cout << "\n Letzter Zeitschritt " << _iter << " " << TIME 
	       << "\t" << TIME*tref/60/60/24 << " Tage" << endl;
	}
      
      // Extrapolation
      GetMultiLevelSolver()->GetSolver()->Equ(extu,1.0+DT/DT, u);
      GetMultiLevelSolver()->GetSolver()->Add(extu,-DT/DT, oldu);
      GetMultiLevelSolver()->Equ(oldu,1.0,u);
      GetMultiLevelSolver()->Equ(oldoldu,1.0,oldu);
      GetMultiLevelSolver()->Equ(u,1.0,extu); // besserer Startwert fuer newton

      
       /// Masslumping
      GetMultiLevelSolver()->Equ(oldh,1.0,h);
      
      GetMultiLevelSolver()->SetProblem("test");


      int NSUB = 800;

     

      for (int sub=0;sub<NSUB;++sub)
	{
	 double dt=DT/NSUB;
	  GetMultiLevelSolver()->GetSolver()->Equ(u,( (double) (NSUB-sub))/ ( (double) NSUB),oldu); // besserer Startwert fuer newton
	  GetMultiLevelSolver()->GetSolver()->Add(u,( (double) (sub))/ ( (double) NSUB),extu);      // besserer Startwert fuer newton
      
	  SolveTransport(dt, h,uh,ul,oldh,oldu,f);
	 
	
	}
      
      GlobalVector& H = GetMultiLevelSolver()->GetSolver()->GetGV(h);
      for (int i=0;i<H.n();++i)
	{
	  H(i,1) = std::min(1.0, H(i,1));
	  H(i,1) = std::max(0.0, H(i,1));
	  H(i,0) = std::max(0.0, H(i,0));
	}
   
      
      
      // Eis-Problem
      
        GetMultiLevelSolver()->SetProblem("seaice");
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      GetMultiLevelSolver()->AddNodeVector("oldu", oldu);
      GetMultiLevelSolver()->AddNodeVector("oldoldu", oldoldu);
      GetMultiLevelSolver()->AddNodeVector("extu", extu);
      GetMultiLevelSolver()->AddNodeVector("H", h);
      GetMultiLevelSolver()->AddNodeVector("W", w);

      start = clock();
      //    cout <<"Solve Seaice"<< endl;
       res = Solve(u,f,"Results/u");
         assert(res == "converged");
      
      end = clock();
      
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
      //   cout<< " cpu_time_used"<<  cpu_time_used << endl;
      //cout<< " start"<<  start << "end"<< end << endl;
     
      //   functionals = Functionals(u,f);
     
      GetMultiLevelSolver()->DeleteNodeVector("oldu");
      GetMultiLevelSolver()->DeleteNodeVector("extu");
      GetMultiLevelSolver()->DeleteNodeVector("H");
      GetMultiLevelSolver()->DeleteNodeVector("oldoldu");
      GetMultiLevelSolver()->DeleteNodeVector("W");
      
      
      
      // Other-Problem
      if (0)
	{
	  GetMultiLevelSolver()->SetProblem("other");
	  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
	  GetMultiLevelSolver()->AddNodeVector("U", u);
	  GetMultiLevelSolver()->AddNodeVector("H", h);
	  Solve(other,f,"Results/o");
	  GetMultiLevelSolver()->DeleteNodeVector("U");
	  GetMultiLevelSolver()->DeleteNodeVector("H");
	}
      
      


      //  OUTF << TIME << " " << functionals <<  cpu_time_used <<endl;
      
      
      double time_in_hours = TIME*tref/60/60;
      if (time_in_hours>writenext)
	{
	  
	  writenext += writeeveryhour;
	  ++writeiter;
	  
	   GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,writeiter);
	  GetMultiLevelSolver()->GetSolver()->Visu("Results/h",h,writeiter);
	  //  GetMultiLevelSolver()->GetSolver()->Visu("Results/o",other,writeiter);
	}
      
      
      if (TIME >= endtime*24*60*60/tref) break;
      
    }
  OUTF.close();
 }


void Loop::runwater(const std::string& problemlabel)
{
  VectorInterface f("f"),  w("w");

  PrintMeshInformation();
  // initialize problem, solver and vectors 
  GetMultiLevelSolver()->ReInit(problemlabel);
  GetMultiLevelSolver()->ReInitVector(w);
  GetMultiLevelSolver()->ReInitVector(f);

  if (0)
    for (int i=0;i<3;++i)
      {
      
	const GascoigneMesh2d* M = dynamic_cast<const GascoigneMesh2d*> (GetMultiLevelSolver()->GetSolver()->GetMesh());
	assert(M);

	const nvector<int>& ref5 = *(M->VertexOnBoundary(5));
	const nvector<int>& ref6 = *(M->VertexOnBoundary(6));
	nvector<int> ref;
	for (int n=0;n<ref5.size();++n) ref.push_back(ref5[n]);
	for (int n=0;n<ref6.size();++n) ref.push_back(ref6[n]);

	GetMeshAgent()->refine_nodes(ref);

      
	GetMultiLevelSolver()->ReInit(problemlabel);
	GetMultiLevelSolver()->ReInitVector(w);
	GetMultiLevelSolver()->ReInitVector(f);
      }




  cout << "------------------------------" << endl;
  cout << "Wasser" << endl;
  cout << "------------------------------" << endl << endl;

  //WaterProblem();
  GetMultiLevelSolver()->SetProblem("water");
  // wasserproblem loesen
  Solve(w,f);
  GetMultiLevelSolver()->GetSolver()->Visu("water",w,0);

  const GlobalVector& W = GetMultiLevelSolver()->GetSolver()->GetGV(w);
  GlobalVector Wohnedruck; Wohnedruck.ncomp()=2; Wohnedruck.resize(W.n());
  Wohnedruck.CompEq(0,1.0,1,W);
  Wohnedruck.CompEq(1,1.0,2,W);
  WriteBackUp(Wohnedruck,"water");
  
}








