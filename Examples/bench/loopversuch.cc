#include "loop.h"
#include "gascoignemesh2d.h"
#include <time.h>
#include  "backup.h"
#include  "stdmultilevelsolver.h"
#include  "sparseblockmatrix.h"
#include  "fmatrixblock.h"
using namespace Gascoigne;
using namespace std;




double TIME, DT, DTSUB, DELTAMIN,DTUL;

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
  SparseStructure SA;
  GetMultiLevelSolver()->GetSolver()->GetDiscretization()->Structure(&SA);
  MM.ReInit(&SA);
  MM.zero();
  
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
}

void Loop::SolveDIV(VectorInterface& div,VectorInterface& vel,VectorInterface& f)
{
  // Rechte Seite
  GetMultiLevelSolver()->SetProblem("div");
  GetMultiLevelSolver()->Zero(f);
  GetMultiLevelSolver()->GetSolver()->AddNodeVector("vel",vel);
  GetMultiLevelSolver()->GetSolver()->Rhs(f);
  GetMultiLevelSolver()->GetSolver()->DeleteNodeVector("vel");
  
  // mit LMM loesen
  GlobalVector& DIV    = GetMultiLevelSolver()->GetSolver()->GetGV(div);
  GlobalVector& F      = GetMultiLevelSolver()->GetSolver()->GetGV(f);
  for (int i=0;i<DIV.n();++i)
    DIV(i,0) = F(i,0)/LMM[i];
}
  

void Loop::SolveTransport(VectorInterface& h,
			  VectorInterface& div,
			  VectorInterface& vel,
			  VectorInterface& tmp,
			  VectorInterface& f)
{
  ///// Low Order Loesung

  // Rechte Seite
  GetMultiLevelSolver()->SetProblem("tg");
  GetMultiLevelSolver()->Zero(f);
  GetMultiLevelSolver()->GetSolver()->AddNodeVector("oldh",h); 
  GetMultiLevelSolver()->GetSolver()->AddNodeVector("vel",vel);
  GetMultiLevelSolver()->GetSolver()->Rhs(f);
  GetMultiLevelSolver()->GetSolver()->DeleteNodeVector("oldh");
  GetMultiLevelSolver()->GetSolver()->DeleteNodeVector("vel");
  
  // Low-Order Loesung
  GlobalVector& F    = GetMultiLevelSolver()->GetSolver()->GetGV(f);
  GlobalVector& H    = GetMultiLevelSolver()->GetSolver()->GetGV(h);
  GlobalVector& DIV  = GetMultiLevelSolver()->GetSolver()->GetGV(div);
   for (int c=0;c<H.ncomp();++c)
  for (int i=0;i<H.n();++i)
    //  H(i,c) = F(i,c) / ( LMM[i] * (1.0 + DTSUB * std::max(0.0,DIV(i,0)) ) );
    H(i,c) = F(i,c) / ( LMM[i] );

  // Low-Order Loesung
  GlobalVector& UL  = GetMultiLevelSolver()->GetSolver()->GetGV(tmp);
  UL.zero();

 for (int c=0;c<H.ncomp();++c)
  for (int i=0;i<H.n();++i)
    
    UL(i,c)= H(i,c);


  // High-Order Loesung
  GlobalVector& UH  = GetMultiLevelSolver()->GetSolver()->GetGV(tmp);
  UH.zero();
  TimePattern TP(2,2); TP.zero(); TP(0,0)=1.0; TP(1,1) = 1.0;
  PrecondCGMass(UH, F, TP, 1.0);
  // for (int c=0;c<H.ncomp();++c)
    // for (int i=0;i<H.n();++i)
    //    H(i,c)=UH(i,c);

  // fluxidux
  const ColumnDiagStencil* ST = dynamic_cast<const ColumnDiagStencil* >(MM.GetStencil());
  assert(ST);
  
  nvector<double> maxi(2),mini(2), Pminus(2), Pplus(2), Qminus(2), Qplus(2);
  GlobalVector Rplus(H.ncomp(),H.n()), Rminus(H.ncomp(),H.n());
  for (int row = 0; row<ST->n(); ++row)
    {
      // max und min von low-order , Pppus, Pminus
      mini[0]=maxi[0] = H(row,0);
      mini[1]=maxi[1] = H(row,1);
      Pminus=-1.e-10; Pplus = 1.e-10;
      for (int pos = ST->start(row); pos<ST->stop(row); ++pos)
	{
	  int col = ST->col(pos);
	  mini[0] = min(mini[0],H(col,0));  mini[1] = min(mini[1],H(col,1));
	  maxi[0] = max(maxi[0],H(col,0));  maxi[1] = max(maxi[1],H(col,1));

	  Pplus[0]  += max(0.0, MM.GetValue(pos) * (UH(row,0) - UH(col,0)));
	  Pminus[0] += min(0.0, MM.GetValue(pos) * (UH(row,0) - UH(col,0)));

	  Pplus[1]  += max(0.0, MM.GetValue(pos) * (UH(row,1) - UH(col,1)));
	  Pminus[1] += min(0.0, MM.GetValue(pos) * (UH(row,1) - UH(col,1)));
	}

      
      // Q + R bestimmen
      Qplus[0]  = LMM[row] * (maxi[0]-H(row,0));
      Qplus[1]  = LMM[row] * (maxi[1]-H(row,1));
      Qminus[0] = LMM[row] * (mini[0]-H(row,0));
      Qminus[1] = LMM[row] * (mini[1]-H(row,1));
      
      Rplus(row,0) = min(1.0, Qplus[0]/Pplus[0]);
      Rplus(row,1) = min(1.0, Qplus[1]/Pplus[1]);
      Rminus(row,0) = min(1.0, Qminus[0]/Pminus[0]);
      Rminus(row,1) = min(1.0, Qminus[1]/Pminus[1]);
    }

  // Korrektur
    for (int row = 0; row<ST->n(); ++row)
    {
      for (int pos = ST->start(row); pos<ST->stop(row); ++pos)
	{
	  int col = ST->col(pos);
	  double fij0 = MM.GetValue(pos)*(UH(row,0)-UH(col,0));
	  double fij1 = MM.GetValue(pos)*(UH(row,1)-UH(col,1));
	  if(fij0*(UL(col,0)-UL(row,0))>0)
	    {
	    cout<<fij0<<"positiv0"<<endl;
	  fij0=0;
	    }
	  if(fij1*(UL(col,1)-UL(row,1))>0)
	    {   cout<<fij1<<"positiv1"<<endl;
	    fij1=0;
	    }
	  if (fij0 >= 0) 
	    {
	    H(row,0) += 1.0/LMM[row] * min(Rplus(row,0), Rminus(col,0)) * fij0;
	   }
	  else 
	   
	    H(row,0) += 1.0/LMM[row] * min(Rminus(row,0), Rplus(col,0)) * fij0;

	  if (fij1 >=0) 
	    H(row,1) += 1.0/LMM[row] * min(Rplus(row,1), Rminus(col,1)) * fij1;
	  else 
	    H(row,1) += 1.0/LMM[row] * min(Rminus(row,1), Rplus(col,1)) * fij1;
	 
	}
    }

    if(
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
  VectorInterface u("u"), h("h"), f("f"), oldu("oldu"), tmp("tmp"), extu("extu"), div("div"), 



  


  PrintMeshInformation();
  // initialize problem, solver and vectors 
  GetMultiLevelSolver()->ReInit("seaice");
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(h);
  GetMultiLevelSolver()->ReInitVector(div);
  GetMultiLevelSolver()->ReInitVector(oldu);
  GetMultiLevelSolver()->ReInitVector(extu);
  GetMultiLevelSolver()->ReInitVector(tmp);
  GetMultiLevelSolver()->ReInitVector(f);
    
  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();


  
  cout << "------------------------------" << endl;
  cout << "See-Eis" << endl;
  cout << "------------------------------" << endl << endl;
  
  // See-Eis
  GetMultiLevelSolver()->SetProblem("seaice");
  InitSolution(u);
  InitSolution(oldu);
  
  // Transport
  GetMultiLevelSolver()->SetProblem("tg");
  InitSolution(h);

  
  GetMultiLevelSolver()->SetProblem("seaice");
  // Masse-Matrix
  AssembleMassMatrix();

  

  // set up system matrix
  double writeeveryhour=0.0; // 1.0;

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
  
  for (_iter=1; _iter<=_niter; _iter++)
    {

       if ((TIME>20)&&(timeinc==0)) { DT*=2.0; ++timeinc; }
      if ((TIME>40)&&(timeinc==1)) { DT*=2.0; ++timeinc; }
       //  if ((TIME>60)&&(timeinc==2)) { DT*=2.0; ++timeinc; }
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

      int NSUB = 20;
      for (int sub=0;sub<NSUB;++sub)
	{
	  DTSUB = DT/NSUB;
	  GetMultiLevelSolver()->GetSolver()->Equ(u,( (double) (NSUB-sub))/ ( (double) NSUB),oldu); // besserer Startwert fuer newton
	  GetMultiLevelSolver()->GetSolver()->Add(u,( (double) (sub))/ ( (double) NSUB),extu);      // besserer Startwert fuer newton
	 
	  SolveDIV(div,u,f);
	  SolveTransport(h,div, u, tmp, f);
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
      GetMultiLevelSolver()->AddNodeVector("extu", extu);
      GetMultiLevelSolver()->AddNodeVector("H", h);
      
      start = clock();
      //    cout <<"Solve Seaice"<< endl;
       res = Solve(u,f,"Results/u");
         assert(res == "converged");
      
      end = clock();
      
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
      //   cout<< " cpu_time_used"<<  cpu_time_used << endl;
      //cout<< " start"<<  start << "end"<< end << endl;
     
        functionals = Functionals(u,f);
     
      GetMultiLevelSolver()->DeleteNodeVector("oldu");
      GetMultiLevelSolver()->DeleteNodeVector("extu");
      GetMultiLevelSolver()->DeleteNodeVector("H");
      
      
      
        OUTF << TIME << " " << functionals <<  cpu_time_used <<endl;
      
      
      double time_in_hours = TIME*tref/60/60;
      if (time_in_hours>writenext)
	{
	  
	  writenext += writeeveryhour;
	  ++writeiter;
	  
	   GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,writeiter);
	  GetMultiLevelSolver()->GetSolver()->Visu("Results/h",h,writeiter);
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








