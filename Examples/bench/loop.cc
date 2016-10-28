#include "loop.h"
#include "gascoignemesh2d.h"
#include <time.h>
#include  "backup.h"
#include  "stdmultilevelsolver.h"
#include  "sparseblockmatrix.h"
#include  "fmatrixblock.h"
#include  "compose_name.h"

using namespace Gascoigne;
using namespace std;
 


extern ofstream ELLIPSE_OUT;

double TIME, DT, DTSUB, DELTAMIN,DTUL;
extern double STEUERUNG_MU;

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
			  VectorInterface& hl,
			  VectorInterface& hh,
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
  GlobalVector& UL  = GetMultiLevelSolver()->GetSolver()->GetGV(hl);
  GlobalVector& DIV  = GetMultiLevelSolver()->GetSolver()->GetGV(div);
  
  for (int c=0;c<H.ncomp();++c)
    for (int i=0;i<H.n();++i)
      //  H(i,c) = F(i,c) / ( LMM[i] * (1.0 + DTSUB * std::max(0.0,DIV(i,0)) ) );
      H(i,c) = F(i,c) / ( LMM[i] );
  
  UL=H;


  // High-Order Loesung
  GlobalVector& UH  = GetMultiLevelSolver()->GetSolver()->GetGV(hh);
  UH.zero();
  TimePattern TP(2,2); TP.zero(); TP(0,0)=1.0; TP(1,1) = 1.0;
  PrecondCGMass(UH, F, TP, 1.0);
  
 
   
  // fluxidux
  const ColumnDiagStencil* ST = dynamic_cast<const ColumnDiagStencil* >(MM.GetStencil());
  assert(ST);
  
  nvector<double> maxi(2),mini(2),meanmini(2), meanmaxi(2), Pminus(2), Pplus(2), Qminus(2), Qplus(2);
  GlobalVector Rplus(H.ncomp(),H.n()), Rminus(H.ncomp(),H.n());
  for (int row = 0; row<ST->n(); ++row)
    {
      // max und min von low-order , Pppus, Pminus
      mini[0]=maxi[0] = H(row,0);
      mini[1]=maxi[1] = H(row,1);

      meanmini[0]=meanmaxi[0] = 0.0;
      meanmini[1]=meanmaxi[1] = 0.0;
      Pminus=-1.e-10; Pplus = 1.e-10;
      for (int pos = ST->start(row); pos<ST->stop(row); ++pos)
	{
	  int col = ST->col(pos);
	  mini[0] = min(mini[0],H(col,0));  mini[1] = min(mini[1],H(col,1));
	  maxi[0] = max(maxi[0],H(col,0));  maxi[1] = max(maxi[1],H(col,1));

	  meanmini[0]+= min(0,H(col,0)); meanmini[1]+= min(0,H(col,1));
	  meanmaxi[0]+= max(0,H(col,0)); meanmaxi[1]+= max(0,H(col,1));

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
	  // Prelimiting
	  
	  if(fij0*(UL(col,0)-UL(row,0))>0)
	    {
	      fij0=0;
	    }
	  if(fij1*(UL(col,1)-UL(row,1))>0)
	    { 
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
}


void Loop::run(const std::string& problemlabel)
{
  TIME=0.0;
  double tref;
  double dtmax, endtime;
  int prerefine;
  string _reloadu,_reloadh,_reloadoldu;
  if (1)
    {
      DataFormatHandler DFH;
      DFH.insert("dt",   &DT, 0.);
      DFH.insert("dtmax",   &dtmax, 0.);
      DFH.insert("time", &TIME, 0.);
      DFH.insert("endtime", &endtime, 0.);
      DFH.insert("reloadu", &_reloadu);
      DFH.insert("reloadh", &_reloadh);
      DFH.insert("reloadoldu", &_reloadoldu);
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
  VectorInterface u("u"), h("h"), f("f"), oldu("oldu"),hl("hl"), hh("hh"), extu("extu"), div("div"),oldh("oldh"),oldoldu("oldoldu"),other("other"); 



  PrintMeshInformation();
  // initialize problem, solver and vectors 
  GetMultiLevelSolver()->ReInit("seaice");
  GetMultiLevelSolver()->ReInitVector(u);
  GetMultiLevelSolver()->ReInitVector(h);
  GetMultiLevelSolver()->ReInitVector(hl);
  GetMultiLevelSolver()->ReInitVector(div);
  GetMultiLevelSolver()->ReInitVector(oldu);
  GetMultiLevelSolver()->ReInitVector(extu);
  GetMultiLevelSolver()->ReInitVector(hh);
  GetMultiLevelSolver()->ReInitVector(f);
  GetMultiLevelSolver()->ReInitVector(oldoldu); 
  GetMultiLevelSolver()->ReInitVector(oldh);
  GetMultiLevelSolver()->ReInitVector(other);
    
  GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
  GetMultiLevelSolver()->GetSolver()->OutputSettings();




  cout << "------------------------------" << endl;
  cout << "See-Eis" << endl;
  cout << "------------------------------" << endl << endl;
  
  // See-Eis
  GetMultiLevelSolver()->SetProblem("seaice");
  if (_initial=="reload")	
    {
      GetMultiLevelSolver()->GetSolver()->Read(u,_reloadu);
      GetMultiLevelSolver()->GetSolver()->Read(oldu,_reloadoldu);    
    }
  else 
    {
      InitSolution(u);     
      GetMultiLevelSolver()->Equ(oldu,1.0,u);
    }
  // Transport
  GetMultiLevelSolver()->SetProblem("tg");
  if (_initial=="reload")
    { GetMultiLevelSolver()->GetSolver()->Read(h,_reloadh);
      }
  else 
    InitSolution(h);

  GetMultiLevelSolver()->SetProblem("seaice");
  // Masse-Matrix
  AssembleMassMatrix();

  // set up system matrix
  double writeeveryhour=0.0; // 1.0;
  double stepback=0.0;
  double writenext= writeeveryhour;
  int writeiter = 0;
  string res;
  
  nvector<double> functionals;
  stringstream str;
  str << "func_" << prerefine << "_" << dtmax << ".txt";
  ofstream OUTF(str.str().c_str());

  int timeinc=0;
  
  clock_t start, end;
  double cpu_time_used,a;

  

  GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,writeiter);
  GetMultiLevelSolver()->GetSolver()->Visu("Results/h",h,writeiter);
  GetMultiLevelSolver()->GetSolver()->Visu("Results/o",other,writeiter);


  
  for (_iter=1; _iter<=_niter; _iter++)
    {

      /*     
	     if(timeinc>0)
	     {   timeinc++;

	     if (timeinc==33)
	     {DT=DT*2.0;
	     }

	     if (timeinc==53)
	     {DT=DT*2.0;

	     timeinc=0.0;
	     cout<<"Doppeltzeit"<<endl;
	     }

	     }

	     if(stepback==0)
	     TIME += DT;

	     if(stepback>0.0)
	     {
	     cout<<"REDUCE"<<endl;
	     TIME-=DT;
	     DT=DT/4;
	     stepback=0;
	     TIME += DT;
	     } 
      */

      TIME += DT;


      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

      cout << "\n Zeitschritt " << _iter << " " << TIME-DT << " -> " 
	   << TIME<< "  [" << DT << "]\t"  
	   << "\t" << TIME*tref/60/60/24 << " Tage" << endl;

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


      // Transoprtgleichung 
      /////////////////////// 
      /////////////////////// 
      /////////////////////// 
      /////////////////////// 
      /////////////////////// 
      GetMultiLevelSolver()->GetSolver()->Equ(oldh,1.0, h);
      int NSUB = 20;
      for (int sub=0;sub<NSUB;++sub)
	{
	  DTSUB = DT/NSUB;
	  GetMultiLevelSolver()->GetSolver()->Equ(u,( (double) (NSUB-sub)) / ( (double) NSUB),oldu); // besserer Startwert fuer newton
	  GetMultiLevelSolver()->GetSolver()->Add(u,( (double) (sub))      / ( (double) NSUB),extu); // besserer Startwert fuer newton
	 
	  // Solve(div,u,f);
	  SolveTransport(h,div, u,hl, hh, f);

	  GlobalVector& H = GetMultiLevelSolver()->GetSolver()->GetGV(h);
	  for (int i=0;i<H.n();++i)
	    {
	      H(i,1) = std::min(1.0, H(i,1));
	      H(i,1) = std::max(0.0, H(i,1));
	      H(i,0) = std::max(0.0, H(i,0));
	    }
	  GlobalVector& UL = GetMultiLevelSolver()->GetSolver()->GetGV(hl);
	  for (int i=0;i<UL.n();++i)
	    {
	      UL(i,1) = std::min(1.0, UL(i,1));
	      UL(i,1) = std::max(0.0, UL(i,1));
	      UL(i,0) = std::max(0.0, UL(i,0));
	    }
	}
        
 
      // Eis-Problem
 
      GetMultiLevelSolver()->SetProblem("seaice");
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      GetMultiLevelSolver()->AddNodeVector("oldu", oldu);
      GetMultiLevelSolver()->AddNodeVector("extu", extu);
      GetMultiLevelSolver()->AddNodeVector("H", hl);     ////// Die Low-Order-Loesung geht in die Momentengleichung ein.

      GetMultiLevelSolver()->GetSolver()->Equ(oldu,1.0, u);
      GetMultiLevelSolver()->GetSolver()->Equ(oldoldu,1.0, oldu);
   
      GetMultiLevelSolver()->AssembleMatrix(u);
      GetMultiLevelSolver()->ComputeIlu(u);  
      res = Solve(u,f,"Results/u");



      if (res!="converged")
	{
	  abort();
	  a=1.0;

	  /*
	    timeinc++;
	    stepback++;

 
	    GetMultiLevelSolver()->GetSolver()->Equ(u,1.0,oldu);
	    GetMultiLevelSolver()->GetSolver()->Equ(oldu,1.0,oldoldu);
     
	    GetMultiLevelSolver()->GetSolver()->Equ(hl,1.0, oldh);
   
    

	    GetMultiLevelSolver()->DeleteNodeVector("oldu");
	    GetMultiLevelSolver()->DeleteNodeVector("extu");
	    GetMultiLevelSolver()->DeleteNodeVector("H");

      
	    continue;
	  */
	}
      else
	{
	  stepback=0;
	  a=0.0;
	 
	}
  

      
      string ell_name="ellipse";
      compose_name(ell_name,_iter);
      ELLIPSE_OUT.open(ell_name.c_str());
      functionals = Functionals(u,f);
      ELLIPSE_OUT.close();
     
      GetMultiLevelSolver()->DeleteNodeVector("oldu");
      GetMultiLevelSolver()->DeleteNodeVector("extu");
      GetMultiLevelSolver()->DeleteNodeVector("H");
      
      // Other-Problem
      GetMultiLevelSolver()->SetProblem("other");
      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      GetMultiLevelSolver()->AddNodeVector("U", u);
      GetMultiLevelSolver()->AddNodeVector("H", h);
      Solve(other,f,"Results/o");
      GetMultiLevelSolver()->DeleteNodeVector("U");
      GetMultiLevelSolver()->DeleteNodeVector("H");
      
      OUTF << TIME << " " << functionals  <<a<<endl;
      //  jetzt ist gerechnet und es gibt umittel, hmittel


      double time_in_hours = TIME*tref/60/60;
      if (time_in_hours>writenext)
	{
	  
	  writenext += writeeveryhour;
	  ++writeiter;
        
	  GetMultiLevelSolver()->GetSolver()->Visu("Results/u",u,writeiter);
	  
 
	  GetMultiLevelSolver()->GetSolver()->Visu("Results/h",h,writeiter);
	  
	  GetMultiLevelSolver()->GetSolver()->Visu("Results/hh",hh,writeiter);
	  
	  GetMultiLevelSolver()->GetSolver()->Visu("Results/hl",hl,writeiter);
    	  GetMultiLevelSolver()->GetSolver()->Visu("Results/o",other,writeiter);

	  string name="Results/h";
	  compose_name(name,writeiter);
	  GetMultiLevelSolver()->GetSolver()->Write(h,name);

	  cout << "[" << name << ".bup]";
	}
     
      
      if (TIME >= endtime*24*60*60/tref) break;
      
    }
  OUTF.close();
}










