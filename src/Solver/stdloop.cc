#include  "stdloop.h"
#include  "adaptordata.h"
#include  "diplomantenadaptor.h"
#include  "malteadaptor.h"
#include  "monitoring.h"
#include  "filescanner.h"
#include  "stdmultilevelsolver.h"
#include  "meshagent.h"
#include  "stdsolver.h"
#include  "energyestimator.h"


using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
StdLoop::StdLoop() : BasicLoop()//, _paramfile(NULL)
{
  _estimator = _extrapolate = "none";
  _FV.resize(0);
}

/*-----------------------------------------*/

StdLoop::~StdLoop()
{
}

/*-----------------------------------------*/

void StdLoop::ClockOutput() const
{
  BasicLoop::ClockOutput();
  cout << "  Functionals\t\t" << _clock_functionals.read() << endl;
  cout << "  Estimate\t\t" << _clock_estimate.read() << endl;
}

/*-----------------------------------------*/

void StdLoop::BasicInit(const ParamFile* paramfile)
{
  BasicLoop::BasicInit(paramfile);

  DataFormatHandler DFH;
  DFH.insert("nmin",               &_nmin,              1000);
  DFH.insert("nmax",               &_nmax,              100000);
  DFH.insert("p",                  &_p,                 0.1);
  DFH.insert("random_coarsening",  &_random_coarsening, 0);
  DFH.insert("coarse",             &_coarse,            0);
  DFH.insert("refiner",            &_refiner,           "global");
  DFH.insert("estimator",          &_estimator,         "none");
  DFH.insert("extrapolate",        &_extrapolate,       "no");
  DFH.insert("resultsdir",         &_s_resultsdir,      "Results");
  FileScanner FS(DFH);
  FS.NoComplain();
  FS.readfile(_paramfile,"Loop");
}

/*-------------------------------------------------------*/

DoubleVector StdLoop::ComputeFunctionals(MultiLevelGhostVector& f, MultiLevelGhostVector& u, const vector<const Functional*>& J) const
{
  int n = J.size(); 
  DoubleVector j(n,0.);
  if (n==0) return j;

  cout << "\nFunctionals: ";
  for(int i=0; i<n; i++)
    {
      j[i] = GetMultiLevelSolver()->ComputeFunctional(f,u,J[i]);
      cout << J[i]->GetName() << " ";
    }
  cout << endl;
  return j;
} 

/*-----------------------------------------*/

void StdLoop::EtaVisu(string name, int i, const DoubleVector& eta) const
{
  Visualization Visu;
  Visu.format("vtk");
  Visu.set_name(name);
  Visu.step(i);
  VisuDataInfo     VDI(1);
  VisuDataNVector  VD(eta);
  
  Visu.SetPointData(&VD);
  Visu.SetMesh(GetMeshAgent()->GetMesh(0));
  Visu.SetPointDataInfo(&VDI);
  Visu.write();
}

/*-------------------------------------------------*/

DoubleVector StdLoop::GetExactValues() const
{
  int n = _FV.size(); 
  DoubleVector j(n,0.);
  for(int i=0; i<n; i++)
    {
      j[i] = _FV[i]->ExactValue();
    }
  return j; 
}

/*-------------------------------------------------*/

DoubleVector StdLoop::Functionals(MultiLevelGhostVector& u, MultiLevelGhostVector& f)
{
  DoubleVector J  = ComputeFunctionals(f,u,_FV);
  _JErr.resize(J.size());
  if (J.size())
    {
      cout << "\nvalue";
      cout.precision(16);
      for(int i=0; i<J.size(); i++) 
        {
          cout << "\t" << J[i];
        }
      cout << "\nerror";
      for(int i=0; i<J.size(); i++) 
        {
          _JErr[i] = GetExactValues()[i] - J[i];
          cout << "\t" << _JErr[i] ;
        }
      cout << endl;
      if(_extrapolate=="yes")
        {
          Extra.NewValues(J);
          Extra.Print();
        }
      cout << endl;
    }
  return J;
}

/*-------------------------------------------------*/

double StdLoop::Estimator(DoubleVector& eta, MultiLevelGhostVector& u, MultiLevelGhostVector& f)
{
  double est = 0.;
  if (_estimator=="energy")
    {
      StdSolver* S = dynamic_cast<StdSolver*>(GetMultiLevelSolver()->GetSolver());
      assert(S);
      
      S->GetHierarchicalMeshPointer() = dynamic_cast<MeshAgent*>(GetMeshAgent())->GetHierarchicalMesh();

      EnergyEstimator E(*S);
      est = E.Estimator(eta,u,f);
      //est = GetMultiLevelSolver()->GetSolver()->EnergyEstimator(eta, u, f);
      EtaVisu(_s_resultsdir+"/eta",_iter,eta);
    }
  else 
    {
      cout << "Estimator\t"; 
      cout << "Baustelle!!!!!\n";
    }
//   double eff = 0.;
//   if      ((_estimator=="weighted") && (_JErr.size()>0))      eff = est/_JErr[0];
//   else if ((_estimator=="energy")   && (_GlobalErr.size()>0)) eff = est/_GlobalErr(1,0);
//   if  (eff!=0.) cout << " @ " << eff << endl; 
  
//   cout << endl; 
  return est;
}

/*-------------------------------------------------*/

void StdLoop::AdaptMesh(const DoubleVector& eta,std::string refine_or_coarsen_step)
{
  //das gleichzeitige vergroebern und verfeinern FUNKTIONIERT nicht
  //wer das machen moechte, muss stattdessen und zwei getrennten laeufen 
  //das gitter vergroebern, reinit+interpolate und dann das gitter verfeinern
  if(refine_or_coarsen_step=="refine") ;
  else if(refine_or_coarsen_step=="coarsen") ;
  else {
    std::cerr<<"the variable refine_or_coarsen_step has to be set, either to 'refine' or 'coarsen'"<<std::endl;
    assert(0);
  }

  if     (_refiner=="global") {
    if(refine_or_coarsen_step=="refine"){
      GetMeshAgent()->global_refine(1);
    }
  }  
  else if(_refiner=="none") GetMeshAgent()->global_refine(0);
  else if(_refiner=="random") 
    {
      if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
      if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;

      if(refine_or_coarsen_step=="refine"){
        GetMeshAgent()->random_patch_refine(_p,0);
      }
      if(refine_or_coarsen_step=="coarsen"){
        if(_random_coarsening){
          GetMeshAgent()->random_patch_coarsen(_p,0);
        }
      }
    }
  else if(_refiner=="random_refine") 
    {
      if(refine_or_coarsen_step=="refine"){
        if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
        if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;
        GetMeshAgent()->random_patch_refine(_p,0);
      }
    }
  else if(_refiner=="random_coarsen") 
    {
      if(refine_or_coarsen_step=="coarsen"){
        if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
        if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;
        GetMeshAgent()->random_patch_coarsen(_p,0);
      }
    }
  else if(_refiner=="eta") 
    {
      IntVector refnodes, coarsenodes,dummynodes;

      MalteAdaptor A(_paramfile,eta);
      A.refine(refnodes,coarsenodes);

      if(refine_or_coarsen_step=="coarsen") GetMeshAgent()->refine_nodes(dummynodes,coarsenodes);
      if(refine_or_coarsen_step=="refine")  GetMeshAgent()->refine_nodes(refnodes,dummynodes);
    }
  else if(_refiner=="dip") 
    {
      IntVector refnodes, coarsenodes;

      AdaptorData info;
      info.rfactor() = 1.; 
      DiplomandenAdaptor A(info,eta);
      A.refine(refnodes);
      if(refine_or_coarsen_step=="refine")  GetMeshAgent()->refine_nodes(refnodes,coarsenodes);
    }
  else assert(0);
}

/*-------------------------------------------------*/

void StdLoop::AdaptMesh(const DoubleVector& eta)
{
  //das gleichzeitige vergroebern und verfeinern FUNKTIONIERT nicht
  //wer das machen moechte, sollte stattdessen zwei getrennte laeufe durchfuehren:
  //das gitter vergroebern, reinit+interpolate und dann das gitter verfeinern
  //das entsprechend die methode AdaptMesh(eta,refine_or_coarsen_step) aufrufen
  if     (_refiner=="global") GetMeshAgent()->global_refine(1);
  else if(_refiner=="none")   GetMeshAgent()->global_refine(0);
  else if(_refiner=="random") 
    {
      if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
      if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;
      if(_random_coarsening){
        std::cerr<<"das gleichzeitige vergroebern und verfeinern FUNKTIONIERT nicht"<<std::endl;
        std::cerr<<"fuehren Sie stattdessen zwei getrennte laeufe durch: random_refine, random_coarsen"<<std::endl;
        std::cerr<<"rufen Sie dazu AdaptMesh(eta,refine_or_coarsen_step) auf"<<std::endl;
        assert(0);
      }
      GetMeshAgent()->random_patch_refine(_p,_random_coarsening);
    }
  else if(_refiner=="random_refine") 
    {
      if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
      if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;
      GetMeshAgent()->random_patch_refine(_p,0);
    }
  else if(_refiner=="random_coarsen") 
    {
      if (GetMeshAgent()->nnodes()>_nmax) _p *= 0.5;
      if (GetMeshAgent()->nnodes()<_nmin) _p *= 1.1;
      GetMeshAgent()->random_patch_coarsen(_p,0);
    }
  else if(_refiner=="eta") 
    {
      IntVector refnodes, coarsenodes;

      MalteAdaptor A(_paramfile,eta);
      A.refine(refnodes,coarsenodes);

      if(refnodes.size()>0 && coarsenodes.size()>0){
        std::cerr<<"das gleichzeitige vergroebern und verfeinern FUNKTIONIERT nicht"<<std::endl;
        std::cerr<<"fuehren Sie stattdessen zwei getrennte laeufe durch, einmal vergroebern, einmal verfeinern"<<std::endl;
        std::cerr<<"rufen Sie dazu AdaptMesh(eta,refine_or_coarsen_step) auf"<<std::endl;
        assert(0);
      }

      GetMeshAgent()->refine_nodes(refnodes,coarsenodes);
    }
  else if(_refiner=="dip") 
    {
      IntVector refnodes, coarsenodes;

      AdaptorData info;
      info.rfactor() = 1.; 
      DiplomandenAdaptor A(info,eta);
      A.refine(refnodes);
      GetMeshAgent()->refine_nodes(refnodes,coarsenodes);
    }
  else assert(0);
}

/*-------------------------------------------------*/

void StdLoop::run(const ProblemDescriptorInterface* PD)
{
  MultiLevelGhostVector u("u"), f("f");
  u.SetMultiLevelSolver(GetMultiLevelSolver());
  f.SetMultiLevelSolver(GetMultiLevelSolver());
  GlobalVector  ualt;

  GetMultiLevelSolver()->RegisterVector(u);
  GetMultiLevelSolver()->RegisterVector(f);
  
  Monitoring Moning;
  
  for (_iter=1; _iter<=_niter; _iter++)
    {
      cout << "\n================== " << _iter << " ================";
      PrintMeshInformation();
      Moning.SetMeshInformation(_iter,GetMeshAgent()->nnodes(),GetMeshAgent()->ncells());
      
      _clock_newmesh.start();

      GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
      GetMultiLevelSolver()->ReInit(*PD);
      GetMultiLevelSolver()->InterpolateSolution(u,ualt);
      GetMultiLevelSolver()->GetSolver()->Visu(_s_resultsdir+"/interpolate",u,_iter);

      _clock_newmesh.stop();

      if (_iter==1) 
        {
          GetMultiLevelSolver()->GetSolver()->OutputSettings();
          InitSolution(u);
          Moning.BasicInit(GetExactValues());
        }

      Solve(u,f);
      ComputeGlobalErrors(u);
      
      _clock_functionals.start();
      DoubleVector juh = Functionals(u,f);
      _clock_functionals.stop();

      DoubleVector eta;

      _clock_estimate.start();
      if (_estimator!="none")
        {
          double est = Estimator(eta,u,f);
          Moning.SetSolutionInformation(_iter,juh,est);
        }
      if (_iter<_niter) 
        {
          CopyVector(ualt,u);

          AdaptMesh(eta);
          //  wenn gleichzeitig verfeinert und vergroebert werden soll, dann rufen Sie
          //  AdaptMesh zweimal folgendermassen auf:
          //     
          //    CopyVector(ualt,u);
          //    AdaptMesh(eta,"refine");
          //    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
          //    GetMultiLevelSolver()->ReInit(*PD);
          //    GetMultiLevelSolver()->InterpolateSolution(u,ualt);
          //    CopyVector(ualt,u);
          //    AdaptMesh(eta,"coarsen");
          //   
        }
      _clock_estimate.stop();
     }
}
}

/*-------------------------------------------------*/
