/*----------------------------   loop.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __loop_H
#define __loop_H
/*----------------------------   loop.h     ---------------------------*/

#include "stdloop.h"
#include "meshagent.h"
#include "boundaryfunction.h"
#include <fstream>
#include "matrix.h"

#include "gascoigne.h"
#include "gascoignemesh2d.h"

namespace Gascoigne
{
  
extern Timer GlobalTimer;
}


/// @brief Describes circle for curved boundary
class Disc : public Gascoigne::BoundaryFunction<2>
{
  double      _r2;
  Gascoigne::Vertex<2> _midpoint;

public:
  /// Initialize midpoint and square of radius
  void BasicInit(Gascoigne::Vertex<2> m, double r2)
  {
    _midpoint = m; _r2 = r2;
  }

  /// The boundary is the zero-contour of the given function
  double operator()(const Gascoigne::Vertex<2> &x) const
  {
    Gascoigne::Vertex<2> c = x;
    c.add(-1.,_midpoint);
    return _r2 - c*c;
  }
};

/// @brief Describes cylinder for curved boundary in 3d
class CylinderZ : public Gascoigne::BoundaryFunction<3>
{
  double      _r2;
  Gascoigne::Vertex<3> _midpoint;

public:
  /// Initialize midpoint and square of radius
  void BasicInit(Gascoigne::Vertex<3> m, double r2)
  {
    _midpoint = m; _r2 = r2;
  }

  /// The boundary is the zero-contour of the given function
  double operator()(const Gascoigne::Vertex<3> &x) const
  {
    Gascoigne::Vertex<3> c = x;
    c.add(-1.,_midpoint);
    return _r2 - c.x()*c.x() - c.y()*c.y();
  }
};

/// \brief Mesh agent that initialized the curved shape
class BenchmarkMeshAgent : public Gascoigne::MeshAgent
{
protected:
  Disc      disc;
  CylinderZ cyl;
  
public:
  BenchmarkMeshAgent(std::string config) : Gascoigne::MeshAgent()
  {
    if ( (config == "bench2d1") || (config == "bench2d2") )
      {
	// initialize disc of radius 0.05 at (0.2,0.2)
	disc.BasicInit(Gascoigne::Vertex2d(0.2,0.2),0.05*0.05);
	AddShape(80,&disc);
      }
    else if ( (config == "bench3d1") || (config == "bench3d2") || (config == "bench3d3") )
      {
	// initialize disc of radius 0.05 at (0.2,0.2)
	cyl.BasicInit(Gascoigne::Vertex3d(0.5,0.2,0.0),0.05*0.05);
	AddShape(80,&cyl);
      }
  }
};



/// \brief Derived loop to setup the curved geometry and to control algorithms
class Loop : public Gascoigne::StdLoop
{
public:

  // initialize our own mesh agent
  void BasicInit(std::string config,
		 const Gascoigne::ParamFile& paramfile,
		 const Gascoigne::ProblemContainer* PC,
		 const Gascoigne::FunctionalContainer* FC) 
    {
      GetMeshAgentPointer() = new BenchmarkMeshAgent(config);
      StdLoop::BasicInit(paramfile, PC, FC);
      
    }

  
  /**
   *  Loop for a stationary problem on a sequence of uniformly refined meshes
   */
  void stationaryrun(std::string label)
  {
    Gascoigne::Vector u("u"), f("f");
    Gascoigne::GlobalVector ucoarse;
    
    Gascoigne::Matrix A("A");

    Gascoigne::GlobalTimer.reset();
    
    for (_iter = 1; _iter <= _niter; _iter++)
      {
	Gascoigne::GlobalTimer.start("iteration");
	std::cout << "\n================== " << _iter << " ================";
	PrintMeshInformation();

	Gascoigne::GlobalTimer.start("init");
	GetMultiLevelSolver()->ReInit();
	GetMultiLevelSolver()->SetProblem(label);
	GetMultiLevelSolver()->ReInitMatrix(A);
	GetMultiLevelSolver()->ReInitVector(u);
	GetMultiLevelSolver()->ReInitVector(f);

	// Init the solution or interpolate it from the previous mesh level
	if (_iter>1)
	  GetMultiLevelSolver()->InterpolateSolution(u, ucoarse);
	else
	  {
	    GetMultiLevelSolver()->GetSolver()->OutputSettings();
	    InitSolution(u);
	  }
	GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
	Gascoigne::GlobalTimer.stop("init");

	Gascoigne::GlobalTimer.start("solve");
	Solve(A,u, f);
	Gascoigne::GlobalTimer.stop("solve");

	Gascoigne::GlobalTimer.start("functionals");
	Gascoigne::DoubleVector juh = Functionals(u, f);
	Gascoigne::GlobalTimer.stop("functionals");

	
	Gascoigne::GlobalTimer.start("refine");
	if (_iter < _niter)
	  {
	    CopyVector(ucoarse, u);
	    Gascoigne::DoubleVector eta;
	    AdaptMesh(eta);
	  }
	Gascoigne::GlobalTimer.stop("refine");
	Gascoigne::GlobalTimer.stop("iteration");

	Gascoigne::GlobalTimer.print();
	Gascoigne::GlobalTimer.reset();
      }

    if (_runtime_statistics)
      {
	std::cout << std::endl;
	std::cout.precision(4);
	Gascoigne::GlobalTimer.print100();
      }
  }
  
  

  /**
   *  Loop for a nonstationary problem on a fixed mesh with fixed step size
   */
  void timerun(std::string label)
  {
    /// stores solution, right hand side and solution of previous step
    Gascoigne::Vector u("u"), f("f"), old("old");
    /// stores the system matrix
    Gascoigne::Matrix A("A");
    
    
    double time, dt, stoptime;
    Gascoigne::DataFormatHandler DFH;
    DFH.insert("starttime", &time,     0.);
    DFH.insert("stoptime",  &stoptime, 0.);
    DFH.insert("dt",  &dt, 0.);
    Gascoigne::FileScanner FS(DFH, _paramfile, "Equation");
    assert(dt>0);
      

    _niter = (stoptime-time+1.e-10)/dt;
    if (fabs(stoptime - time - _niter*dt)>1.e-8)
      {
	std::cerr << "The length of the time interval "
		  << time << " to " << stoptime
		  << " is no multiple of the step size " << dt << std::endl;
      }
    
    /// Initialization of the problem
    PrintMeshInformation();
    GetMultiLevelSolver()->ReInit();
    GetMultiLevelSolver()->SetProblem(label);
    GetMultiLevelSolver()->ReInitMatrix(A);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(old);
    GetMultiLevelSolver()->ReInitVector(f);
    InitSolution(u);

    /// Matrix must always be build before the first step
    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;

    std::ofstream LOG("functionals.txt");
    
    
    for (_iter = 1; _iter <= _niter; _iter++)
      {
	std::cout << std::endl << "----------------------------------" << std::endl;
	std::cout << "time step " << _iter << " \t "
		  << time << " -> " << time + dt << "\t[" << dt << "]" << std::endl;
	time += dt;

	// save old solution
	GetMultiLevelSolver()->Equ(old,1.,u);

	// and pass it to the discretization
	GetMultiLevelSolver()->AddNodeVector("OLD",old);

	// Solve the time step
	Solve(A, u, f);

	// compute Functionals and write status file
	Gascoigne::DoubleVector juh = Functionals(u, f);
	LOG << time << "\t" << juh << std::endl;


	
	// delete pointer to old solution
	GetMultiLevelSolver()->DeleteNodeVector("OLD");
      }
    LOG.close();
  }









  /**
   *  Loop for a nonstationary problem on a fixed mesh with fixed step size 
   *  using BDF2 and explicit convection handling
   */
  void timerunbdf2(std::string label)
  {
    /// stores solution, right hand side and solution of previous step
    Gascoigne::Vector u("u"), f("f"), old("old"), oldold("oldold");
    /// stores the system matrix
    Gascoigne::Matrix A("A");
    
    
    double time, dt, stoptime;
    Gascoigne::DataFormatHandler DFH;
    DFH.insert("starttime", &time,     0.);
    DFH.insert("stoptime",  &stoptime, 0.);
    DFH.insert("dt",  &dt, 0.);
    Gascoigne::FileScanner FS(DFH, _paramfile, "Equation");
    assert(dt>0);
    
    
    _niter = (stoptime-time+1.e-10)/dt;
    if (fabs(stoptime - time - _niter*dt)>1.e-8)
      {
	std::cerr << "The length of the time interval "
		  << time << " to " << stoptime
		  << " is no multiple of the step size " << dt << std::endl;
      }
    
    /// Initialization of the problem
    PrintMeshInformation();
    GetMultiLevelSolver()->ReInit();
    GetMultiLevelSolver()->SetProblem(label);
    GetMultiLevelSolver()->ReInitMatrix(A);
    GetMultiLevelSolver()->ReInitVector(u);
    GetMultiLevelSolver()->ReInitVector(old);
    GetMultiLevelSolver()->ReInitVector(oldold);
    GetMultiLevelSolver()->ReInitVector(f);
    InitSolution(u);
    GetMultiLevelSolver()->Equ(oldold, 1., u);
    GetMultiLevelSolver()->Equ(old,    1., u);
    
    /// Matrix must always be build before the first step
    GetSolverInfos()->GetNLInfo().control().matrixmustbebuild() = 1;
    
    std::ofstream LOG("functionals.txt");
    
    
    for (_iter = 1; _iter <= _niter; _iter++)
      {
	std::cout << std::endl << "----------------------------------" << std::endl;
	std::cout << "time step " << _iter << " \t "
		  << time << " -> " << time + dt << "\t[" << dt << "]" << std::endl;
	time += dt;

	// save old solution
	GetMultiLevelSolver()->Equ(oldold,1.,old);
	GetMultiLevelSolver()->Equ(old,1.,u);

	// and pass it to the discretization
	GetMultiLevelSolver()->AddNodeVector("OLD",old);
	GetMultiLevelSolver()->AddNodeVector("OLDOLD",oldold);

	// Solve the time step
	Solve(A, u, f);

	// compute Functionals and write status file
	Gascoigne::DoubleVector juh = Functionals(u, f);
	LOG << time << "\t" << juh << std::endl;

	
	// delete pointer to old solution
	GetMultiLevelSolver()->DeleteNodeVector("OLD");
	GetMultiLevelSolver()->DeleteNodeVector("OLDOLD");
      }
    LOG.close();
  }
};


/*----------------------------   loop.h     ---------------------------*/
/* end of #ifndef __loop_H */
#endif
/*----------------------------   loop.h     ---------------------------*/
