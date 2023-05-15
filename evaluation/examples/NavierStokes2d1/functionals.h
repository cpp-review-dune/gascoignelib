/*----------------------------   functionals.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __functionals_H
#define __functionals_H
/*----------------------------   functionals.h     ---------------------------*/


#include "domainfunctional.h"
#include "dirichletdatabycolor.h"
#include "domainfunctional.h"
#include "residualfunctional.h"
#include "weightedpointfunctional.h"
#include "boundaryfunctional.h"

#include "navierstokes.h"


// Defines the drag as ResidualFunctional, i.e. by using the 
// Babuska-Miller trick for evaluation
class Drag : public virtual Gascoigne::ResidualFunctional
{
public:
  Drag(std::string config) : ResidualFunctional()
  {
    __comps.push_back(1);
    __cols.insert(80);

    if (config=="bench2d1")
      {
	__scales.push_back(500.0);
	// exact value from extrapolation on very fine meshes
	ExactValue() = 5.579532446386295;
      }
    else if (config=="bench2d2")
      {
	__scales.push_back(20.0);
      }
    else if (config=="bench3d1")
      {
	__scales.push_back(500.0/0.41);
	// Taken from Braack & Richter: Solutions of 3D Navier-Stokes benchmarkproblems with adaptive finite elements
	ExactValue() = 5.579532446386295;
      }
    else if ( (config=="bench3d2") || (config=="bench3d3") )
      {
	__scales.push_back(20.0/0.41);
	// no exact values. The original paper by Schaefer/Turek did not consider the nonstationary 3d case with circular obstacle.
	ExactValue() = 0.0;
      }
    else
      {
	std::cerr << "Configuration " << config << " not known!" << std::endl;
	abort();
      }
    
    __DD = new Gascoigne::DirichletDataByColor(GetComps(), GetColors(), GetScales());
  }

  std::string GetName() const
  {
    return "Drag";
  }
};
// Defines the lift as ResidualFunctional, i.e. by using the 
// Babuska-Miller trick for evaluation
class Lift : public virtual Gascoigne::ResidualFunctional
{
public:
  Lift(std::string config) : ResidualFunctional()
  {
    __comps.push_back(2);
    __cols.insert(80);

    if (config=="bench2d1")
      {
	__scales.push_back(500.0);
	// exact value from extrapolation on very fine meshes
	ExactValue() = 0.01061908795852302;
      }
    else if (config=="bench2d2")
      {
	__scales.push_back(20.0);
      }
    else if (config=="bench3d1")
      {
	__scales.push_back(500.0/0.41);
	// Taken from Braack & Richter: Solutions of 3D Navier-Stokes benchmarkproblems with adaptive finite elements
	ExactValue() = 6.88130e-2;
      }
    else if ( (config=="bench3d2") || (config=="bench3d3"))
      {
	__scales.push_back(20.0/0.41);
	// no exact values. The original paper by Schaefer/Turek did not consider the nonstationary 3d case with circular obstacle.
	ExactValue() = 0.0;
      }
    else
      {
	std::cerr << "Configuration " << config << " not known!" << std::endl;
	abort();
      }

    __DD = new Gascoigne::DirichletDataByColor(GetComps(), GetColors(), GetScales());
  }

  std::string GetName() const
  {
    return "Lift";
  }
};

/**
 * measures p(x1) - p(x0), where x0,x1 are two points left and
 * right of the obstacle
 */
class PDiff : public Gascoigne::WeightedPointFunctional
{
public:
  PDiff(std::string config) : WeightedPointFunctional()
  {
    // We add two points:
    if ( (config=="bench2d1") || (config=="bench2d2") )
      {
	_v2d.push_back(Gascoigne::Vertex2d(0.15,0.2)); // coordinate
	_comps.push_back(0);                           // comp 0 = Pressure
	_weights.push_back(1.0);                       // weight 
	
	_v2d.push_back(Gascoigne::Vertex2d(0.25,0.2)); // coordinate
	_comps.push_back(0);                           // comp 0 = Pressure
	_weights.push_back(-1.0);                      // weight

	// Exact value taken from:
	// Nabh: On High Order Methods for the Stationary Incompressible Navier-Stokes Equations; University of Heidelberg; Preprint 42/98, 1998
	// https://ganymed.math.uni-heidelberg.de/Paper/Preprint1998-14.pdf
	ExactValue() = 0.11752016697;
      }
    else if ( (config=="bench3d1") || (config=="bench3d2") || (config=="bench3d3") )
      {
	_v3d.push_back(Gascoigne::Vertex3d(0.15,0.2,0.2)); // coordinate
	_comps.push_back(0);                           // comp 0 = Pressure
	_weights.push_back(1.0);                       // weight 
	
	_v3d.push_back(Gascoigne::Vertex3d(0.25,0.2,0.2)); // coordinate
	_comps.push_back(0);                           // comp 0 = Pressure
	_weights.push_back(-1.0);                      // weight	
      }
    


  }

  std::string GetName() const
  {
    return "pressure drop";
  }


};
  

// Computes the squared L2-Norm of the divergence of the velocity field 
template<int DIM>
class Divergence : public Gascoigne::DomainFunctional
{
public:
  double J(const Gascoigne::FemFunction& U, const Gascoigne::Vertex<DIM>& v) const
  {
    double div = 0;
    for (int i=0;i<DIM;++i)
      div += U[i+1][i+1];
    return div*div;
  }

  std::string GetName() const
  {
    return "Divergence";
  }
};

// Computes the kinematic Energy 
template<int DIM>
class Energy : public Gascoigne::DomainFunctional
{
public:
  double J(const Gascoigne::FemFunction& U, const Gascoigne::Vertex<DIM>& v) const
  {
    double en = 0;
    for (int i=0;i<DIM;++i)
      en += U[i+1].m() * U[i+1].m();
    return en;
  }

  std::string GetName() const
  {
    return "Energy";
  }
};


// Defines the drag, evaluated as boundary integral
class BoundaryDrag : public Gascoigne::BoundaryFunctional
{
  const Gascoigne::NavierStokesData& data;

  double weight;
  
public:

  // copy the problem data
  BoundaryDrag(const Gascoigne::NavierStokesData& dat,
	       const std::string& config) : BoundaryFunctional(), data(dat)
  {
    if (config=="bench2d1")
      {
	weight = 500.0;
	ExactValue() = 5.579532446386295;
      }
    else if (config=="bench3d1")
      {
	weight = 500.0/0.41;
	ExactValue() = 5.579532446386295;
      }
    else if ( (config=="bench3d2") || (config=="bench3d3") )
      {
	weight = 20.0/0.41;
	ExactValue() = 0.0;
      }
    else if (config=="bench2d2")
      weight = 20.0;
    else abort();
  }

  std::string GetName() const
  {
    return "Boundary-Drag";
  }


  // evaluates:
  // [ nu * (nabla u + nabla u^T) n - pn ] * ex   where ex = (1,0) (2d) or ex = (1,0,0) (3d)
  double J(const Gascoigne::FemFunction& U, const Gascoigne::Vertex2d& v,
	   const Gascoigne::Vertex2d& n, int color) const
  {
    return - weight *
      ( data.visc * (2.0*U[1].x()*n.x() + (U[1].y()+U[2].x()) * n.y())
	- U[0].m()*n.x() );
  }
  double J(const Gascoigne::FemFunction& U, const Gascoigne::Vertex3d& v,
	   const Gascoigne::Vertex3d& n, int color) const
  {
    return - weight *
      ( data.visc * (2.0*U[1].x()*n.x() +
		     (U[1].y()+U[2].x()) * n.y() +
		     (U[1].z() + U[3].x()) * n.z())
	- U[0].m()*n.x());
  }
};

// Defines the lift, evaluated as boundary integral
class BoundaryLift : public Gascoigne::BoundaryFunctional
{
  const Gascoigne::NavierStokesData& data;

  double weight;
  
public:

  // copy the problem data
  BoundaryLift(const Gascoigne::NavierStokesData& dat,
	       const std::string& config) : BoundaryFunctional(), data(dat)
  {
    if (config=="bench2d1")
      {
	weight = 500.0;
	ExactValue() = 0.01061908795852302;
      }
    else if (config=="bench3d1")
      {
	weight = 500.0/0.41;
	ExactValue() = 6.88130e-2;
      }
    else if ( (config=="bench3d2") || (config=="bench3d3") )
      {
	weight = 20.0/0.41;
	ExactValue() = 0;
      }
    else if (config=="bench2d2")
      weight = 20.0;
    else abort();
  }

  std::string GetName() const
  {
    return "Boundary-Lift";
  }


  // evaluates:
  // [ nu * (nabla u + nabla u^T) n - pn ] * ey   where ey = (0,1) (2d) or ex = (0,1,0) (3d)
  double J(const Gascoigne::FemFunction& U, const Gascoigne::Vertex2d& v,
	   const Gascoigne::Vertex2d& n, int color) const
  {
    return - weight *
      (data.visc * ( (U[2].x() + U[1].y())*n.x() + 2.0*U[2].y() * n.y())
       - U[0].m()*n.y());
  }
  double J(const Gascoigne::FemFunction& U, const Gascoigne::Vertex3d& v,
	   const Gascoigne::Vertex3d& n, int color) const
  {
    return - weight *
      (data.visc * ((U[2].x() + U[1].y()) * n.x() +
		    2.0*U[2].y() * n.y() + 
		    (U[2].z() + U[3].y()) * n.z())
       - U[0].m()*n.y() );
  }
};



/*----------------------------   functionals.h     ---------------------------*/
/* end of #ifndef __functionals_H */
#endif
/*----------------------------   functionals.h     ---------------------------*/
