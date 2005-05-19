#ifndef __Dwrfem_h
#define __Dwrfem_h

#include "q22d.h"
#include "q23d.h"
#include "baseq1patch.h"
#include "baseq13dpatch.h"
#include "transformation2d.h"
#include "transformation3d.h"
#include "finiteelement.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class DwrFem2d : public Q22d
{
 protected:
  typedef Transformation2d<BaseQ12d>         TransQ1;
  FiniteElement<2,1,TransQ1,BaseQ12dPatch>   LowOrderFem;
  
  void TransformationQ1(FemInterface::Matrix& T, int iq) const;

 public:
  DwrFem2d();
  ~DwrFem2d() {}

  void BasicInit(const ParamFile* paramfile);
};

/*---------------------------------------------------*/

class DwrFem3d : public Q23d
{
 protected:
  typedef Transformation3d<BaseQ13d>         TransQ1;
  FiniteElement<3,2,TransQ1,BaseQ13dPatch>   LowOrderFem;
  
  void TransformationQ1(FemInterface::Matrix& T, int iq) const;

 public:
  DwrFem3d();
  ~DwrFem3d() {}

  void BasicInit(const ParamFile* paramfile);
};

/*---------------------------------------------------*/
/*---------------------------------------------------*/

class DwrFemQ1Q22d : public DwrFem2d
{
 protected:
  void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex2d& p0, int i, double s) const;

 public:
  DwrFemQ1Q22d() : DwrFem2d() {}
  ~DwrFemQ1Q22d() {}

  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
  void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
  void BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const;
  
  void MassMatrix(MatrixInterface& M) const;
};

/*---------------------------------------------------*/

class DwrFemQ1Q23d : public DwrFem3d
{
 protected:
  void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex3d& p0, int i, double s) const;

 public:
  DwrFemQ1Q23d() : DwrFem3d() {}
  ~DwrFemQ1Q23d() {}

  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
  void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
  void BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const;
  
  void MassMatrix(MatrixInterface& M) const;
};

/*---------------------------------------------------*/
/*---------------------------------------------------*/

class DwrFemQ2Q12d : public DwrFem2d
{
 protected:
  void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex2d& p0, int i, double s) const;

 public:
  DwrFemQ2Q12d() : DwrFem2d() {}
  ~DwrFemQ2Q12d() {}

  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
  void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
  void BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const;
 
  void MassMatrix(MatrixInterface& M) const;
};

/*---------------------------------------------------*/

class DwrFemQ2Q13d : public DwrFem3d
{
 protected:
  void DiracRhsPoint(GlobalVector& f, const DiracRightHandSide& DRHS, const Vertex3d& p0, int i, double s) const;

 public:
  DwrFemQ2Q13d() : DwrFem3d() {}
  ~DwrFemQ2Q13d() {}

  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
  void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
  void BoundaryRhs(GlobalVector& f, const IntSet& Colors, const BoundaryRightHandSide& NRHS, double s) const;
  
  void MassMatrix(MatrixInterface& M) const;
};
}

/*---------------------------------------------------*/

#endif 
