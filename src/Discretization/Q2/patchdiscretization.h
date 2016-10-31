/**
 *
 * Copyright (C) 2004, 2005, 2006, 2008, 2011 by the Gascoigne 3D authors
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


#ifndef  __PatchDiscretization_h
#define  __PatchDiscretization_h

#include  "basicdiscretization.h"
#include  "patchmesh.h"
#include  "gascoignemesh.h"
#include  "hnstructureq1.h"


namespace Gascoigne
{

  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments PatchDiscretization

  ///
  ///
  /////////////////////////////////////////////

  class PatchDiscretization : public BasicDiscretization
  {
  protected:

    // Data

    HNStructureQ1*          HN;
    FemInterface*         __FEM;
    IntegratorInterface*  __INT;


    // Access to data

    const FemInterface* GetFem()               const {return __FEM;}
    const IntegratorInterface* GetIntegrator() const {return __INT;}
    IntegratorInterface*& GetIntegratorPointer()     {return __INT;}
    FemInterface*& GetFemPointer()                   {return __FEM;}


    // Init
    void InitColoring() ;
    virtual void Transformation(FemInterface::Matrix& T, int iq) const;
    virtual double compute_element_mean_matrix(int iq, EntryMatrix& E) const;
    virtual int GetPatchNumber(const Vertex2d& p0, Vertex2d& p) const {
      std::cerr << "\"PatchDiscretization::GetPatchNumber\" not written!" << std::endl;
      abort();
    }
    virtual int GetPatchNumber(const Vertex3d& p0, Vertex3d& p) const {
      std::cerr << "\"PatchDiscretization::GetPatchNumber\" not written!" << std::endl;
      abort();
    }
  
    const PatchMesh* GetPatchMesh() const {
      const PatchMesh* MP = dynamic_cast<const PatchMesh*>(GetMesh());
      assert(MP);
      return MP;
    }

    const GascoigneMeshBase*  GetGascoigneMesh() const {
      const GascoigneMeshBase* MP = dynamic_cast<const GascoigneMeshBase*>(GetMesh());
      assert(MP);
      return MP;
    }
  
    nvector<int> GetLocalIndices(int iq) const {
      return *GetPatchMesh()->IndicesOfPatch(iq);
    }
  
    // DOF's

    void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const
    {
      nvector<int> indices = GetLocalIndices(iq);
      HN->CondenseHanging(E,indices);
      nvector<int>::const_iterator  start = indices.begin();
      nvector<int>::const_iterator  stop  = indices.end();
      A.entry(start,stop,E,s);
    }
  


    // Computations

    double ComputePointValue(const GlobalVector& u, const Vertex2d& p0,int comp) const;
    double ComputePointValue(const GlobalVector& u, const Vertex3d& p0,int comp) const; 

    virtual void DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex2d& p0,int i,double s) const;
    virtual void DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex3d& p0,int i,double s) const;

    void GlobalToLocalCell(LocalVector& U, const GlobalVector& u, int iq) const;

  public:

    //
    ///  Constructor 
    //
    PatchDiscretization() : BasicDiscretization(), HN(NULL), __FEM(NULL), __INT(NULL) {}
    ~PatchDiscretization(){
      assert(HN==NULL);
      if(__FEM) {delete __FEM; __FEM=NULL;}
      if(__INT) {delete __INT; __INT=NULL;}
    }

    void ReInit(const MeshInterface* MP)
    {
      BasicDiscretization::ReInit(MP);
      HN->ReInit(MP);
    }


    // Access
  

    std::string GetName() const {return "PatchDiscretization";}
    int n() const  { return GetMesh()->nnodes(); }
    int nc() const { return GetMesh()->ncells(); }
    int n_withouthanging() const { return GetMesh()->nnodes()-HN->nhnodes(); }
  


    // DOF's

    void Structure(SparseStructureInterface* S) const;
    void Interpolate(GlobalVector& u, const DomainInitialCondition& U) const
    {
      std::cerr << "\"Q2::Interpolate\" not written!" << std::endl;
      abort();
    }
    void InterpolateSolution(GlobalVector& u, const GlobalVector& uold)const;


    // Boundary

    void StrongDirichletMatrix(MatrixInterface& A, int col, const std::vector<int>& comp) const;
    void StrongDirichletMatrixOnlyRow(MatrixInterface& A, int col, const std::vector<int>& comp) const;
    void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp, double d) const;
    void StrongPeriodicVector(GlobalVector& u, const PeriodicData& BF, int col, const std::vector<int>& comp, double d) const;
    void StrongDirichletVectorZero(GlobalVector& u, int col, const std::vector<int>& comp) const;
    



    

    // Integrate


    void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
    void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
    void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
    void Matrix(MatrixInterface& A, const GlobalVector& u, const ProblemDescriptorInterface* PD, double) const;
    void BoundaryMatrix(MatrixInterface& A, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
    void MassMatrix(MatrixInterface& M) const;
    void MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const;

    void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const;

    void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
    void DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const;
  
    void BoundaryRhs(GlobalVector& f, const IntSet& Colors,  const BoundaryRightHandSide& NRHS, double s) const;

    void InitFilter(nvector<double>&) const;

    // Functionals
    double ComputeBoundaryFunctional(const GlobalVector& u, const IntSet& Colors, const BoundaryFunctional& BF) const;
    double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const;
    double ComputePointFunctional(const GlobalVector& u, const PointFunctional& FP) const;

    void EvaluateParameterRightHandSide(GlobalVector& f, const DomainRightHandSide& CF, double d) const;
    void EvaluateBoundaryParameterRightHandSide(GlobalVector& f,const IntSet& Colors, const BoundaryRightHandSide& CF, double d) const;

    virtual nmatrix<double> GetLocalInterpolationWeights(int iq) const {
      std::cerr << "\"PatchDiscretization::GetLocalInterpolationWeights\" not written!" << std::endl;
      abort();
    }

    /// Hanging Nodes
    void HNAverage   (GlobalVector& x) const
    { HN->Average(x); }
    void HNDistribute(GlobalVector& x) const
    { HN->Distribute(x); }
    void HNZero      (GlobalVector& x) const
    { HN->Zero(x); }
    bool HNZeroCheck (const GlobalVector& x) const
    { return HN->ZeroCheck(x); }
  
  

  };
}

#endif
