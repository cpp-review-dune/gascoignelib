/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2011 by the Gascoigne 3D authors
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


#ifndef  __CellDiscretization_h
#define  __CellDiscretization_h

#include  "basicdiscretization.h"
#include  "hierarchicalmesh.h"
#include  "hnstructureinterface.h"
#include  "edgeinfocontainer.h"

#include  "transformation2d.h"
#include  "transformation3d.h"
#include  "baseq1.h"
#include  "finiteelement.h"

namespace Gascoigne
{

  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments CellDiscretization
  ///
  ///
  /////////////////////////////////////////////

  template<int DIM>
  class FemType
  { };
  template<>
  class FemType<2> 
  {
  public:
    typedef Transformation2d<BaseQ1<2> >           TransQ1;
    typedef FiniteElement<2,1,TransQ1,BaseQ1<2> >  FFF;
  };
  template<>
  class FemType<3> 
  {
  public:
    typedef Transformation3d<BaseQ1<3>  >           TransQ1;
    typedef FiniteElement<3,2,TransQ1,BaseQ1<3> >  FFF;
  };


    
    

  template <int DIM>
  class Q1 : virtual public BasicDiscretization
  {
  private:
    
  protected:


    FemInterface*         __FEM;
    IntegratorInterface*  __INT;

    const FemInterface* GetFem() const               { assert(__FEM); return __FEM; }
    const IntegratorInterface* GetIntegrator() const { assert(__INT); return __INT; }
    IntegratorInterface*& GetIntegratorPointer()     { return __INT; }
    FemInterface*& GetFemPointer()                   { return __FEM; }

    virtual void Transformation(FemInterface::Matrix& T, int iq) const;
    
    double ComputePointValue(const GlobalVector& u, const Vertex<DIM>& p0,int comp) const;

    virtual int GetCellNumber(const Vertex<DIM>& p0, Vertex<DIM>& p, int c0=0) const {
      std::cerr << "\"CellDiscretization::GetCellNumber\" not written!" << std::endl;
      abort();
    }

    // Haengende Knoten aus Q1 alt
    HNStructureInterface*    HN;
    
    IntVector GetLocalIndices(int iq) const {
      return GetMesh()->IndicesOfCell(iq);
    }
    void LocalToGlobal(MatrixInterface& A, EntryMatrix& E, int iq, double s) const;
    
    
    virtual HNStructureInterface* NewHNStructure();
    

    /////
    void Transformation_HM(FemInterface::Matrix& T, const HierarchicalMesh* HM, int iq) const;
    void GlobalToLocal_HM(LocalVector& U, const GlobalVector& u, const HierarchicalMesh* HM, int iq) const;
    void swapIndices(IntVector& indices) const;
  
    virtual void DiracRhsPoint(GlobalVector& f,const DiracRightHandSide& DRHS,const Vertex<DIM>& p0,int i,double s) const;
    

  public:

    //    void BasicInit(const ParamFile* pf);
    

    



    void ConstructInterpolator(MgInterpolatorInterface* I, const MeshTransferInterface* MT);    

    // Threads
    void InitColoring();
    int                            _col_graph_ncol;
    std::vector<std::vector<int> > _col_graph;

    //
    ////  Constructor 
    //

    Q1() : BasicDiscretization(), __FEM(NULL), __INT(NULL), HN(NULL) {}
    ~Q1(){
      if(__FEM) {delete __FEM; __FEM=NULL;}
      if(__INT) {delete __INT; __INT=NULL;}
      if (HN) delete HN; HN  = NULL;
      
    }

    std::string GetName() const {return "Q1";}

    void Structure(SparseStructureInterface* S) const;

    void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
    void AdjointForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
    void BoundaryForm(GlobalVector& f, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
    void Matrix(MatrixInterface& A, const GlobalVector& u, const ProblemDescriptorInterface* PD, double) const;
    void BoundaryMatrix(MatrixInterface& A, const GlobalVector& u, const IntSet& Colors, const BoundaryEquation& BE, double d) const;
    void MassMatrix(MatrixInterface& M) const;
    void BoundaryMassMatrix(MatrixInterface& A, const IntSet& Colors) const;
    void MassForm(GlobalVector& f, const GlobalVector& u, const TimePattern& TP, double s) const;

    void ComputeError(const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const;
    void AssembleError(GlobalVector& eta, const GlobalVector& u, LocalVector& err, const ExactSolution* ES) const;

    void Rhs(GlobalVector& f, const DomainRightHandSide& RHS, double s) const;
    void DiracRhs(GlobalVector& f, const DiracRightHandSide& DRHS, double s) const;

    void BoundaryRhs(GlobalVector& f, const IntSet& Colors,  const BoundaryRightHandSide& NRHS, double s) const;

    void InitFilter(DoubleVector&) const;

    // Functionals
    double ComputeBoundaryFunctional(const GlobalVector& u, const IntSet& Colors, const BoundaryFunctional& BF) const;
    double ComputeDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const;
    double ComputeErrorDomainFunctional(const GlobalVector& u, const DomainFunctional& F) const;

    double ComputePointFunctional(const GlobalVector& u, const PointFunctional& FP) const;

    void EvaluateCellRightHandSide(GlobalVector& f, const DomainRightHandSide& CF, double d) const;
    void EvaluateBoundaryCellRightHandSide(GlobalVector& f,const IntSet& Colors, const BoundaryRightHandSide& CF, double d) const;
    void EvaluateParameterRightHandSide(GlobalVector& f, const DomainRightHandSide& CF, double d) const;
    void EvaluateBoundaryParameterRightHandSide(GlobalVector& f,const IntSet& Colors, const BoundaryRightHandSide& CF, double d) const;

    void InterpolateDomainFunction(GlobalVector& f, const DomainFunction& DF) const;
    void InterpolateCellDomainFunction(GlobalVector& f, const DomainFunction& DF) const;

    virtual nmatrix<double> GetLocalInterpolationWeights() const {
      std::cerr << "\"CellDiscretization::GetLocalInterpolationWeights\" not written!" << std::endl;
      abort();
    }

    void GetVolumes(DoubleVector& a) const;
    void GetAreas(DoubleVector& a, const IntSet& Colors) const;
    void GetMassDiag(DoubleVector& a) const;
    void GetBoundaryMassDiag(DoubleVector& a) const;      



    //// Aus Q1

    
    int n() const {
      return GetMesh()->nnodes();
    }
    int nc() const {
      return GetMesh()->ncells();
    }
    //  const HNStructureQ1* GetHNStructure() const { return HN;}
    
    void ReInit   (const MeshInterface* MP)
    {
      BasicDiscretization::ReInit(MP);
      assert(HN);
      HN->ReInit(MP);
    }
    
    void StrongDirichletVector(GlobalVector& u, const DirichletData& BF, int col, const std::vector<int>& comp, double d) const ;
    
    void StrongDirichletMatrix       (MatrixInterface& A, int col, const std::vector<int>& comp) const;
    void StrongDirichletMatrixOnlyRow(MatrixInterface& A, int col, const std::vector<int>& comp) const;
    void StrongDirichletVectorZero(GlobalVector& u, int col, const std::vector<int>& comp) const;
    void InterpolateSolution(GlobalVector& u, const GlobalVector& uold)const;
    virtual void InterpolateSolutionByPatches(GlobalVector& u, const GlobalVector& uold) const {}
    void HNAverage   (GlobalVector& x) const;
    void HNDistribute(GlobalVector& x) const;
    void HNZero      (GlobalVector& x) const;
    bool HNZeroCheck (const GlobalVector& x) const;

    virtual void EnergyEstimator(EdgeInfoContainerInterface& EIC, DoubleVector& eta, const GlobalVector& u, const Equation& EQ, const DomainRightHandSide* RHS, const std::string & s_energytype, double d_visc) const{ abort(); }
    
    


  };
}

#endif
