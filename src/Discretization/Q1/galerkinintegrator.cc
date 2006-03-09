#include  "galerkinintegrator.h"

using namespace std;

/* ----------------------------------------- */

namespace Gascoigne
{
template<int DIM>
GalerkinIntegrator<DIM>::GalerkinIntegrator() : BasicIntegrator(),
  IFF(NULL), IFE(NULL), IFB(NULL), IFM(NULL)
{
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::BasicInit()
{
  if (DIM==2)
    {
      if (!FormFormulaPointer())     FormFormulaPointer() = new QuadGauss4;
      if (!ErrorFormulaPointer())    ErrorFormulaPointer() = new QuadGauss9;
      if (!BoundaryFormulaPointer()) BoundaryFormulaPointer() = new LineGauss2;
      if (!MassFormulaPointer())     MassFormulaPointer() = new QuadGauss4;;
//      if (!MassFormulaPointer())     MassFormulaPointer() = new QuadTrapez;
    }
  else if (DIM==3)
    {
      if (!FormFormulaPointer())     FormFormulaPointer() = new HexGauss8;
      if (!ErrorFormulaPointer())    ErrorFormulaPointer() = new HexGauss27;
      if (!BoundaryFormulaPointer()) BoundaryFormulaPointer() = new QuadGauss4;
      if (!MassFormulaPointer())     MassFormulaPointer() = new HexGauss8;
//      if (!MassFormulaPointer())     MassFormulaPointer() = new HexTrapez;
    }
  assert(FormFormulaPointer());
  assert(ErrorFormulaPointer());
  assert(BoundaryFormulaPointer());
  assert(MassFormulaPointer());
}

/* ----------------------------------------- */

template<int DIM>
GalerkinIntegrator<DIM>::~GalerkinIntegrator<DIM>()
{
  if(FormFormulaPointer()){delete FormFormulaPointer();}
  if(ErrorFormulaPointer()){delete ErrorFormulaPointer();}
  if(BoundaryFormulaPointer()){delete BoundaryFormulaPointer();}
  if(MassFormulaPointer()){delete MassFormulaPointer();}
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::Rhs(const DomainRightHandSide& f, LocalVector& F, const FemInterface& FEM, 
    const LocalNodeData& Q, const LocalCellData& QC) const
{
  F.ReInit(f.GetNcomp(),FEM.n());

  BasicIntegrator::universal_point(_QCH,QC);
  f.SetCellData(_QCH);

  const IntegrationFormulaInterface& IF = *FormFormula();

  F.zero();
  Vertex<DIM> x, xi;

  for (int k=0; k<IF.n(); k++)
    {
		IF.xi(xi,k);
		FEM.point(xi);
		double vol = FEM.J();
		double weight  = IF.w(k) * vol;
		BasicIntegrator::universal_point(FEM,_QH,Q);
		f.SetFemData(_QH);
		FEM.x(x);
		for (int i=0;i<FEM.n();i++)
			{
			FEM.init_test_functions(_NN,weight,i);
			f(F.start(i),_NN,x);
			}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::BoundaryRhs(const BoundaryRightHandSide& f, LocalVector& F, const FemInterface& FEM, 
    int ile, int col, const LocalNodeData& Q, const LocalCellData& QC) const
{
  F.ReInit(f.GetNcomp(),FEM.n());

  BasicIntegrator::universal_point(_QCH,QC);
  f.SetCellData(_QCH);

  const IntegrationFormulaInterface& IF = *BoundaryFormula();

  F.zero();
  Vertex<DIM> x, n;
  Vertex<DIM-1> xi;

  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point_boundary(ile,xi);
      BasicIntegrator::universal_point(FEM,_QH,Q);
      f.SetFemData(_QH);
      FEM.x(x);
      FEM.normal(n);
      double  h = FEM.G();
      double  weight = IF.w(k)*h;
      for (int i=0;i<FEM.n();i++)
	{
	  FEM.init_test_functions(_NN,weight,i);
	  f(F.start(i),_NN,x,n,col);
	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::Form(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
    const LocalNodeData& Q, const LocalCellData& QC) const
{
  F.ReInit(EQ.GetNcomp(),FEM.n());

  BasicIntegrator::universal_point(_QCH,QC);
  EQ.SetCellData(_QCH);

  const IntegrationFormulaInterface& IF = *FormFormula();

  F.zero();
  Vertex<DIM> x, xi;

  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,_UH,U);
      BasicIntegrator::universal_point(FEM,_QH,Q);
      FEM.x(x);
      EQ.SetFemData(_QH);
      EQ.point(h,_UH,x);
      for (int i=0;i<FEM.n();i++)
	{
	  FEM.init_test_functions(_NN,weight,i);
	  EQ.Form(F.start(i),_UH,_NN);
	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::AdjointForm(const Equation& EQ, LocalVector& F, const FemInterface& FEM, const LocalVector& Z, 
    const LocalNodeData& Q, const LocalCellData& QC) const
{
  F.ReInit(EQ.GetNcomp(),FEM.n());

  BasicIntegrator::universal_point(_QCH,QC);
  EQ.SetCellData(_QCH);

  const IntegrationFormulaInterface& IF = *FormFormula();

  F.zero();
  Vertex<DIM> x, xi;

  _NNN.resize(FEM.n());
  EntryMatrix E;
  E.SetDimensionDof(FEM.n(),FEM.n());
  E.SetDimensionComp(Z.ncomp(),Z.ncomp());
  E.resize();
  E.zero();

  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,_UH,Z);
      BasicIntegrator::universal_point(FEM,_QH,Q);
      FEM.x(x);
      // EQ.pointmatrix(h,_QH["u"],_QH,x);
      EQ.pointmatrix(h,_QH["u"],x);
      double sw = sqrt(weight);
      for (int i=0; i<FEM.n(); i++)
	{
	  FEM.init_test_functions(_NNN[i],sw,i);
	}
      for (int j=0; j<FEM.n(); j++)
	{
	  for (int i=0; i<FEM.n(); i++)
	    {
	      E.SetDofIndex(j,i);
	      EQ.Matrix(E,_QH["u"],_NNN[i],_NNN[j]);
	    }
	}
    }
  for (int i=0; i<FEM.n(); i++)
    {
      for (int c=0; c<Z.ncomp(); c++)
	{
	  double sum = 0.;
	  for (int j=0; j<FEM.n(); j++)
	    {
	      for (int d=0; d<Z.ncomp(); d++)
		{
		  sum += E(j,i,d,c)*Z(j,d);
		}
	    }
	  F(i,c) += sum;
	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::BoundaryForm(const BoundaryEquation& BE, LocalVector& F, const FemInterface& FEM, const LocalVector& U, 
    int ile, int col, const LocalNodeData& Q, const LocalCellData& QC) const
{
  F.ReInit(BE.GetNcomp(),FEM.n());

  BasicIntegrator::universal_point(_QCH,QC);
  BE.SetCellData(_QCH);

  const IntegrationFormulaInterface& IF = *BoundaryFormula();

  F.zero();
  Vertex<DIM> x,n;
  Vertex<DIM-1> xi;

  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point_boundary(ile,xi);
      BasicIntegrator::universal_point(FEM,_UH,U);
      BasicIntegrator::universal_point(FEM,_QH,Q);
      FEM.x(x);
      FEM.normal(n);
      double  h = FEM.G();
      double  weight = IF.w(k)*h;
      BE.SetFemData(_QH);
      BE.pointboundary(h,_UH,x,n);
      for (int i=0;i<FEM.n();i++)
        {
          FEM.init_test_functions(_NN,weight,i);
          BE.Form(F.start(i),_UH,_NN,col);
        }
    }
}

/* ----------------------------------------- */

template<int DIM>
double GalerkinIntegrator<DIM>::MassMatrix(EntryMatrix& E, const FemInterface& FEM) const
{
  _NNN.resize(FEM.n());
  E.SetDimensionDof(FEM.n(),FEM.n());
  int ncomp=1;
  E.SetDimensionComp(ncomp,ncomp);
  E.resize();
  E.zero();
  const IntegrationFormulaInterface& IF = *MassFormula();

  Vertex<DIM> x, xi;
  double omega = 0.;
  for (int k=0; k<IF.n(); k++)
    {
		IF.xi(xi,k);
		FEM.point(xi);
		double vol = FEM.J();
		double weight  = IF.w(k) * vol;
		omega += weight;
		for (int i=0;i<FEM.n();i++)
			{
			FEM.init_test_functions(_NNN[i],1.,i);
			}
		for (int i=0;i<FEM.n();i++)
			{
			for (int j=0;j<FEM.n();j++)
				{
				E.SetDofIndex(i,j);
				E(0,0) += weight * _NNN[j].m()*_NNN[i].m();
				}
			}
    }
  return omega;
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::Matrix(const Equation& EQ, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
    const LocalNodeData& Q, const LocalCellData& QC) const
{
  _NNN.resize(FEM.n());
  E.SetDimensionDof(FEM.n(),FEM.n());
  E.SetDimensionComp(U.ncomp(),U.ncomp());
  E.resize();
  E.zero();

  BasicIntegrator::universal_point(_QCH,QC);
  EQ.SetCellData(_QCH);

  const IntegrationFormulaInterface& IF = *FormFormula();

  Vertex<DIM> x, xi;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double h  = Volume2MeshSize(vol);
      double weight  = IF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,_UH,U);
      BasicIntegrator::universal_point(FEM,_QH,Q);
      FEM.x(x);
      EQ.SetFemData(_QH);
      EQ.pointmatrix(h,_UH,x);

      double sw = sqrt(weight);
      for (int i=0;i<FEM.n();i++)
				{
				FEM.init_test_functions(_NNN[i],sw,i);
				}
      for (int j=0; j<FEM.n(); j++)
				{
				for (int i=0; i<FEM.n(); i++)
					{
					E.SetDofIndex(i,j);
					EQ.Matrix(E,_UH,_NNN[j],_NNN[i]);
					}
				}
    }
}

/*-----------------------------------------------------------*/

template<int DIM>
void GalerkinIntegrator<DIM>::BoundaryMatrix (const BoundaryEquation& BE, EntryMatrix& E, const FemInterface& FEM, const LocalVector& U, 
    int ile, int col, const LocalNodeData& Q, const LocalCellData& QC) const
{
  _NNN.resize(FEM.n());
  E.SetDimensionDof(FEM.n(),FEM.n());
  E.SetDimensionComp(U.ncomp(),U.ncomp());
  E.resize();
  E.zero();

  BasicIntegrator::universal_point(_QCH,QC);
  BE.SetCellData(_QCH);

  const IntegrationFormulaInterface& IF = *BoundaryFormula();

  Vertex<DIM> x,n;
  Vertex<DIM-1> xi;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point_boundary(ile,xi);
      BasicIntegrator::universal_point(FEM,_UH,U);
      BasicIntegrator::universal_point(FEM,_QH,Q);
      FEM.x(x);
      FEM.normal(n);
      double  h = FEM.G();
      double  weight = IF.w(k)*h;
      BE.SetFemData(_QH);
      BE.pointmatrixboundary(h,_UH,x,n);
      double sw = sqrt(weight);
      for (int i=0;i<FEM.n();i++)
        {
          FEM.init_test_functions(_NNN[i],sw,i);
        }
      for (int j=0;j<FEM.n();j++)
        {
          for (int i=0;i<FEM.n();i++)
            {
              E.SetDofIndex(i,j);
              BE.Matrix(E,_UH,_NNN[j],_NNN[i],col);
            }
        }
    }
}

/*-----------------------------------------------------------*/

template<int DIM>
void GalerkinIntegrator<DIM>::RhsPoint
(LocalVector& F, const FemInterface& E, const Vertex<DIM>& p, int comp) const
{
  F.zero();

  E.point(p);
  for (int i=0; i<E.n(); i++)
    {
      E.init_test_functions(_NN,1.,i);
      F(i,comp) += _NN.m();
    }
}

/* ----------------------------------------- */
template<int DIM>
void GalerkinIntegrator<DIM>::DiracRhsPoint(LocalVector& b, const FemInterface& E, const Vertex<DIM>& p, const DiracRightHandSide& DRHS, int j, 
    const LocalNodeData& Q, const LocalCellData& QC) const
{
  b.zero();

  BasicIntegrator::universal_point(_QCH,QC);
  DRHS.SetCellData(_QCH);

  Vertex<DIM> x;
  E.point(p);     
  E.x(x);
  BasicIntegrator::universal_point(E,_QH,Q);
  DRHS.SetFemData(_QH);

  for (int i=0; i<E.n(); i++)
    {
      E.init_test_functions(_NN,1.,i);
      DRHS.operator()(j,b.start(i),_NN,x);
    }
}

/* ----------------------------------------- */
template<int DIM>
double GalerkinIntegrator<DIM>::ComputePointValue(const FemInterface& E, const Vertex<DIM>& p, const LocalVector& U, int comp) const
{
  E.point(p);
  BasicIntegrator::universal_point(E,_UH,U);

  return _UH[comp].m();
}

/* ----------------------------------------- */

template<int DIM>
double GalerkinIntegrator<DIM>::ComputeBoundaryFunctional(const BoundaryFunctional& F, const FemInterface& FEM, int ile, 
    int col, const LocalVector& U, const LocalNodeData& Q, const LocalCellData& QC) const
{
  const IntegrationFormulaInterface& IF = *BoundaryFormula();

  Vertex<DIM>   x, n;
  Vertex<DIM-1> xi;


  double j = 0.;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point_boundary(ile,xi);
      double h = FEM.G();
      double weight  = IF.w(k) * h;
      BasicIntegrator::universal_point(FEM,_UH,U);
      BasicIntegrator::universal_point(FEM,_QH,Q);
      F.SetFemData(_QH);

      FEM.x(x);
      FEM.normal(n);
      // FEM.normal(n);
      j += weight * F.J(_UH,x,n,col);
    }
  return j;
}

/* ----------------------------------------- */

template<int DIM>
double GalerkinIntegrator<DIM>::ComputeDomainFunctional(const DomainFunctional& F, const FemInterface& FEM, const LocalVector& U, 
    const LocalNodeData& Q, const LocalCellData& QC) const
{
  BasicIntegrator::universal_point(_QCH,QC);
  F.SetCellData(_QCH);

  const IntegrationFormulaInterface& IF = *FormFormula();

  Vertex<DIM> x, xi;
  double j = 0.;
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double weight  = IF.w(k) * vol;
      BasicIntegrator::universal_point(FEM,_UH,U);
      BasicIntegrator::universal_point(FEM,_QH,Q);
      FEM.x(x);
      F.SetFemData(_QH);
      j += weight * F.J(_UH,x);
    }
  return j;
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::EvaluateCellRightHandSide(LocalVector& F, const DomainRightHandSide& CF,const FemInterface& FEM, 
    const LocalNodeData& Q, const LocalCellData& QC) const
{
  F.ReInit(CF.GetNcomp(),1);

  BasicIntegrator::universal_point(_QCH,QC);
  CF.SetCellData(_QCH);

  const IntegrationFormulaInterface& IF = *FormFormula();
  Vertex<DIM> x, xi;

  F.zero();
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double weight  = IF.w(k) * vol;

      BasicIntegrator::universal_point(FEM,_QH,Q);
      FEM.x(x);
      CF.SetFemData(_QH);

      _NN.zero();
      _NN.m() = weight;

      CF(F.start(0),_NN,x);
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::ErrorsByExactSolution(LocalVector& dst, const FemInterface& FE, const ExactSolution& ES, const LocalVector& U, 
    const LocalNodeData& Q, const LocalCellData& QC) const
{
  BasicIntegrator::universal_point(_QCH,QC);
  ES.SetCellData(_QCH);

  const IntegrationFormulaInterface& IF = *ErrorFormula();

  Vertex<DIM> x, xi;
  dst.zero();
  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FE.point(xi);
      BasicIntegrator::universal_point(FE,_UH,U);
      double vol = FE.J();
      double weight = IF.w(k) * vol;

      FE.x(x);
      for(int c=0; c<U.ncomp(); c++)
	{
	  _UH[c].m() -= ES(c,x);
	  _UH[c].x() -= ES.x(c,x);
	  _UH[c].y() -= ES.y(c,x);
	  if (DIM==3) _UH[c].z() -= ES.z(c,x);
	}
      for(int c=0; c<U.ncomp(); c++) 
	{
	  // L2 Norm
	  dst(0,c) += weight * _UH[c].m() * _UH[c].m();
	  // H1 Seminorm
	  double a = _UH[c].x()*_UH[c].x() + _UH[c].y()*_UH[c].y();
	  if (DIM==3) a += _UH[c].z() * _UH[c].z();
	  dst(1,c) += weight * a;
	  // L8 Norm
	  dst(2,c) = Gascoigne::max(dst(2,c),fabs(_UH[c].m()));
	}
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::IntegrateMassDiag(DoubleVector& F, const FemInterface& FEM) const 
{
  F.resize(FEM.n());

  const IntegrationFormulaInterface& IF = *FormFormula();

  F.zero();
  Vertex<DIM> x, xi;

  for (int k=0; k<IF.n(); k++)
    {
      IF.xi(xi,k);
      FEM.point(xi);
      double vol = FEM.J();
      double weight  = IF.w(k) * vol;
      for (int i=0;i<FEM.n();i++)
      {
	  FEM.init_test_functions(_NN,1.,i);
	  F[i] += weight*_NN.m()*_NN.m();
      }
    }
}

/* ----------------------------------------- */

template<int DIM>
void GalerkinIntegrator<DIM>::IntegrateBoundaryMassDiag(DoubleVector& F, const FemInterface& FEM, int ile, int col) const 
{
    F.resize(FEM.n());

    const IntegrationFormulaInterface& IF = *BoundaryFormula();

    F.zero();
    Vertex<DIM> x,n;
    Vertex<DIM-1> xi;

    for (int k=0; k<IF.n(); k++)
    {
	IF.xi(xi,k);
	FEM.point_boundary(ile,xi);
	
      	double  h = FEM.G();
	double  weight = IF.w(k)*h;
      
	for (int i=0;i<FEM.n();i++)
        {
	    FEM.init_test_functions(_NN,1.,i);
	    F[i] += _NN.m()*_NN.m()*weight;
        }
    }
}

/* ----------------------------------------- */

template class GalerkinIntegrator<2>;
template class GalerkinIntegrator<3>;
}
