#include "edgeintegrator.h"

namespace Gascoigne
{

  template<int DIM>
  void EdgeIntegrator<DIM>::EdgeForm(const EdgeEquation& EQ, const DGEdge& edge,
				     LocalVector& F1, LocalVector& F2,
				     const FemInterface& FEM1, const FemInterface& FEM2,
				     const LocalVector& U1, const LocalVector& U2,
				     const LocalData& Q, const LocalData& QC) const
  {
    abort();
    
//     assert(FEM1.n() == FEM2.n());
//     assert(U1.ncomp() == U2.ncomp());
    
//     // first, only interiour edges
//     assert(edge.master()!=-1);
//     assert(edge.slave()!=-1);
    
//     F1.ReInit(U1.ncomp(),FEM1.n()); F1.zero(); 
//     F2.ReInit(U2.ncomp(),FEM2.n()); F2.zero(); 
    
//     // no cell-data on edge... edge data? Or do we need two different cell-data's???
// #warning cell data oin EdgeForm
//     // BasicIntegrator::universal_point(_QCH,QC);
//     // EQ.SetCellData(_QCH);
    
//     //const IntegrationFormulaInterface& IF = *FormFormula();
//     //    Edgeformula...
    
    
//     Vertex<DIM> x,x1,x2, xi1, xi2, n,n1,n2;

//     TestFunction N1,N2;
//     FemFunction UH1,UH2;

    
    
//     for (int k=0; k<IF.n(); k++)
//       {
// 	IF.xi(xi1,k); // IF on edge. 
// 	IF.xi(xi2,IF.n()-k-1); // we need to, as two separate edges are considered, one running backwards
	
// 	FEM1.point_boundary(xi1, edge.local_master());
// 	FEM2.point_boundary(xi2, edge.local_slave());
	
// 	double vol1 = FEM1.J();
// 	double vol2 = FEM2.J();
// 	double h1  = Volume2MeshSize(vol1);
// 	double h2  = Volume2MeshSize(vol2);
// 	assert(fabs(h1-h2)<1.e-14);
	
// 	double weight  = IF.w(k) * vol1;
	
// 	BasicIntegrator::universal_point(FEM1,UH1,U1);
// 	BasicIntegrator::universal_point(FEM2,UH2,U2);
	
// #warning femdata!
// 	//      BasicIntegrator::universal_point(FEM,_QH,Q);
// 	//EQ.SetFemData(_QH);
	
// 	FEM1.x(x1);
// 	FEM2.x(x2);
// 	x = x1-x2; assert(x.norm()<1.e-13);
	
// 	FEM1.normal(n1);
// 	FEM2.normal(n2);
// 	n = n1+n2; assert(n.norm()<1.e-13);
	
// 	EQ.point_edge(h1,_UH1,_UH2,n1);
	
// 	for (int i=0;i<FEM1.n();i++)
// 	  {
// 	    FEM1.init_test_functions(N1,weight,i);
// 	    FEM2.init_test_functions(N2,weight,i);
// 	    EQ.EdgeForm(F1.start(i),F2.start(i),UH1,UH2,N1,N2);
// 	  }
//       }

  }

  template class EdgeIntegrator<2>;
  template class EdgeIntegrator<3>;
}
