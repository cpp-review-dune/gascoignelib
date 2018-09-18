#include "dgsolver.h"
#include "compose_name.h"
#include "gascoignemesh2d.h"
#include <array>

namespace Gascoigne
{

  void DGSolver::Visu(const std::string &name,
                      const VectorInterface &gu,
                      int i) const
  {
    assert(gu.GetType() == "node");
    std::string filename = name;
    compose_name(filename, i);
    filename += ".vtk";
    std::cout << filename << std::endl;

    std::ofstream VTKOUT(filename.c_str());
    ////////////////////////////////////////////////// Header
    VTKOUT << "# vtk DataFile Version 2.0" << std::endl
           << "Gascoigne VTK-Output Discontinuous Galerkin" << std::endl
           << "ASCII" << std::endl;

    ////////////////////////////////////////////////// Grid
    const GascoigneMesh2d *M = dynamic_cast<const GascoigneMesh2d *>(GetMesh());
    assert(M);

    // DG1 / DG2 ?
    const GlobalVector& U = GetGV(gu);
    int dgtype = 0;
    if (U.n() == M->ncells()*4) dgtype = 1;
    if (U.n() == M->ncells()*9) dgtype = 2;
    assert(dgtype>0);

    if (dgtype == 1) ////////////////////////////////////////////////// DG1
      {
	// points // DG ... stupid? multiple times...
	VTKOUT << "DATASET UNSTRUCTURED_GRID" << std::endl
	       << "POINTS " << M->ncells()*4 << " DOUBLE" << std::endl;
	int itoj[4]= {0,1,3,2};
	for (int c = 0; c < M->ncells(); ++c)
	  {
	    const IntVector& ioc = M->IndicesOfCell(c);
	    for (int i = 0; i < ioc.size();++i)
	      VTKOUT << M->vertex2d(ioc[itoj[i]]) << " 0" << std::endl;
	  }
	VTKOUT << std::endl << std::endl;

	// cells
	VTKOUT << "CELLS " << M->ncells() << " " << M->ncells() * 5 << std::endl;
	for (int c = 0; c < M->ncells(); ++c)
	  {
	    VTKOUT << "4 ";
	    for (int i=0;i<4;++i) VTKOUT << c*4+i << " ";
	    VTKOUT << std::endl;
	  }
	VTKOUT << std::endl << std::endl;
	
	VTKOUT << "CELL_TYPES " << M->ncells() << std::endl;
	for (int i = 0; i < M->ncells(); ++i)
	  VTKOUT << "9 ";
	VTKOUT << std::endl << std::endl;
	
	// data
	assert(U.n() == M->ncells() * 4);
	VTKOUT << "POINT_DATA " << M->ncells() * 4 << std::endl;
	for (int nc=0;nc<U.ncomp();++nc)
	  {
	    VTKOUT << "SCALARS u" << nc << " DOUBLE" << std::endl
		   << "LOOKUP_TABLE DEFAULT"  << std::endl;
	    for (int c=0;c<M->ncells();++c)
	      for (int i=0;i<4;++i)
		VTKOUT << U(4*c+itoj[i],nc) << std::endl;
	    VTKOUT << std::endl;
	  }
      }
    else if (dgtype == 2) //////////////////////////////////////// DG 2
      { 
	////////////////////////////////////////////////// split quad into 4 smaller ones
	VTKOUT << "DATASET UNSTRUCTURED_GRID" << std::endl
	       << "POINTS " << M->ncells()*9 << " DOUBLE" << std::endl;
	int itop[4] = {0,2,6,8}; // corners of patch
	int itoj[4]= {0,1,3,2}; // counter-clockwise / lexico
	int iinp[4][4] = {{0,1,4,3},{1,2,5,4},{3,4,7,6},{4,5,8,7}}; // cells in patch
	int avge[4][3] = {{1,0,2},{5,2,8},{3,0,6},{7,6,8}}; //edges
	for (int c = 0; c < M->ncells(); ++c)
	  {
	    const IntVector& ioc = M->IndicesOfCell(c);
	    assert(ioc.size()==4);
	    std::array<Vertex2d,9> vop;  // vertices of patch
	
	    // corners
	    for (int i = 0; i < ioc.size();++i)
	      vop[itop[i]] = M->vertex2d(ioc[i]);
	    // middle
	    vop[4].zero(); 
	    for (int i = 0; i < ioc.size();++i)
	      vop[4].add(0.25,vop[itop[i]]);
	    // edges
	    for (int e=0;e<4;++e)
	      {
		vop[avge[e][0]].zero();
		for (int i=0;i<2;++i)
		  vop[avge[e][0]].add(0.5, vop[avge[e][i+1]]);
	      }

	    // print
	    for (int i=0;i<9;++i)
	      {
		VTKOUT << vop[i]  << " 0" << std::endl;
	      }
	  }
	VTKOUT << std::endl << std::endl;

	// cells
	VTKOUT << "CELLS " << 4*M->ncells() << " " << 4*M->ncells() * 5 << std::endl;
	for (int c = 0; c < M->ncells(); ++c)
	  {
	    for (int cc=0;cc<4;++cc)
	      {
		VTKOUT << "4 ";
		for (int i=0;i<4;++i) VTKOUT << c*9+iinp[cc][i] << " ";
		VTKOUT << std::endl;
	      }
	  }
	VTKOUT << std::endl << std::endl;
	
	VTKOUT << "CELL_TYPES " << 4*M->ncells() << std::endl;
	for (int i = 0; i < 4*M->ncells(); ++i)
	  VTKOUT << "9 ";
	VTKOUT << std::endl << std::endl;
	
	// data
	assert(U.n() == M->ncells() * 9); ///// !!!!!! 
	VTKOUT << "POINT_DATA " << M->ncells() * 9 << std::endl;

	for (int nc=0;nc<U.ncomp();++nc)
	  {
	    VTKOUT << "SCALARS u" << nc << " DOUBLE" << std::endl
		   << "LOOKUP_TABLE DEFAULT"  << std::endl;
	    for (int c=0;c<M->ncells();++c)
	      {
		for (int i=0;i<9;++i)
		  VTKOUT << U(9*c+i,nc) << std::endl;
	      }
	    VTKOUT << std::endl;
	  }




	// //// SIMPLE : Only outer points!
	// VTKOUT << "DATASET UNSTRUCTURED_GRID" << std::endl
	//        << "POINTS " << M->ncells()*4 << " DOUBLE" << std::endl;
	// int itoj[4]= {0,1,3,2};
	// for (int c = 0; c < M->ncells(); ++c)
	//   {
	//     const IntVector& ioc = M->IndicesOfCell(c);
	//     for (int i = 0; i < ioc.size();++i)
	//       VTKOUT << M->vertex2d(ioc[itoj[i]]) << " 0" << std::endl;
	//   }
	// VTKOUT << std::endl << std::endl;

	// // cells
	// VTKOUT << "CELLS " << M->ncells() << " " << M->ncells() * 5 << std::endl;
	// for (int c = 0; c < M->ncells(); ++c)
	//   {
	//     VTKOUT << "4 ";
	//     for (int i=0;i<4;++i) VTKOUT << c*4+i << " ";
	//     VTKOUT << std::endl;
	//   }
	// VTKOUT << std::endl << std::endl;
	
	// VTKOUT << "CELL_TYPES " << M->ncells() << std::endl;
	// for (int i = 0; i < M->ncells(); ++i)
	//   VTKOUT << "9 ";
	// VTKOUT << std::endl << std::endl;
	
	// // data
	// assert(U.n() == M->ncells() * 9); ///// !!!!!! 
	// VTKOUT << "POINT_DATA " << M->ncells() * 4 << std::endl;
	// int itojQ2[4]= {0,2,8,6};
	// for (int nc=0;nc<U.ncomp();++nc)
	//   {
	//     VTKOUT << "SCALARS u" << nc << " DOUBLE" << std::endl
	// 	   << "LOOKUP_TABLE DEFAULT"  << std::endl;
	//     for (int c=0;c<M->ncells();++c)
	//       for (int i=0;i<4;++i)
	// 	VTKOUT << U(9*c+itojQ2[i],nc) << std::endl;
	//     VTKOUT << std::endl;
	//   }
      }
    
	
	
    VTKOUT.close();
  }
}
