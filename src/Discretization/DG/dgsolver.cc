#include "dgsolver.h"
#include "compose_name.h"
#include "gascoignemesh2d.h"

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
    // points // DG ... stupid? multiple times...
    VTKOUT << "DATASET UNSTRUCTURED_GRID" << std::endl
           << "POINTS " << M->ncells()*4 << " DOUBLE" << std::endl;
    int itoj[4]= {0,1,3,2};
    for (int c = 0; c < M->ncells(); ++c)
      {
	const IntVector& ioc = M->IndicesOfCell(c);
	for (int i = 0; i < ioc.size();++i)
	  {
	    VTKOUT << M->vertex2d(ioc[itoj[i]]) << " 0" << std::endl;
	  }
	
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
    const GlobalVector& U = GetGV(gu);
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
    
    
    

    VTKOUT.close();
  }
}
