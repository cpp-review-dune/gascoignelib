/*----------------------------   cgdofhandler.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __cgdofhandler_H
#define __cgdofhandler_H
/*----------------------------   cgdofhandler.h     ---------------------------*/


/**
 *
 * Copyright (C) 2019, 2020 by the Gascoigne 3D authors
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


#include "dofhandler.h"

/*-----------------------------------------*/

namespace Gascoigne
{
  /**
   * Handler for degrees of freedom for discretizations
   * with globally continuous functions, e.g. standard
   * lagrange discs. 
   * 
   * The DofHandler collects elements, each element has M^DIM
   * dofs, i.e. in 2D:
   *
   * M=2
   *
   * 2 - 3
   * |   |
   * 0 - 1
   *
   * M=3
   *
   * 6 - 7 - 8
   * |       |
   * 3   4   5
   * |       |
   * 0 - 1 - 2
   *
   * and general
   *
   * M(M-1) ..    M^2-1
   * .                .
   * .                .
   * M - ...      (2M-1)
   * 0 - 1 - ... - (M-1)
   *
   * In 3D, sorting in x-y-z direction
   *
   *   6----7
   *  /|   /|    y  z
   * 2----3 |    | /
   * | 4--|-5    |/
   * |/   |/     +---x
   * 0 -- 1
   *
   *
   * HANGING NODES
   * 
   * 2:1 mesh refinement is allowed.
   *
   *
   ***/
   
  /// DIM is spatial dimension, M dof's per element in each direction
  template <int DIM, int M>
  class CGDofHandler : public DofHandlerBase
  {
  protected:
    /// geometric coordinates for all dof's
    std::vector<Vertex<DIM> > nx;

  public:
    CGDofHandler()  {};
    ~CGDofHandler() {};

    /**
     * Initialize from mesh-structure in Base DofHandler<DIM> = GascoigneMesh<DIM>
     *
     * Each cell of the GascoigneMesh is filled with M^DIM dof's
     * dof's on nodes / edges / faces are shared
     **/
    void InitFromGascoigneMesh(const DofHandler<DIM>& GM);

    /// General Access
    int dimension() const        { return DIM; }
    int dofs_per_element() const { return pow(M,DIM); }
    int ndofs() const            { return nx.size(); }
    int nelements() const        { return nc.size() / dofs_per_element(); }
    int nhanging() const         { return HangingHandler.GetStructure()->size(); }

    IntVector GetElement(int iq) const
    {
      abort();
      
      IntVector tmp;
      int start = dofs_per_element()*iq;
      for (int i=0;i<dofs_per_element();++i)
	tmp.push_back(nc[start+i]);
      return tmp;
    }

    
    /// Access to coordinates    
    std::vector<Vertex<DIM> >  &GetVertexVector()            { return nx; }
    const std::vector<Vertex<DIM> > &GetVertexVector() const { return nx; }

    const Vertex<DIM> &vertex(int i) const { return nx[i]; }
    virtual const Vertex<2> &vertex2d(int i) const {abort();}
    virtual const Vertex<3> &vertex3d(int i) const {abort();}
    
    ////// Boundary
    const IntVector *ElementOnBoundary(int color) const
    {
      return &(BoundaryHandler.Cells(color));
    }
    const IntVector *ElementLocalOnBoundary(int color) const
    {
      return &(BoundaryHandler.Localind(color));
    }
    const IntVector *ElementOnBoundary(int degree, int color) const
    {
      std::cerr << "Element on Boundary with degree not used!" << std::endl;
      abort();
    }
    const IntVector *ElementLocalOnBoundary(int degree,int color) const
    {
      std::cerr << "ElementLocal on Boundary with degree not used!" << std::endl;
      abort(); 
    }

    int VtkType(int i) const
    {
      return (DIM==2)?9:12;
    }
        
    int vertex_of_cell(int i, int ii) const
    {
      return nc[dofs_per_element() * i + ii];
    }
    
    /// Dummy? Old Interface
    int ncells()              const { abort(); }
    int nelements(int degree) const { abort(); }
    int nodes_per_cell(int i) const { abort(); } 
    int nnodes()              const { abort(); }
    
    IntVector IndicesOfCell(int iq) const
    { std::cerr << "CGDofHandler: Use GetElement" << std::endl; assert(0); }

    IntVector GetElement(int degree, int iq) const
    { std::cerr << "GetElement with degree not used" << std::endl; abort(); }

  };

  template<> inline const Vertex<2>& CGDofHandler<2,2>::vertex2d(int i) const { return vertex(i); }
  template<> inline const Vertex<2>& CGDofHandler<2,3>::vertex2d(int i) const { return vertex(i); }
  template<> inline const Vertex<2>& CGDofHandler<2,5>::vertex2d(int i) const { return vertex(i); }
  template<> inline const Vertex<3>& CGDofHandler<3,2>::vertex3d(int i) const { return vertex(i); }
  template<> inline const Vertex<3>& CGDofHandler<3,3>::vertex3d(int i) const { return vertex(i); }
  template<> inline const Vertex<3>& CGDofHandler<3,5>::vertex3d(int i) const { return vertex(i); }
  
} // namespace Gascoigne


/*----------------------------   dofhandler.h     ---------------------------*/
/* end of #ifndef __dofhandler_H */
#endif
/*----------------------------   dofhandler.h     ---------------------------*/
