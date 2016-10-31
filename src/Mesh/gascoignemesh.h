/**
 *
 * Copyright (C) 2004, 2005, 2006, 2007, 2011 by the Gascoigne 3D authors
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


#ifndef  __GascoigneMesh_h
#define  __GascoigneMesh_h

#include  "patchmesh.h"
#include  "patchindexhandler.h"
#include  "boundaryindexhandler.h"
#include  "hangingindexhandler.h"
#include  "gascoigne.h"

/*-----------------------------------------*/

namespace Gascoigne
{

  class GascoigneMeshBase : public PatchMesh
  {
  protected:

    IntVector              nc, vo2n;
    PatchIndexHandler      PatchHandler;
    BoundaryIndexHandler   BoundaryHandler;
    HangingIndexHandler    HangingHandler;
  

  public:

    GascoigneMeshBase()          {};
    virtual ~GascoigneMeshBase() {};

    void BasicInit(const ParamFile* pf) {
      std::cerr << "\"GascoigneMesh:BasicInit\" not written!" << std::endl;
      abort();
    }

    const IntVector&               GetCellVector()          const  { return nc;}
    const PatchIndexHandler&    GetPatchIndexHandler()   const  { return PatchHandler;}
    const BoundaryIndexHandler&    GetBoundaryIndexHandler()const  { return BoundaryHandler;}
    const HangingIndexHandler&  GetHangingIndexHandler() const  { return HangingHandler;}
    const IntVector*               Vertexo2n()              const  { return &vo2n;}

    IntVector&               GetCellVector()           { return nc;}
    PatchIndexHandler&    GetPatchIndexHandler()    { return PatchHandler;}
    BoundaryIndexHandler&    GetBoundaryIndexHandler() { return BoundaryHandler;}
    HangingIndexHandler&  GetHangingIndexHandler()  { return HangingHandler;}
    IntVector*               Vertexo2n()               { return &vo2n;}

    int  nodes_per_patch()   const { return PatchHandler.nodes_per_patch();}
    int  nodes_per_q4patch() const { return PatchHandler.nodes_per_q4patch();}
    bool HasPatch()          const { return PatchHandler.HasPatch();}
    bool HasQ4Patch()        const { return PatchHandler.HasQ4Patch();}
    int  npatches()          const { return PatchHandler.npatches(); }
    int  nq4patches()        const { return PatchHandler.nq4patches(); }

    const IntVector* IndicesOfPatch    (int i)     const { return &(PatchHandler.IndicesOfPatch(i));}
    const IntVector* IndicesOfQ4Patch  (int i)     const { return &(PatchHandler.IndicesOfQ4Patch(i));}
    const IntVector* VertexOnBoundary(int color) const { return &(BoundaryHandler.Verteces(color)); }
    const IntVector* CellOnBoundary  (int color) const { return &(BoundaryHandler.Cells(color));    }
    const IntVector* LocalOnBoundary (int color) const { return &(BoundaryHandler.Localind(color)); }
    const IntVector* PatchOnBoundary  (int color) const { return &(BoundaryHandler.Patches(color));    }
    const IntVector* LocalPatchOnBoundary (int color) const { return &(BoundaryHandler.LocalPatchind(color)); }
  
    std::set<int> GetColors()             const { return BoundaryHandler.GetColors();}
    IntVector  Q2IndicesOfQ4Patch(int i)  const { return PatchHandler.Q2IndicesOfQ4Patch(i);}
    IntVector  CoarseIndices(int iq)   const { return PatchHandler.CoarseIndices(iq);}
    IntVector  CoarseIndicesQ4(int iq) const { return PatchHandler.CoarseIndicesQ4(iq);}

    virtual IntVector  IndicesOfCell(int iq) const{ 
      std::cerr << "\"GascoigneMesh:IndicesOfCell\" not written!" << std::endl;
      abort();
    }
    
    virtual bool CellIsCurved(int iq) const { return 0;}
    virtual int  nedges()    const { return 0; }
  };




  template<int DIM>
    class GascoigneMesh : public GascoigneMeshBase
    {
    protected:
      
      // basic
      std::vector<Vertex<DIM> >   nx;
      
    public:
  
      GascoigneMesh<DIM>()  {};
      ~GascoigneMesh<DIM>() {};
      
      std::string GetName() const {return "GascoigneMesh-DIM";}
      
      std::vector<Vertex<DIM> >& GetVertexVector()       {return nx;}
      const std::vector<Vertex<DIM> >& GetVertexVector() const {return nx;}
      
      int  dimension() const {return DIM; }
      int  nnodes()    const {return nx.size();}
      int  ncells()    const {return nc.size()/4;}
      int  nhanging()  const { return HangingHandler.GetStructure()->size(); }
      
      int  nodes_per_cell(int i)  const { return (DIM==2)?4:8;  }
      int  VtkType(int i)         const { return (DIM==2)?9:12; }
      
      const Vertex<DIM>& vertex(int i) const { return nx[i];}
      int  vertex_of_cell(int i, int ii) const { return nc[( (DIM==2)?4:8 )*i+ii]; }
      
      IntVector  IndicesOfCell(int iq) const;
    };

#define GascoigneMesh2d GascoigneMesh<2>
#define GascoigneMesh3d GascoigneMesh<3>
  
}




#endif
