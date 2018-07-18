/**
*
* Copyright (C) 2004, 2005, 2006, 2007 by the Gascoigne 3D authors
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


#ifndef __finehierarchicalmesh3d_h
#define __finehierarchicalmesh3d_h

#include  "vertex.h" 
#include  "boundaryquad.h" 
#include  "hex.h" 
#include  "hexlawandorder.h" 
#include  "boundaryfunction.h"
#include  "hangcontainer3d.h" 
#include  "hierarchicalmesh3d.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class FineHierarchicalMesh3d : public HierarchicalMesh3d
{
  protected :
  
  /*  typedef  */

  typedef  std::vector<Vertex3d>       VertexVec3d;
  typedef  BoundaryCell<4>          BoundaryQuad;
  typedef  std::vector<Hex>            HexVec;
  typedef  std::vector<BoundaryQuad>   BQuadVec;
  typedef  BoundaryFunction<3>      BoundaryFunction3d;


  /*  Data  */
  VertexVec3d        vertexs3d; 

  /* info fuer interpolation auf neues gitter */
  HexVec             hexs;
  BQuadVec           Bquads;
  HexLawAndOrder     HexLaO;


  /*  Functionen  */

  


  std::pair<bool,tint>  check_inp(const std::string&);

 
  public:


  //FineHierarchicalMesh3d(const FineHierarchicalMesh3d& H);
  //FineHierarchicalMesh3d& operator=(const FineHierarchicalMesh3d& H);
  FineHierarchicalMesh3d(const std::string gridname);
  ~FineHierarchicalMesh3d() {  }

  std::string GetName() const {return "FineHierarchicalMesh3d";}

  /*  Zugriff  */

  int  dimension()            const { return 3;}

  int  nnodes   ()            const { return vertexs3d.size();}
  int  ncells   ()            const { return hexs.size();}
  int  nbquads  ()            const { return Bquads.size();}

  int  nodes_per_cell(int i)  const { return 8;}
  int  VtkType(int i) const { return 12;}



  const VertexVec3d& GetVertexVector() const {return vertexs3d; }
  VertexVec3d& GetVertexVector() {return vertexs3d; }

  const Vertex3d& vertex3d(int i)         const { return vertexs3d[i];}

  const Hex&   hex (int i)                const { return hexs[i];}
  const BoundaryQuad&  bquad(int i)       const { return Bquads[i];}

  int  vertex_of_cell (int i, int ii)      const { return hexs[i].vertex(ii);}
  int  vertex_of_bquad(int i, int ii)      const { return Bquads[i].vertex(ii);}
  int  face_of_hex    (int i, int ii)      const { return hexs[i].edge(ii); }
  int  level(int i)                        const { return hexs[i].level();}
  bool sleep(int i)                        const { return hexs[i].sleep();}

  const std::vector<BoundaryQuad>& quad_list() const { return Bquads; }

  const VertexVec3d&        vertex3d() const {return vertexs3d;}
  const HexVec&             hex     () const {return hexs;}
  const BQuadVec&           bquad   () const {return Bquads;}



  /*  Functionen  */

  void   read_inp (const std::string&);


};
}

/*---------------------------------------------------*/

#endif
