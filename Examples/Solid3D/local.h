
#ifndef __local_h
#define __local_h

#include "boundarysolid.h"
#include "componentinformationbase.h"
#include "dirichletdata.h"
#include "domainrighthandside.h"
#include "problemdescriptorbase.h"
#include "solid.h"
#include <stdio.h>
#include <stdlib.h>

#include "finehierarchicalmesh3d.h"
#include "hierarchicalmesh.h"
#include "hierarchicalmesh2d.h"
#include "hierarchicalmesh3d.h"
#include "localhierarchicalmesh3d.h"
#include "localmeshagent.h"
#include "stdloop.h"
#include <array>

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

template<int DIM>
class FSI_CI : public ComponentInformationBase
{
public:
  void BasicInit(const ParamFile* pf) {}

  std::string GetName() const { return "solid CI"; }

  const int GetNScalars() const { return DIM; }

  void GetScalarName(int i, std::string& s_name) const
  {
    if (DIM == 3) {
      if (i == 0)
        s_name = "vx";
      if (i == 1)
        s_name = "vy";
      if (i == 2)
        s_name = "vz";
    }
  }

  const int GetNVectors() const { return 1; }
  void GetVectorName(int i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "U_Vector";
    else
      abort();
  }
  void GetVectorIndices(int i, array<int, 3>& fa_vectorindices) const
  {
    if (i == 0) {
      fa_vectorindices[0] = 0;
      fa_vectorindices[1] = 1;
      fa_vectorindices[2] = 2;
    }
  }
};

class MyDD3d : public DirichletData
{
protected:
public:
  MyDD3d(const ParamFile* pf)
  {
    // DataFormatHandler DFH;
    // DFH.insert("vmean" ,    &vmean , 0.0);
    // FileScanner FS(DFH, pf, "Equation");
  }

  std::string GetName() const { return "MyDD"; }

  void operator()(DoubleVector& b, const Vertex3d& v, int color) const
  {
    b.zero();
  }
};

// -----------------------------------------

class ProblemDescriptor3d : public ProblemDescriptorBase
{
public:
  std::string GetName() const { return "solid"; }
  void BasicInit(const ParamFile* pf)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new Solid<3>(GetParamFile());
    GetBoundaryEquationPointer() = new BoundarySolid<3>(GetParamFile());
    GetDirichletDataPointer() = new MyDD3d(GetParamFile());

    ProblemDescriptorBase::BasicInit(pf);

    GetComponentInformationPointer() = new FSI_CI<3>;
  }
};

/* ----------------------------------------- */

class DistanceToFineMesh : public BoundaryFunction<3>
{

protected:
  int col;
  FineHierarchicalMesh3d* FHM;
  mutable double r, rr, dx;
  mutable vector<int> quads_of_closest_node;

public:
  DistanceToFineMesh(){};
  ~DistanceToFineMesh(){};

  std::string GetName() const { return "DistanceToFineMesh"; }
  void BasicInit(FineHierarchicalMesh3d* FHM_loaded, int color)
  {
    col = color;
    FHM = FHM_loaded;
  }

  void newton(Vector& dst) const
  {
    // Projektion of the node dst on surface with color col
    // surface described by boundaryelements of fine mesh FHM->nbquads()
    Vertex3d c = dst;
    Vertex3d NewPoint;

    r = 1.0e100;
    quads_of_closest_node.clear();
    // Step1: Compute Node on surface with minimal distance to the point c

    // Search of minimal distance to boundary surface
    // slope over all(!!!) boundary quads
    // Find node closest to point c
    for (int i = 0; i < FHM->nbquads(); i++) {
      if (FHM->bquad(i).material() == col) {
        for (int ii = 0; ii < 4; ii++) {
          rr = 0;
          Vertex3d _c = FHM->vertex3d(FHM->vertex_of_bquad(i, ii));
          for (int iii = 0; iii < 3; iii++) {
            dx = c[iii] - _c[iii];
            rr += dx * dx;
          }
          // Falls Knoten schon auf Interface breche ich ab
          // if(rr<tol) abort();
          if (rr > r - 1.0e-5 && rr < r + 1.0e-5)
            quads_of_closest_node.push_back(i);
          else if (rr < r - 1.0e-5) {
            quads_of_closest_node.clear();
            quads_of_closest_node.push_back(i);
            r = rr;
          }
        }
      }
    }
    r = 1.0e100;
    // find point close to the closest node on surrounding quads
    for (int i = 0; i < quads_of_closest_node.size(); i++) {

      if (FHM->bquad(quads_of_closest_node[i]).material() == col)
      // Division of quad in two triangles
      // Computation of distance of Point c to triangle
      {
        rr = 0;
        Vertex3d NewPoint_zwischen;
        rr = DistancePointTriangle(
          NewPoint_zwischen,
          c,
          FHM->vertex3d(FHM->vertex_of_bquad(quads_of_closest_node[i], 0)),
          FHM->vertex3d(FHM->vertex_of_bquad(quads_of_closest_node[i], 1)),
          FHM->vertex3d(FHM->vertex_of_bquad(quads_of_closest_node[i], 2)));
        if (rr < r) {
          r = rr;
          NewPoint = NewPoint_zwischen;
        }
        rr = DistancePointTriangle(
          NewPoint_zwischen,
          c,
          FHM->vertex3d(FHM->vertex_of_bquad(quads_of_closest_node[i], 0)),
          FHM->vertex3d(FHM->vertex_of_bquad(quads_of_closest_node[i], 2)),
          FHM->vertex3d(FHM->vertex_of_bquad(quads_of_closest_node[i], 3)));
        if (rr < r) {
          r = rr;
          NewPoint = NewPoint_zwischen;
        }
      }
    }
    dst = NewPoint;
  }
  double operator()(const Vertex3d& c) const
  {
    cout << "operator() in DistanceToFineMesh not written" << endl;
    abort();
    return 0;
  }

  // David Eberly, Geometric Tools, Redmond WA 98052
  // Copyright (c) 1998-2017
  // Distributed under the Boost Software License, Version 1.0.
  // http://www.boost.org/LICENSE_1_0.txt
  // http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
  // File Version: 3.0.0 (2016/06/19)
  double Dot(const Vertex3d& U, const Vertex3d& V) const
  {
    double dot = 0;
    for (int i = 0; i < 3; ++i) {
      dot += U[i] * V[i];
    }
    return dot;
  }

  double DistancePointTriangle(Vertex3d& new_point,
                               const Vertex3d& point,
                               const Vertex3d& triangle1,
                               const Vertex3d& triangle2,
                               const Vertex3d& triangle3) const
  {
    Vertex3d diff;
    Vertex3d edge0;
    Vertex3d edge1;
    for (int i = 0; i < 3; i++) {
      diff[i] = point[i] - triangle1[i];
      edge0[i] = triangle2[i] - triangle1[i];
      edge1[i] = triangle3[i] - triangle1[i];
    }
    double a00 = Dot(edge0, edge0);
    double a01 = Dot(edge0, edge1);
    double a11 = Dot(edge1, edge1);
    double b0 = -Dot(diff, edge0);
    double b1 = -Dot(diff, edge1);
    double const zero = 0;
    double const one = 1;
    double det = a00 * a11 - a01 * a01;
    double t0 = a01 * b1 - a11 * b0;
    double t1 = a01 * b0 - a00 * b1;

    if (t0 + t1 <= det) {
      if (t0 < zero) {
        if (t1 < zero) // region 4
        {
          if (b0 < zero) {
            t1 = zero;
            if (-b0 >= a00) // V1
            {
              t0 = one;
            } else // E01
            {
              t0 = -b0 / a00;
            }
          } else {
            t0 = zero;
            if (b1 >= zero) // V0
            {
              t1 = zero;
            } else if (-b1 >= a11) // V2
            {
              t1 = one;
            } else // E20
            {
              t1 = -b1 / a11;
            }
          }
        } else // region 3
        {
          t0 = zero;
          if (b1 >= zero) // V0
          {
            t1 = zero;
          } else if (-b1 >= a11) // V2
          {
            t1 = one;
          } else // E20
          {
            t1 = -b1 / a11;
          }
        }
      } else if (t1 < zero) // region 5
      {
        t1 = zero;
        if (b0 >= zero) // V0
        {
          t0 = zero;
        } else if (-b0 >= a00) // V1
        {
          t0 = one;
        } else // E01
        {
          t0 = -b0 / a00;
        }
      } else // region 0, interior
      {
        double invDet = one / det;
        t0 *= invDet;
        t1 *= invDet;
      }
    } else {
      double tmp0, tmp1, numer, denom;

      if (t0 < zero) // region 2
      {
        tmp0 = a01 + b0;
        tmp1 = a11 + b1;
        if (tmp1 > tmp0) {
          numer = tmp1 - tmp0;
          denom = a00 - (2) * a01 + a11;
          if (numer >= denom) // V1
          {
            t0 = one;
            t1 = zero;
          } else // E12
          {
            t0 = numer / denom;
            t1 = one - t0;
          }
        } else {
          t0 = zero;
          if (tmp1 <= zero) // V2
          {
            t1 = one;
          } else if (b1 >= zero) // V0
          {
            t1 = zero;
          } else // E20
          {
            t1 = -b1 / a11;
          }
        }
      } else if (t1 < zero) // region 6
      {
        tmp0 = a01 + b1;
        tmp1 = a00 + b0;
        if (tmp1 > tmp0) {
          numer = tmp1 - tmp0;
          denom = a00 - (2) * a01 + a11;
          if (numer >= denom) // V2
          {
            t1 = one;
            t0 = zero;
          } else // E12
          {
            t1 = numer / denom;
            t0 = one - t1;
          }
        } else {
          t1 = zero;
          if (tmp1 <= zero) // V1
          {
            t0 = one;
          } else if (b0 >= zero) // V0
          {
            t0 = zero;
          } else // E01
          {
            t0 = -b0 / a00;
          }
        }
      } else // region 1
      {
        numer = a11 + b1 - a01 - b0;
        if (numer <= zero) // V2
        {
          t0 = zero;
          t1 = one;
        } else {
          denom = a00 - (2) * a01 + a11;
          if (numer >= denom) // V1
          {
            t0 = one;
            t1 = zero;
          } else // 12
          {
            t0 = numer / denom;
            t1 = one - t0;
          }
        }
      }
    }
    for (int i = 0; i < 3; i++) {
      diff[i] = point[i] - triangle1[i] - t0 * edge0[i] - t1 * edge1[i];
      new_point[i] = triangle1[i] + t0 * edge0[i] + t1 * edge1[i];
    }
    return Dot(diff, diff);
  }
};
/* ----------------------------------------- */
#include "filescanner.h"
class ProjectionOnFineMeshAgent : public LocalMeshAgent
{
protected:
  vector<DistanceToFineMesh> DTFM;
  FineHierarchicalMesh3d* FHM;

public:
  ProjectionOnFineMeshAgent(const ParamFile* paramfile)
    : FHM(NULL)
  {
    LocalMeshAgent();

    string fine_grid_name;
    IntVector col_projec;

    DataFormatHandler DFH;
    DFH.insert("finegridname", &fine_grid_name);
    DFH.insert("colprojec", &col_projec);
    FileScanner FS(DFH, paramfile, "BoundaryManager");

    assert(FHM == NULL);
    cout << "DistanceToFineMesh --- Reading Mesh" << endl;
    FHM = new FineHierarchicalMesh3d(fine_grid_name);
    // FHM =new  FineHierarchicalMesh3d("Verzweigung_FSI_bound_neu_solid.inp");
    cout << "DistanceToFineMesh --- Finished Reading Mesh" << endl;
    // Initialisation of DistanceToFineMesh

    // Boundary Colors which should be projected

    cout << "fine_grid_name" << fine_grid_name << endl;
    DTFM.resize(col_projec.size());
    for (int i = 0; i < col_projec.size(); i++) {
      DTFM[i].BasicInit(FHM, col_projec[i]);
      AddShape(col_projec[i], &DTFM[i]);
    }
  }

  ~ProjectionOnFineMeshAgent()
  {
    if (FHM != NULL) {
      delete FHM;
      FHM = NULL;
    }
  };

  void BasicInit(const ParamFile* paramfile)
  {

    assert(HMP == NULL);
    int dim = 0;

    {
      DataFormatHandler DFH;
      DFH.insert("dimension", &dim);
      // um die zuordnung alte GMNr. -> GMNr. an/abzuschalten
      DFH.insert("cellnumtrans", &_goc2nc, false);
      FileScanner FS(DFH);
      FS.NoComplain();
      FS.readfile(paramfile, "Mesh");
    }
    {
      DataFormatHandler DFH;
      DFH.insert("periodic", &_periodicCols);
      FileScanner FS(DFH);
      FS.NoComplain();
      FS.readfile(paramfile, "BoundaryManager");
    }

    if (dim == 3) {
      HMP = new LocalHierarchicalMesh3d;
      for (map<int, BoundaryFunction<3>*>::const_iterator p = _curved3d.begin();
           p != _curved3d.end();
           p++) {
        HMP->AddShape(p->first, p->second);
      }
    } else {
      cout << "dimension of Mesh ? " << dim << endl;
    }
    assert(HMP);
    HMP->BasicInit(paramfile);

    GMG = NewMultiGridMesh();

    ReInit();
  }

  /*-----------------------------------------*/

  void BasicInit(const string& gridname,
                 int dim,
                 int patchdepth,
                 int epatcher,
                 bool goc2nc)
  {
    assert(HMP == NULL);
    _goc2nc = goc2nc;
    if (dim == 3) {
      HMP = new LocalHierarchicalMesh3d;
      for (map<int, BoundaryFunction<3>*>::const_iterator p = _curved3d.begin();
           p != _curved3d.end();
           p++) {
        HMP->AddShape(p->first, p->second);
      }
    } else {
      cout << "dimension of Mesh ? " << dim << endl;
    }
    assert(HMP);
    HMP->SetParameters(gridname, patchdepth, epatcher);

    GMG = NewMultiGridMesh();

    ReInit();
  }

  /*---------------------------------------------------*/
};

#endif
