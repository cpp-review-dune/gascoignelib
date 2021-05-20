
#ifndef __local_h
#define __local_h

#include "boundaryfsi.h"
#include "componentinformationbase.h"
#include "dirichletdata.h"
#include "domainrighthandside.h"
#include "fsi.h"
#include "problemdescriptorbase.h"
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

extern double __TIME;

template<int DIM>
class FSI_CI : public ComponentInformationBase
{
public:
  void BasicInit(const ParamFile* pf) {}

  std::string GetName() const { return "FSI CI"; }

  const int GetNScalars() const { return 2 * DIM + 1; }

  void GetScalarName(int i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "p";
    if (DIM == 2) {
      if (i == 1)
        s_name = "vx";
      if (i == 2)
        s_name = "vy";
      if (i == 3)
        s_name = "ux";
      if (i == 4)
        s_name = "uy";
    }
    if (DIM == 3) {
      if (i == 1)
        s_name = "vx";
      if (i == 2)
        s_name = "vy";
      if (i == 3)
        s_name = "vz";
      if (i == 4)
        s_name = "ux";
      if (i == 5)
        s_name = "uy";
      if (i == 6)
        s_name = "uz";
    }
  }

  const int GetNVectors() const { return 2; }
  void GetVectorName(int i, std::string& s_name) const
  {
    if (i == 0)
      s_name = "V";
    else if (i == 1)
      s_name = "U";
    else
      abort();
  }
  void GetVectorIndices(int i, array<int, 3>& fa_vectorindices) const
  {
    if (i == 0) {
      fa_vectorindices[0] = 1;
      fa_vectorindices[1] = 2;
      fa_vectorindices[2] = -1;
      if (DIM == 3)
        fa_vectorindices[3] = 3;
    } else {
      fa_vectorindices[0] = DIM + 1;
      fa_vectorindices[1] = DIM + 2;
      fa_vectorindices[2] = -1;
      if (DIM == 3)
        fa_vectorindices[3] = DIM + 3;
    }
  }
};

class MyDD : public DirichletData
{
protected:
  double vmean;

public:
  MyDD(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("vmean", &vmean, 0.0);
    FileScanner FS(DFH, pf, "Equation");
  }

  std::string GetName() const { return "MyDD"; }

  void operator()(DoubleVector& b, const Vertex2d& v, int color) const
  {
    b.zero();

    double t = __TIME;
    double sc = 1.0;
    if (t < 1.0)
      sc = 0.5 - 0.5 * cos(M_PI * t);

    double veff = vmean * sc;

    if (color == 0)
      b[1] += v.y() * (0.41 - v.y()) / 0.205 / 0.205 * veff * 1.5;

    //    b.zero();
    /*
      double t= __TIME;


      double x = v.x();  double y = v.y(); double z = v.z();
      double fact_time=0;
      if (__TIME >= 0.0 && __TIME < 0.2)
      {

      fact_time =  (-11292 * __TIME*__TIME+ 2493.4 * __TIME);

      //           b[1] =  (-9291.6667 * __TIME * __TIME
      //                + 1893.3333 * __TIME
      //                + 40)
      //             *__EINSTROMGESCHWINDIGKEIT* v.y() * (1.61-v.y())
      *4/1.61/1.61;
      }
      else if (__TIME >= 0.9 && __TIME < 1.1)
      {
      fact_time =  (-9291.6667 * (__TIME-0.9) * (__TIME-0.9)
      + 1893.3333 * (__TIME-0.9)
      + 40);
      }
      else if (__TIME >= 1.8 && __TIME < 2.0)
      {
      fact_time=  (-9291.6667 * (__TIME-1.8) * (__TIME-1.8)
      + 1893.3333 * (__TIME-1.8)
      + 40);
      }
      else if (__TIME >= 2.7 && __TIME < 2.9)
      {
      fact_time =  (-9291.6667 * (__TIME-2.7) * (__TIME-2.7)
      + 1893.3333 * (__TIME-2.7)
      + 40);
      }
      else if (__TIME >= 0.2 && __TIME < 0.25)
      {
      fact_time = (120 * (__TIME - 0.2) + 47);
      }
      else if (__TIME >= 1.1 && __TIME < 1.15)
      {
      fact_time =  (120 * (__TIME - 1.1) + 47);
      }
      else if (__TIME >= 2.0 && __TIME < 2.05)
      {
      fact_time=  (120 * (__TIME - 2.0) + 47);
      }
      else if (__TIME >= 2.9 && __TIME < 2.95)
      {
      fact_time =  (120 * (__TIME - 2.9) + 47);
      }
      else if (__TIME >= 0.25 && __TIME < 0.3)
      {
      fact_time = (-220 * (__TIME - 0.25) + 53);
      }
      else if (__TIME >= 1.15 && __TIME < 1.2)
      {
      fact_time = (-220 * (__TIME - 1.15) + 53);
      }
      else if (__TIME >= 2.05 && __TIME < 2.1)
      {
      fact_time =  (-220 * (__TIME - 2.05) + 53);
      }
      else if (__TIME >= 2.95 && __TIME < 3.0)
      {
      fact_time = (-220 * (__TIME - 2.95) + 53);
      }
      else if (__TIME >= 0.3 && __TIME < 0.37)
      {
      fact_time = (114.2857142 * (__TIME - 0.3) + 42);
      }
      else if (__TIME >= 1.2 && __TIME < 1.27)
      {
      fact_time = (114.2857142 * (__TIME - 1.2) + 42);
      }
      else if (__TIME >= 2.1 && __TIME < 2.17)
      {
      fact_time = (114.2857142 * (__TIME - 2.1) + 42);
      }
      else if (__TIME >= 3.0 && __TIME < 3.07)
      {
      fact_time = (114.2857142 * (__TIME - 3.0) + 42);
      }
      else if (__TIME >= 0.37 && __TIME < 0.9)
      {
      fact_time = (50*exp(-0.4210255685 * (__TIME - 0.37)));
      }
      else if (__TIME >= 1.27 && __TIME < 1.8)
      {
      fact_time = (50*exp(-0.4210255685 * (__TIME - 1.27)));
      }
      else if (__TIME >= 2.17 && __TIME < 2.7)
      {
      fact_time = (50*exp(-0.4210255685 * (__TIME - 2.17)));
      }
      else if (__TIME >= 3.07 )
      {
      fact_time = (50*exp(-0.4210255685 * (__TIME - 3.07)));
      }



      if (color==7)
      {
      Vertex3d center; center[0]=3.78223304703; center[1]=-3.882233047;
      center[2]=1.0928932; Vertex3d normal;  normal[0]=  0.5000  ; normal[1]=
      -0.4998 ; normal[2]=   0.7072; double diam=0.3; b[1] =
      -normal[0]*fact_time * (
      (x-center[0])*(x-center[0])+(y-center[1])*(y-center[1])+(z-center[2])*(z-center[2])-diam*diam)/(-diam*diam);
      b[2] = -normal[1] *fact_time *(
      (x-center[0])*(x-center[0])+(y-center[1])*(y-center[1])+(z-center[2])*(z-center[2])-diam*diam)/(-diam*diam);
      b[3] = -normal[2]*fact_time * (
      (x-center[0])*(x-center[0])+(y-center[1])*(y-center[1])+(z-center[2])*(z-center[2])-diam*diam)/(-diam*diam);
      }
    */

    // double teff = __TIME;
    // int Ti=0;
    // while (teff>5.0)
    //   {
    // 	++Ti;
    // 	teff -=5.0;
    //   }

    // // transition by 0.1
    // double vnew = 1.0+0.1*Ti;
    // double vold = 1.0+0.1*(Ti-1.0);
    // if (vold<1.0) vold = 1.0;
    // double sc = 1.0;
    // if (teff<0.1)
    //   sc = 0.5*(1.0-cos(M_PI*teff/0.1));

    // //    sc = 1.0;

    // //    ////// Std-Werte v=2
    // vnew = vold = 1.15; // re = 130

    // double veff = vold + sc*(vnew-vold);
    // //    cout << veff << endl;

    // if (color==0)
    //   b[1] += v.y()*(0.41-v.y())/0.205/0.205 * veff * 1.5;

    // return;

    // //    double sc = 1.0;
    // // transition from zero to vmean at time t=0 (2 sec)
    // if (__TIME<2.0)
    //   {
    // 	sc = 0.5*(1.0-cos(M_PI*__TIME/2.0));
    // 	if (color==0)
    // 	  b[1] += v.y()*(0.41-v.y())/0.205/0.205 * vmean * 1.5 * sc;
    //   }
    // else if (__TIME<10.0)
    //   {
    // 	sc = 1.0;
    // 	if (color==0)
    // 	  b[1] += v.y()*(0.41-v.y())/0.205/0.205 * vmean * 1.5 * sc;
    //   }
    // else if (__TIME<11.0) // transition from 2.0 to vmean in t in [10,11]
    //   {
    // 	sc = 2.0 + 0.5*(1.0-cos(M_PI*(__TIME-10.0)))*(vmean-2.0);

    // 	if (color==0)
    // 	  b[1] += v.y()*(0.41-v.y())/0.205/0.205 *  1.5 * sc;

    //   }
    // else
    //   {
    // 	sc = 1.0;
    // 	if (color==0)
    // 	  b[1] += v.y()*(0.41-v.y())/0.205/0.205 * vmean * 1.5 * sc;
    //   }

    // // transition from 2.0 to vmean at time t=10 (1 sec)
    // if (__TIME>10.0)
    //   {
    // 	sc = 0.5*(1.0-cos(M_PI*__TIME/2.0));

    //   }
  }
};

class MyDD3d : public DirichletData
{
protected:
  double vmean;

public:
  MyDD3d(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("vmean", &vmean, 0.0);
    FileScanner FS(DFH, pf, "Equation");
  }

  std::string GetName() const { return "MyDD"; }

  void operator()(DoubleVector& b, const Vertex3d& v, int color) const
  {
    b.zero();
    /*
      double sc = 1.0;
      if (__TIME<2.0) sc = 0.5*(1.0-cos(M_PI*__TIME/2.0));

      if (color==0)
      b[1] += v.y()*(0.4-v.y())/0.2/0.2 * (0.4-v.z())*(0.4+v.z())/0.4/0.4 *
      vmean * sc * 9.0/8.0;

    */

    double x = v.x();
    double y = v.y();
    double z = v.z();
    double fact_time = 0;
    if (__TIME >= 0.0 && __TIME < 0.2) {

      fact_time = (-11292 * __TIME * __TIME + 2493.4 * __TIME);

      //           b[1] =  (-9291.6667 * __TIME * __TIME
      //                + 1893.3333 * __TIME
      //                + 40)
      //             *__EINSTROMGESCHWINDIGKEIT* v.y() * (1.61-v.y())
      //             *4/1.61/1.61;
    } else if (__TIME >= 0.9 && __TIME < 1.1) {
      fact_time = (-9291.6667 * (__TIME - 0.9) * (__TIME - 0.9) +
                   1893.3333 * (__TIME - 0.9) + 40);
    } else if (__TIME >= 1.8 && __TIME < 2.0) {
      fact_time = (-9291.6667 * (__TIME - 1.8) * (__TIME - 1.8) +
                   1893.3333 * (__TIME - 1.8) + 40);
    } else if (__TIME >= 2.7 && __TIME < 2.9) {
      fact_time = (-9291.6667 * (__TIME - 2.7) * (__TIME - 2.7) +
                   1893.3333 * (__TIME - 2.7) + 40);
    } else if (__TIME >= 0.2 && __TIME < 0.25) {
      fact_time = (120 * (__TIME - 0.2) + 47);
    } else if (__TIME >= 1.1 && __TIME < 1.15) {
      fact_time = (120 * (__TIME - 1.1) + 47);
    } else if (__TIME >= 2.0 && __TIME < 2.05) {
      fact_time = (120 * (__TIME - 2.0) + 47);
    } else if (__TIME >= 2.9 && __TIME < 2.95) {
      fact_time = (120 * (__TIME - 2.9) + 47);
    } else if (__TIME >= 0.25 && __TIME < 0.3) {
      fact_time = (-220 * (__TIME - 0.25) + 53);
    } else if (__TIME >= 1.15 && __TIME < 1.2) {
      fact_time = (-220 * (__TIME - 1.15) + 53);
    } else if (__TIME >= 2.05 && __TIME < 2.1) {
      fact_time = (-220 * (__TIME - 2.05) + 53);
    } else if (__TIME >= 2.95 && __TIME < 3.0) {
      fact_time = (-220 * (__TIME - 2.95) + 53);
    } else if (__TIME >= 0.3 && __TIME < 0.37) {
      fact_time = (114.2857142 * (__TIME - 0.3) + 42);
    } else if (__TIME >= 1.2 && __TIME < 1.27) {
      fact_time = (114.2857142 * (__TIME - 1.2) + 42);
    } else if (__TIME >= 2.1 && __TIME < 2.17) {
      fact_time = (114.2857142 * (__TIME - 2.1) + 42);
    } else if (__TIME >= 3.0 && __TIME < 3.07) {
      fact_time = (114.2857142 * (__TIME - 3.0) + 42);
    } else if (__TIME >= 0.37 && __TIME < 0.9) {
      fact_time = (50 * exp(-0.4210255685 * (__TIME - 0.37)));
    } else if (__TIME >= 1.27 && __TIME < 1.8) {
      fact_time = (50 * exp(-0.4210255685 * (__TIME - 1.27)));
    } else if (__TIME >= 2.17 && __TIME < 2.7) {
      fact_time = (50 * exp(-0.4210255685 * (__TIME - 2.17)));
    } else if (__TIME >= 3.07) {
      fact_time = (50 * exp(-0.4210255685 * (__TIME - 3.07)));
    }

    if (color == 7) {
      Vertex3d center;
      center[0] = 3.78223304703;
      center[1] = -3.882233047;
      center[2] = 1.0928932;
      Vertex3d normal;
      normal[0] = 0.5000;
      normal[1] = -0.4998;
      normal[2] = 0.7072;
      double diam = 0.3;
      b[1] +=
        -normal[0] * 1.0e-1 * fact_time *
        ((x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) +
         (z - center[2]) * (z - center[2]) - diam * diam) /
        (-diam * diam);
      b[2] +=
        -normal[1] * 1.0e-1 * fact_time *
        ((x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) +
         (z - center[2]) * (z - center[2]) - diam * diam) /
        (-diam * diam);
      b[3] +=
        -normal[2] * 1.0e-1 * fact_time *
        ((x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) +
         (z - center[2]) * (z - center[2]) - diam * diam) /
        (-diam * diam);
    }
  }
};

// -----------------------------------------

class ProblemDescriptor2d : public ProblemDescriptorBase
{
public:
  std::string GetName() const { return "fsi"; }
  void BasicInit(const ParamFile* pf)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new FSI<2>(GetParamFile());
    //    GetBoundaryEquationPointer() = new FSI<2>(GetParamFile());
    GetDirichletDataPointer() = new MyDD(GetParamFile());

    ProblemDescriptorBase::BasicInit(pf);

    //    GetComponentInformationPointer() = new FSI_CI<2>;
  }
};

class ProblemDescriptor3d : public ProblemDescriptorBase
{
public:
  std::string GetName() const { return "fsi"; }
  void BasicInit(const ParamFile* pf)
  {
    GetParamFilePointer() = pf;
    GetEquationPointer() = new FSI<3>(GetParamFile());
    GetBoundaryEquationPointer() = new BoundaryFSI(GetParamFile());
    GetDirichletDataPointer() = new MyDD3d(GetParamFile());

    ProblemDescriptorBase::BasicInit(pf);

    //    GetComponentInformationPointer() = new FSI_CI<3>;
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

class ProjectionOnFineMeshAgent : public LocalMeshAgent
{
protected:
  vector<DistanceToFineMesh> DTFM;
  FineHierarchicalMesh3d* FHM;

public:
  ProjectionOnFineMeshAgent()
    : FHM(NULL)
  {
    LocalMeshAgent();

    assert(FHM == NULL);
    cout << "DistanceToFineMesh --- Reading Mesh" << endl;
    FHM = new FineHierarchicalMesh3d("Verzweigung_FSI_fine_bound_neu.inp");
    cout << "DistanceToFineMesh --- Finished Reading Mesh" << endl;
    // Initialisation of DistanceToFineMesh
    DTFM.resize(8);
    DTFM[0].BasicInit(FHM, 1);
    DTFM[1].BasicInit(FHM, 5);
    DTFM[2].BasicInit(FHM, 3);
    DTFM[3].BasicInit(FHM, 6);
    DTFM[4].BasicInit(FHM, 8);
    DTFM[5].BasicInit(FHM, 9);
    DTFM[6].BasicInit(FHM, 10);
    DTFM[7].BasicInit(FHM, 14);
    // Boundary Colors which should be projected
    AddShape(1, &DTFM[0]);
    AddShape(5, &DTFM[1]);
    AddShape(3, &DTFM[2]);
    AddShape(6, &DTFM[3]);
    AddShape(8, &DTFM[4]);
    AddShape(9, &DTFM[5]);
    AddShape(10, &DTFM[6]);
    AddShape(14, &DTFM[7]);
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

    if (dim == 2) {
      HMP = new HierarchicalMesh2d;
      for (map<int, BoundaryFunction<2>*>::const_iterator p = _curved2d.begin();
           p != _curved2d.end();
           p++) {
        HMP->AddShape(p->first, p->second);
      }
    } else if (dim == 3) {
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
    if (dim == 2) {
      HMP = new HierarchicalMesh2d;
      for (map<int, BoundaryFunction<2>*>::const_iterator p = _curved2d.begin();
           p != _curved2d.end();
           p++) {
        HMP->AddShape(p->first, p->second);
      }
    } else if (dim == 3) {
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
  // void inner_vertex_newton3d(const IntVector& vnew,
  // 			     const IntSet& CellRefList,
  // 			     const IntSet& adjustvertex)
  // {
  //   if (GetCurvedShapes().empty()) return;

  //   // baue lokalen set auf um spaeter nicht alle hex zu justieren
  //   //Urspruengliche intention nur Randzellen anzupassen. Will ich hier ja
  //   gerade nicht!!!
  //   //IntSet Hexset;
  //   //for (int i=0; i<Bquads.size(); i++)
  //   //  {
  //   //    int hi = Bquads[i].of_quad();
  //   //    Hexset.insert(hi);
  //   //  }

  //   std::array<int,4> v;
  //   IntSetIt  cp=CellRefList.begin();

  //   //Schleife uber alle Zellen, zaehlen der Nachbarzellen eines
  //   Kantenmittelpunktes
  //   // Im Gebiet normalerweise 4 am Rand 2 (weitere Werte moegloch!)
  //   vector<int> RealBoundaryFaces;
  //   RealBoundaryFaces.clear();
  //   RealBoundaryFaces.resize(vertexs3d.size(),0);

  //   for (int i=0; i<CellRefList.size(); i++)
  //     {
  // 	int hi = co2n[*cp++];
  // 	const Hex& h = hex(hi);
  // 	for (int e=0;e<12;++e)
  // 	  {
  // 	    int ev = HexLaO.edge_vertex(h,e);
  // 	    RealBoundaryFaces[ev]++;
  // 	  }
  //     }

  //   //for (int i=0; i<Bquads.size(); i++)
  //   //  {
  //   //    int hi = Bquads[i].of_quad();
  //   //    Hexset.insert(hi);
  //   //  }
  //   //////////////////////////////////////////////////////////
  //   vector< vector<int> > Knotengegenueber;
  //   Knotengegenueber.resize(12);
  //   Knotengegenueber[0].push_back(2); Knotengegenueber[0].push_back(4);
  //   Knotengegenueber[1].push_back(3); Knotengegenueber[1].push_back(5);
  //   Knotengegenueber[2].push_back(0); Knotengegenueber[2].push_back(6);
  //   Knotengegenueber[3].push_back(7); Knotengegenueber[3].push_back(1);
  //   Knotengegenueber[4].push_back(0); Knotengegenueber[4].push_back(6);
  //   Knotengegenueber[5].push_back(1); Knotengegenueber[5].push_back(7);
  //   Knotengegenueber[6].push_back(4); Knotengegenueber[6].push_back(2);
  //   Knotengegenueber[7].push_back(3); Knotengegenueber[7].push_back(5);
  //   Knotengegenueber[8].push_back(9); Knotengegenueber[8].push_back(11);
  //   Knotengegenueber[9].push_back(8); Knotengegenueber[9].push_back(10);
  //   Knotengegenueber[10].push_back(9); Knotengegenueber[10].push_back(11);
  //   Knotengegenueber[11].push_back(8); Knotengegenueber[11].push_back(10);

  //   IntSet adjustadditionalvertex;
  //   IntSet localadjustvertex=adjustvertex;
  //   /////////////////////////////////////////////////////////
  //   int iter=0;
  //   while(adjustadditionalvertex.size()!=0 ||iter==0)
  //     {
  // 	iter++;
  // 	adjustadditionalvertex.clear();
  // 	cp=CellRefList.begin();
  // 	for (int i=0; i<CellRefList.size(); i++)
  // 	  {
  // 	    int hi = co2n[*cp++];
  // 	    const Hex& h = hex(hi);

  // 	    //if (Hexset.find(hi)==Hexset.end()) continue;
  // 	    //Abfrage ob es sich um eine Randzelle handelt

  // 	    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEU
  // 	    // edges
  // 	    for (int e=0;e<12;++e)
  // 	      {
  // 		int ev = HexLaO.edge_vertex(h,e);
  // 		if (localadjustvertex.find(ev)!=localadjustvertex.end())
  // continue;

  // 		if
  // (adjustadditionalvertex.find(ev)==adjustadditionalvertex.end())
  // 		  {
  // 		    std::array<int,2> fe;
  // 		    HexLaO.globalvertices_of_edge(h,fe,e);
  // 		    vertexs3d[ev]  = vertexs3d[fe[0]];
  // 		    vertexs3d[ev]  += vertexs3d[fe[1]];
  // 		    vertexs3d[ev]	 *=0.5;
  // 		  }
  // 		//uberpruefen ob die Gegeneuberliegendende Kanten Knoten
  // verschoben worden sind 		for(int ii=0;ii<2;ii++)
  // 		  {
  // 		    int ev_gegenueber =
  // HexLaO.edge_vertex(h,Knotengegenueber[e][ii]); 		    if
  // (localadjustvertex.find(ev_gegenueber)!=localadjustvertex.end())
  // 		      {
  // 			//berechnen wie weit der Knoten ev_gegenueber verschoben
  // worden ist 			std::array<int,2> fe;
  // 			HexLaO.globalvertices_of_edge(h,fe,Knotengegenueber[e][ii]);
  // 			Vertex3d verschiebungsvector=vertexs3d[ev_gegenueber];
  // 			verschiebungsvector*=2.0;
  // 			verschiebungsvector -= vertexs3d[fe[0]];
  // 			verschiebungsvector -= vertexs3d[fe[1]];
  // 			verschiebungsvector*=0.5;
  // 			//Verhaltnis verschiebungsvector zu Zellhoehe:
  // 			double lengthverschiebungsvector=0;
  // 			double lengthcell =0;
  // 			for(int iii=0;iii<3;iii++)
  // 			  {
  // 			    lengthverschiebungsvector+=verschiebungsvector[iii]*verschiebungsvector[iii];
  // 			    lengthcell+=(0.5*vertexs3d[fe[0]][iii]+0.5*vertexs3d[fe[1]][iii]-vertexs3d[ev][iii])*(0.5*vertexs3d[fe[0]][iii]+0.5*vertexs3d[fe[1]][iii]-vertexs3d[ev][iii]);
  // 			  }
  // 			//addieren von 0.25* verschiebungsvector auf den Knoten
  // ev(knoten grenzt im normalvall an 4 patches!!
  // 			// bei Randknoten 0.5
  // 			if(lengthverschiebungsvector>lengthcell)
  // 			  verschiebungsvector*=1./RealBoundaryFaces[ev];
  // 			else
  // 			  verschiebungsvector*=1./RealBoundaryFaces[ev]*sqrt(sqrt(lengthverschiebungsvector)/sqrt(lengthcell));

  // 			vertexs3d[ev]+=verschiebungsvector;
  // 			adjustadditionalvertex.insert(ev);
  // 		      }

  // 		  }
  // 	      }
  // 	  }
  // 	localadjustvertex.insert(adjustadditionalvertex.begin(),adjustadditionalvertex.end());
  // 	cout<< "iter: "<< iter<<"new adjust nodes:
  // "<<adjustadditionalvertex.size()<<endl;
  //     }
  //   cp=CellRefList.begin();

  //   for (int i=0; i<CellRefList.size(); i++)
  //     {
  // 	int hi = co2n[*cp++];
  // 	const Hex& h = hex(hi);

  // 	//if (Hexset.find(hi)==Hexset.end()) { cout<<"hi not found"<<endl;
  // continue;}
  // 	//Abfrage ob es sich um eine Randzelle handelt

  // 	// faces
  // 	for (int f=0;f<6;++f)
  // 	  {

  // 	    int fv = HexLaO.face_vertex(h,f);

  // 	    if (localadjustvertex.find(fv)!=localadjustvertex.end()) continue;

  // 	    std::array<int,4> fe;
  // 	    HexLaO.LoadEdgeVerticesOfFace(h,f,fe);
  // 	    std::array<int,4> fe_corner;
  // 	    HexLaO.globalvertices_of_face(h,fe_corner,f);

  // 	    ///////////////////////////////////////////////////////////////////////////////
  // 	    //If face is very anisotropic (||f[3]-f[1]||<<<||f[0]-f[2]|| ) the
  // algorithm will fail for slightly convex quads
  // 	    //*------------------*f[3]--------------*
  // 	    //|		       				 |
  // |
  // 	    //*f[0]--------------*midface-------f[2]*
  // 	    //|		      				 |
  // |
  // 	    //*------------------*f[1]--------------*
  // 	    // Then use for midface the midpoint between f[3] and f[1] which is
  // always in the face!

  // 	    /*fixarray<3,double> v1,v2;
  // 	      v1=vertexs3d[fe[0]]-vertexs3d[fe[2]];
  // 	      v2=vertexs3d[fe[1]]-vertexs3d[fe[3]];
  // 	      vector<double> nv(2,0);
  // 	      for(int i=0;i<4;i++)
  // 	      {
  // 	      nv[0]+=v1[i]*v1[i];
  // 	      nv[1]+=v2[i]*v2[i];
  // 	      }
  // 	      //Test if face is anisotropic
  // 	      if (nv[0]>25*nv[1])
  // 	      {
  // 	      vertexs3d[fv].equ( 0.5, vertexs3d[fe[1]], 0.5, vertexs3d[fe[3]]);
  // 	      }
  // 	      else if (nv[1]>25*nv[0])
  // 	      {
  // 	      vertexs3d[fv].equ( 0.5, vertexs3d[fe[0]],0.5, vertexs3d[fe[2]]);
  // 	      }
  // 	      /////////////////////////////////////////////////////////////////////////////
  // 	      else
  // 	      {
  // 	      vertexs3d[fv].equ(0.25, vertexs3d[fe[0]],
  // 	      0.25, vertexs3d[fe[1]],
  // 	      0.25, vertexs3d[fe[2]],
  // 	      0.25, vertexs3d[fe[3]]);
  // 	      }*/
  // 	    Vertex3d zwischen1,zwischen2;
  // 	    zwischen1.equ(0.5, vertexs3d[fe[0]],
  // 			  0.5, vertexs3d[fe[1]],
  // 			  0.5, vertexs3d[fe[2]],
  // 			  0.5, vertexs3d[fe[3]]);
  // 	    zwischen2.equ(-0.25, vertexs3d[fe_corner[0]],
  // 			  -0.25, vertexs3d[fe_corner[1]],
  // 			  -0.25, vertexs3d[fe_corner[2]],
  // 			  -0.25, vertexs3d[fe_corner[3]])	;
  // 	    vertexs3d[fv].equ(1.0,zwischen1,1.0,zwischen2);

  // 	  }

  // 	// middle
  // 	int fv = HexLaO.middle_vertex(h);
  // 	assert (localadjustvertex.find(fv)==localadjustvertex.end());
  // 	std::array<int,6> fe;
  // 	HexLaO.LoadFaceVertices(h,fe);
  // 	vector<Vertex3d> zwischen;
  // 	zwischen.resize(10);
  // 	for (int f=0;f<6;++f)
  // 	  {
  // 	    std::array<int,4> fe_corner;
  // 	    HexLaO.globalvertices_of_face(h,fe_corner,f);
  // 	    zwischen[f].equ(-0.25/12.0,vertexs3d[fe_corner[0]],
  // 			    -0.25/12.0,vertexs3d[fe_corner[1]],
  // 			    -0.25/12.0,vertexs3d[fe_corner[2]],
  // 			    -0.25/12.0,vertexs3d[fe_corner[3]]);

  // 	    int fface = HexLaO.face_vertex(h,f);

  // 	    //Falls FaceMittelpunkt projeziert -> Korrektur um diesen Vektor
  // 	    if (localadjustvertex.find(fface )!=localadjustvertex.end())
  // 	      {
  // 		std::array<int,> feedge;
  // 		HexLaO.LoadEdgeVerticesOfFace(h,f,feedge);

  // 		Vertex3d zwischen1,zwischen2;
  // 		zwischen1.equ(0.5, vertexs3d[feedge[0]],
  // 			      0.5, vertexs3d[feedge[1]],
  // 			      0.5, vertexs3d[feedge[2]],
  // 			      0.5, vertexs3d[feedge[3]]);
  // 		zwischen2.equ(-0.25, vertexs3d[fe_corner[0]],
  // 			      -0.25, vertexs3d[fe_corner[1]],
  // 			      -0.25, vertexs3d[fe_corner[2]],
  // 			      -0.25, vertexs3d[fe_corner[3]])	;
  // 		zwischen[f].equ(1.0, zwischen[f],
  // 				-1.0/4.0,zwischen1,
  // 				-1.0/4.0,zwischen2,
  // 				+1.0/4.0 ,vertexs3d[fe[f]]);
  // 	      }
  // 	  }

  // 	zwischen[6].equ(1.0,zwischen[0],
  // 			1.0,zwischen[1],
  // 			1.0,zwischen[2],
  // 			1.0,zwischen[3]);
  // 	zwischen[7].equ(1.0,zwischen[6],
  // 			1.0,zwischen[4],
  // 			1.0,zwischen[5]);

  // 	zwischen[8].equ(0.25,vertexs3d[fe[0]],
  // 			0.25,vertexs3d[fe[1]],
  // 			0.25,vertexs3d[fe[2]],
  // 			0.25,vertexs3d[fe[3]]);
  // 	zwischen[9].equ(1.0,zwischen[8],
  // 			0.25,vertexs3d[fe[4]],
  // 			0.25,vertexs3d[fe[5]]);

  // 	vertexs3d[fv].equ(1.0,zwischen[7],1.0,zwischen[9]	);

  // 	///////////////////////////////////////////////////////////////////////////////
  // 	//If face is very anisotropic the algorithm will fail for slightly
  // convex quads
  // 	/*std::array<3,double> v1,v2,v3;
  // 	  v1=vertexs3d[fe[0]]-vertexs3d[fe[5]];
  // 	  v2=vertexs3d[fe[1]]-vertexs3d[fe[3]];
  // 	  v3=vertexs3d[fe[4]]-vertexs3d[fe[2]];
  // 	  vector<double> nv(3,0);
  // 	  int sort[]={0,1,2};
  // 	  for(int i=0;i<4;i++)
  // 	  {
  // 	  nv[0]+=v1[i]*v1[i];
  // 	  nv[1]+=v2[i]*v2[i];
  // 	  nv[2]+=v3[i]*v3[i];
  // 	  }

  // 	  // Sort of the distance
  // 	  if(nv[0]>nv[1]) {double zw=nv[0]; nv[0]=nv[1]; nv[1]=zw; int
  // zwsort=sort[0]; sort[0]=sort[1];sort[1]=zwsort;} 	  if(nv[1]>nv[2])
  // {double zw=nv[1]; nv[1]=nv[2]; nv[2]=zw; int zwsort=sort[1];
  // sort[1]=sort[2];sort[2]=zwsort;} 	  if(nv[0]>nv[1]) {double zw=nv[0];
  // nv[0]=nv[1]; nv[1]=zw; int zwsort=sort[0]; sort[0]=sort[1];sort[1]=zwsort;}

  // 	  Vertex3d zw1;
  // 	  //Test if  quad is anisotropic
  // 	  vector<double> weight(3,0);
  // 	  if(nv[0]*25<nv[2] && nv[1]*25<nv[2])
  // 	  {
  // 	  for(int ii=0;ii<3;ii++) weight[ii]=0.25;
  // 	  weight[sort[2]]=0;
  // 	  }
  // 	  else if(nv[0]*25<nv[2])
  // 	  {
  // 	  for(int ii=0;ii<3;ii++) weight[ii]=0.5;
  // 	  weight[sort[1]]=0;
  // 	  weight[sort[2]]=0;
  // 	  }
  // 	  else
  // 	  {
  // 	  for(int ii=0;ii<3;ii++) weight[ii]=1./6.;
  // 	  }

  // 	  vertexs3d[fv].equ(weight[0],vertexs3d[fe[0]],
  // 	  weight[0],vertexs3d[fe[5]],
  // 	  weight[1],vertexs3d[fe[1]],
  // 	  weight[1],vertexs3d[fe[3]],
  // 	  weight[2],vertexs3d[fe[4]],
  // 	  weight[2],vertexs3d[fe[2]]);
  // 	*/

  // 	//////////////////////////////////////////////////////////////////
  // 	/*else
  // 	  {
  // 	  for (int i=0;i<6;++i) vertexs3d[fv]+=vertexs3d[fe[i]];
  // 	  vertexs3d[fv]*=1./6.;
  // 	  }*/

  // 	// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEU

  // 	// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ALT
  // 	//       for (int face=0; face<6; face++)
  // 	// 	{
  // 	// 	  HexLaO.LoadEdgeVerticesOfFace(h,face,v);
  // 	// 	  int mv = HexLaO.face_vertex(h,face);
  // 	// 	  if (adjustvertex.find(mv)!=adjustvertex.end()) continue;
  // 	// 	  new_face_vertex3d(mv,v);
  // 	// 	}
  // 	// //       std::array<6,int> w;
  // 	// //       int mv = HexLaO.middle_vertex(h);
  // 	// //       HexLaO.LoadFaceVertices(h,w);
  // 	// //       new_vertex3d(mv,w);
  // 	// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ALT
  //     }
  // }

  /*---------------------------------------------------*/
};

#endif
