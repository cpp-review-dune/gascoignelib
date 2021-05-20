
#ifndef __FSI_main_dirichlet_h
#define __FSI_main_dirichlet_h

#include "dirichletdata.h"
#include "filescanner.h"
#include "paramfile.h"

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

extern double __TIME;

class FSI_main_MyDD : public DirichletData
{
protected:
  double vmean;

public:
  FSI_main_MyDD(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("vmean", &vmean, 0.0);
    FileScanner FS(DFH, pf, "Equation");
  }

  std::string GetName() const { return "MyDD"; }

  void operator()(DoubleVector& b, const Vertex2d& v, int color) const
  {
    b.zero();

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

    double sc = 1.0;
    // // transition from zero to vmean at time t=0 (2 sec)
    if (__TIME < 2.0) {
      sc = 0.5 * (1.0 - cos(M_PI * __TIME / 2.0));
      if (color == 0)
        b[1] += v.y() * (0.41 - v.y()) / 0.205 / 0.205 * vmean * 1.5 * sc;
    }
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
    else {
      sc = 1.0;
      if (color == 0)
        b[1] += v.y() * (0.41 - v.y()) / 0.205 / 0.205 * vmean * 1.5 * sc;
    }

    // // transition from 2.0 to vmean at time t=10 (1 sec)
    // if (__TIME>10.0)
    //   {
    // 	sc = 0.5*(1.0-cos(M_PI*__TIME/2.0));

    //   }
  }
};

class FSI_main_MyDD3d : public DirichletData
{
protected:
  double vmean;

public:
  FSI_main_MyDD3d(const ParamFile* pf)
  {
    DataFormatHandler DFH;
    DFH.insert("vmean", &vmean, 0.0);
    FileScanner FS(DFH, pf, "Equation");
  }

  std::string GetName() const { return "MyDD"; }

  void operator()(DoubleVector& b, const Vertex3d& v, int color) const
  {
    b.zero();
    // inflow radius
    double radius = 0.3;
    // inflow in ml/s
    double flow = 100. / 60.;
    double fact_time = 2 * flow / (radius * radius * M_PI);
    // cout<< "fact_time: "<<fact_time<<endl;
    // if (__TIME<1.0) fact_time *= 0.5*(1.0-cos(M_PI*__TIME));

    double x = v.x();
    double y = v.y();
    double z = v.z();
    /*  double fact_time=0;
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
                  fact_time=fact_time*1.0e-1;
  */

    if (color == 7) {
      Vertex3d center;
      center[0] = 3.78223304703;
      center[1] = -3.882233047;
      center[2] = 1.0928932;
      Vertex3d normal;
      normal[0] = 0.5000;
      normal[1] = -0.4998;
      normal[2] = 0.7072;

      b[1] +=
        -normal[0] * fact_time *
        ((x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) +
         (z - center[2]) * (z - center[2]) - radius * radius) /
        (-radius * radius);
      b[2] +=
        -normal[1] * fact_time *
        ((x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) +
         (z - center[2]) * (z - center[2]) - radius * radius) /
        (-radius * radius);
      b[3] +=
        -normal[2] * fact_time *
        ((x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) +
         (z - center[2]) * (z - center[2]) - radius * radius) /
        (-radius * radius);
    }
  }
};

#endif
