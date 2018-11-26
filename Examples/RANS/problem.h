#ifndef PROBLEM_H
#define PROBLEM_H

#include <boundaryfunctional.h>
#include <boundaryfunction.h>
#include <meshagent.h>
#include <navierstokeslps2d.h>
#include <problemdescriptorbase.h>
#include <usefullfunctionsbd.h>
#include "rans_lsq.h"
#include "rans_twomodels.h"
#include "rans_lin.h"
#include "rans_const.h"

using namespace Gascoigne;

/*----------------------------------------------------------------------------*/

class MyDirichletData : public DirichletData
{
public:
    MyDirichletData()
    {
    }

    std::string GetName() const
    {
        return "Wall Mount Problem";
    }

    void operator()(DoubleVector& b, const Vertex2d& v, int color) const
    {
        double y = v.y();

        b.zero();
        if (color == 8)
        {
            double vmax = 3.0;

            double sc = 1.0;
            if (GetTime() < 0.5)
            {
                sc = 0.5 - 0.5 * cos(2 * M_PI * GetTime());
            }

            b[1] = sc * vmax * ParabelFunction(y, 0., 1);
        }
    }

protected:
    double vmax;
};

/*----------------------------------------------------------------------------*/

class Drag : public virtual BoundaryFunctional
{
    double _visc;

public:
    Drag(const ParamFile* pf) : BoundaryFunctional()
    {
        DataFormatHandler DFH;
        DFH.insert("visc", &_visc, 0.1);

        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(pf, "Equation");
    }

    ~Drag()
    {
    }

    std::string GetName() const
    {
        return "Drag";
    }

    double J(const FemFunction& U, const Vertex2d& v, const Vertex2d& n, int color) const
    {
        // Implement the drag functional here
        if (color == 6)
        {
            return (U[1].x() - U[0].m()) * n.x() + (U[1].y() + U[2].x()) * n.y();
        }
        else
            return 0.;
    }
};

// force = (-p*I + 2*nu*D)*n
// D=force[0]*ds(1)
// L=force[1]*ds(1)
//  plot 'functional.dat' using 1:2 title 'lift', 'functional.dat' using 1:3 title 'drag'

/*----------------------------------------------------------------------------*/

class Lift : public virtual BoundaryFunctional
{
    double _visc;

public:
    Lift(const ParamFile* pf) : BoundaryFunctional()
    {
        DataFormatHandler DFH;
        DFH.insert("visc", &_visc, 0.1);

        FileScanner FS(DFH);
        FS.NoComplain();
        FS.readfile(pf, "Equation");
    }

    ~Lift()
    {
    }

    std::string GetName() const
    {
        return "Lift";
    }

    double J(const FemFunction& U, const Vertex2d& v, const Vertex2d& n, int color) const
    {
        // Implement the lift functional here
        if (color == 6)
        {
            return (U[2].y() - U[0].m()) * n.y() + (U[1].y() + U[2].x()) * n.x();
        }
        else
            return 0.;
    }
};

/*----------------------------------------------------------------------------*/

class MyProblemDescriptor : public ProblemDescriptorBase
{
public:
    std::string GetName() const
    {
        return "Benchmark Problem";
    }

    void BasicInit(const ParamFile* pf)
    {
        GetParamFilePointer() = pf;
        // GetEquationPointer()  = new rans_base<2>(GetParamFile());
        GetEquationPointer() = new rans_lsq<2, 10, 3>(GetParamFile());
        // GetBoundaryEquationPointer() = new rans_lsq<2, 10, 3>(GetParamFile());
        GetDirichletDataPointer() = new MyDirichletData;
        ProblemDescriptorBase::BasicInit(pf);
    }
};

class RTMProblemDescriptor : public ProblemDescriptorBase
{
public:
    std::string GetName() const
    {
        return "Benchmark Problem";
    }

    void BasicInit(const ParamFile* pf)
    {
        GetParamFilePointer() = pf;
        // GetEquationPointer()  = new rans_base<2>(GetParamFile());
        GetEquationPointer()      = new rans_twomodels<2, 7, 2>(GetParamFile());
        GetDirichletDataPointer() = new MyDirichletData;
        ProblemDescriptorBase::BasicInit(pf);
    }
};

class NSProblemDescriptor : public ProblemDescriptorBase
{
public:
    std::string GetName() const
    {
        return "Benchmark Problem";
    }

    void BasicInit(const ParamFile* pf)
    {
        GetParamFilePointer()     = pf;
        GetEquationPointer()      = new rans_base<2>(GetParamFile());
        GetDirichletDataPointer() = new MyDirichletData;
        ProblemDescriptorBase::BasicInit(pf);
    }
};

class LinProblemDescriptor : public ProblemDescriptorBase
{
public:
    std::string GetName() const
    {
        return "Benchmark Problem";
    }

    void BasicInit(const ParamFile* pf)
    {
        GetParamFilePointer() = pf;
        // GetEquationPointer()  = new rans_base<2>(GetParamFile());
        GetEquationPointer()      = new rans_twomodels<2, 7, 1>(GetParamFile());
        GetDirichletDataPointer() = new MyDirichletData;
        ProblemDescriptorBase::BasicInit(pf);
    }
};

class ConstProblemDescriptor : public ProblemDescriptorBase
{
public:
    std::string GetName() const
    {
        return "Benchmark Problem";
    }

    void BasicInit(const ParamFile* pf)
    {
        GetParamFilePointer()     = pf;
        GetEquationPointer()      = new rans_const<2>(GetParamFile());
        GetDirichletDataPointer() = new MyDirichletData;
        ProblemDescriptorBase::BasicInit(pf);
    }
};

/*----------------------------------------------------------------------------*/

// class MyDirichletData : public DirichletData
// {
// public:
//   MyDirichletData()
//   {
//   }
//
//   std::string GetName() const
//   {
//     return "Benchmark Problem";
//   }
//
//   void operator()(DoubleVector &b, const Vertex2d &v, int color) const
//   {
//     double y = v.y();
//
//     b.zero();
//     if (color != 80) {
//       double high = 4.1;
//
//       double vmax;
//       if (GetTime()>0.1) vmax = 1.5;
//       else vmax = 15.*GetTime();
//
//       b[1] = vmax * ParabelFunction(y, 0., high);
//     }
//   }
//
// protected:
//   double vmax;
// };
//
//
// /*----------------------------------------------------------------------------*/
//
//
// class Drag : public virtual BoundaryFunctional{
//
//   double _visc;
//
// public:
//
//   Drag(const ParamFile* pf): BoundaryFunctional(){
//
//     DataFormatHandler DFH;
//     DFH.insert("visc" ,&_visc , 0.1);
//
//     FileScanner FS(DFH);
//     FS.NoComplain();
//     FS.readfile(pf,"Equation");
// 	}
//
//   ~Drag(){}
//
//   std::string GetName() const {
//      return "Drag";
//   }
//
//   double J(const FemFunction& U, const Vertex2d& v, const Vertex2d& n, int color) const {
//     //Implement the drag functional here
//     if (color==80) {
//      return  (U[1].x()  - U[0].m())*n.x() + (U[1].y() + U[2].x())*n.y()  ;
//     } else return 0.;
//    	//return 0.;
//   }
// };
//
// //force = (-p*I + 2*nu*D)*n
// //D=force[0]*ds(1)
// //L=force[1]*ds(1)
// //  plot 'functional.dat' using 1:2 title 'lift', 'functional.dat' using 1:3 title 'drag'
//
//
// /*----------------------------------------------------------------------------*/
//
//
// class Lift : public virtual BoundaryFunctional{
//
//   double _visc;
//
//  public:
//
//   Lift(const ParamFile* pf): BoundaryFunctional(){
//
//     DataFormatHandler DFH;
//     DFH.insert("visc", &_visc, 0.1);
//
//     FileScanner FS(DFH);
//     FS.NoComplain();
//     FS.readfile(pf,"Equation");
//
//   }
//
//
//   ~Lift(){}
//
//   std::string GetName() const {
//      return "Lift";
//   }
//
//   double J(const FemFunction& U, const Vertex2d& v, const Vertex2d& n, int color) const {
//    //Implement the lift functional here
//        if (color==80) {
//      return  (U[2].y() - U[0].m())*n.y() + (U[1].y() + U[2].x())*n.x() ;
//     } else return 0.;
//   }
// };
//
//
//
// /*----------------------------------------------------------------------------*/
//
// class MyProblemDescriptor : public ProblemDescriptorBase
// {
// public:
//   std::string GetName() const
//   {
//     return "Benchmark Problem";
//   }
//
//   void BasicInit(const ParamFile *pf)
//   {
//     GetParamFilePointer() = pf;
//     GetEquationPointer() = new rans_lsq(GetParamFile());//7.1
//     //GetEquationPointer() = new MyLpsEquation(GetParamFile());//7.2
//     GetDirichletDataPointer() = new MyDirichletData;
//
//     ProblemDescriptorBase::BasicInit(pf);
//   }
// };
//
// /*----------------------------------------------------------------------------*/
//
// class Disc : public BoundaryFunction<2>
// {
// public:
//   std::string GetName() const
//   {
//     return "Disc";
//   }
//
//   void BasicInit(Vertex2d c, double r)
//   {
//     _c = c;
//     _r = r;
//   }
//
//   double operator()(const Vertex2d &c) const
//   {
//     double r = -_r;
//     for (int i = 0; i < 2; i++) {
//       double dx = c[i] - _c[i];
//       r += dx * dx;
//     }
//     return r;
//   }
//
// private:
//   double _r;
//   Vertex2d _c;
// };
//
// /*----------------------------------------------------------------------------*/
//
// class MyMeshAgent : public MeshAgent
// {
// public:
//   MyMeshAgent() : MeshAgent()
//   {
//     double r = 0.25;
//     Vertex2d v(2., 2.);
//     D.BasicInit(v, r);
//
//     AddShape(80, &D);
//   }
//
// private:
//   Disc D;
// };

#endif /* PROBLEM_H */
