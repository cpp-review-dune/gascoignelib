
#ifndef __local_h
#define __local_h

#include "fsi.h"
#include "div.h"
#include "dirichletdata.h"
#include "problemdescriptorbase.h"
#include "domainrighthandside.h"
#include "componentinformationbase.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace Gascoigne;

/* ----------------------------------------- */

extern double __TIME;

template <int DIM>
class FSI_CI : public ComponentInformationBase
{
public:
    void BasicInit(const ParamFile* pf) override
    {
    }

    std::string GetName() const override final
    {
        return "FSI CI";
    }

    const int GetNScalars() const override final
    {
        return 2 * DIM + 1;
    }

    void GetScalarName(int i, std::string& s_name) const override final
    {
        if (i == 0)
            s_name = "p";
        if (DIM == 2)
        {
            if (i == 1)
                s_name = "vx";
            if (i == 2)
                s_name = "vy";
            if (i == 3)
                s_name = "ux";
            if (i == 4)
                s_name = "uy";
        }
        if (DIM == 3)
        {
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

    const int GetNVectors() const override final
    {
        return 2;
    }
    void GetVectorName(int i, std::string& s_name) const override final
    {
        if (i == 0)
            s_name = "V";
        else if (i == 1)
            s_name = "U";
        else
            abort();
    }
    void GetVectorIndices(int i, std::array<int, 3>& fa_vectorindices) const override final
    {
        if (i == 0)
        {
            fa_vectorindices[0] = 1;
            fa_vectorindices[1] = 2;
            fa_vectorindices[2] = -1;
            if (DIM == 3)
                fa_vectorindices[3] = 3;
        }
        else
        {
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

    std::string GetName() const
    {
        return "MyDD";
    }

    void operator()(DoubleVector& b, const Vertex2d& v, int color) const
    {
        b.zero();

        double t = GetTime();

        double sc = 1.0;
        if (t < 1.0)
            sc = 0.5 - 0.5 * cos(M_PI * t);
        double veff = vmean * sc;

        // veff *= veff;
        // veff *= 0.5;
        // for wall_mount color 8, for benchmark color 0

        // if (color == 0)
        //     b[1] += v.y() * (0.41 - v.y()) / 0.205 / 0.205 * veff * 1.5;

        if (color == 8)
            b[1] += v.y() * (1 - v.y()) / 0.5 / 0.5 * veff * 1.5;
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

    std::string GetName() const
    {
        return "MyDD";
    }

    void operator()(DoubleVector& b, const Vertex3d& v, int color) const
    {
        b.zero();
        double sc = 1.0;
        double t  = GetTime();
        if (t < 2.0)
            sc = 0.5 * (1.0 - cos(M_PI * t / 2.0));

        if (color == 0)
            b[1] += v.y() * (0.4 - v.y()) / 0.2 / 0.2 * (0.4 - v.z()) * (0.4 + v.z()) / 0.4 / 0.4
                    * vmean * sc * 9.0 / 8.0;
    }
};

// -----------------------------------------

class ProblemDescriptor2d : public ProblemDescriptorBase
{
public:
    std::string GetName() const
    {
        return "fsi";
    }
    void BasicInit(const ParamFile* pf)
    {
        GetParamFilePointer() = pf;
        GetEquationPointer()  = new FSI<2>(GetParamFile());
        //    GetBoundaryEquationPointer() = new FSI<2>(GetParamFile());
        GetDirichletDataPointer() = new MyDD(GetParamFile());

        ProblemDescriptorBase::BasicInit(pf);

        GetComponentInformationPointer() = new FSI_CI<2>;
    }
};

class DivergenceDescriptor2d : public ProblemDescriptorBase
{
public:
    std::string GetName() const
    {
        return "div";
    }
    void BasicInit(const ParamFile* pf)
    {
        GetParamFilePointer() = pf;
        GetEquationPointer()  = new DIV_proj<2>(GetParamFile());
        //    GetBoundaryEquationPointer() = new FSI<2>(GetParamFile());
        // ?
        GetDirichletDataPointer() = new MyDD(GetParamFile());

        ProblemDescriptorBase::BasicInit(pf);

        GetComponentInformationPointer() = new FSI_CI<2>;
    }
};

class ProblemDescriptor3d : public ProblemDescriptorBase
{
public:
    std::string GetName() const
    {
        return "fsi";
    }
    void BasicInit(const ParamFile* pf)
    {
        GetParamFilePointer() = pf;
        GetEquationPointer()  = new FSI<3>(GetParamFile());
        //    GetBoundaryEquationPointer() = new FSI<3>(GetParamFile());
        GetDirichletDataPointer() = new MyDD3d(GetParamFile());

        ProblemDescriptorBase::BasicInit(pf);

        GetComponentInformationPointer() = new FSI_CI<3>;
    }
};

#endif
