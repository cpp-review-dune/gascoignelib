#ifndef RESFUNCTIONAL_H
#define RESFUNCTIONAL_H

#include "local.h"
#include "residualfunctional.h"
#include "dirichletdatabycolor.h"

class Drag : public virtual ResidualFunctional
{
    std::string GetName() const
    {
        return "drag";
    }

public:
    Drag()
    {
        __comps.push_back(1);
        __scales.push_back(1.0);
        __cols.insert(80);
        //__cols.insert(81);
        __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
    }
};

class Lift : public virtual ResidualFunctional
{
    std::string GetName() const
    {
        return "lift";
    }

public:
    Lift()
    {
        __comps.push_back(2);
        __scales.push_back(1.0);
        __cols.insert(80);
        //__cols.insert(81);
        __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
    }
};

class Drag3d : public virtual ResidualFunctional
{
    std::string GetName() const
    {
        return "drag";
    }

public:
    Drag3d()
    {
        __comps.push_back(1);
        __scales.push_back(1.0);
        __cols.insert(84);
        __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
    }
};
class Lift3d : public virtual ResidualFunctional
{
    std::string GetName() const
    {
        return "drag";
    }

public:
    Lift3d()
    {
        __comps.push_back(2);
        __scales.push_back(1.0);
        __cols.insert(84);
        __DD = new DirichletDataByColor(GetComps(), GetColors(), GetScales());
    }
};

#endif
