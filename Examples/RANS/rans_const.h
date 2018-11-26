#ifndef RANS_CONST_H
#define RANS_CONST_H

#include <cmath>
#include <limits>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "rans_base.h"
#include "paramfile.h"
#include "polynomial_features.h"

/*-----------------------------------------*/

namespace Gascoigne
{
//! \brief RANS implementation with LPS stabilization.
//!
//! RST is calculated in two parts: linear and nonlinear part. Linear part is computed with linear
//! least squares solution and the nonlinear part as least squares solution with polynomial degree
//! degree and n_features features.
template <int DIM>
class rans_const : public rans_base<DIM>
{
public:
    rans_const();
    rans_const(const ParamFile* pf);

    std::string GetName() const
    {
        return "Rans ML RST with Elastic Net";
    }

    void point(double h, const FemFunction& U, const Vertex2d& v) const;

    ~rans_const();

protected:
    using matrix_rans = typename rans_base<DIM>::matrix_rans;
    using vector_rans = typename rans_base<DIM>::vector_rans;

    void get_rst() const;

private:
    mutable double _turb_visc;
    using rans_base<DIM>::rst_size;
    using rans_base<DIM>::_rst;
    using rans_base<DIM>::_h;
    using rans_base<DIM>::_old;

    mutable matrix_rans _nabla_velocity;
    mutable matrix_rans _strain_rate;
    mutable matrix_rans _rotation_rate;
};

/*-----------------------------------------*/

template <int DIM>
rans_const<DIM>::~rans_const()
{
}

template <int DIM>
rans_const<DIM>::rans_const()
    : rans_base<DIM>()
    , _nabla_velocity(matrix_rans::Zero())
    , _strain_rate(matrix_rans::Zero())
    , _rotation_rate(matrix_rans::Zero())
{
}

template <int DIM>
rans_const<DIM>::rans_const(const ParamFile* pf)
    : rans_base<DIM>(pf)
    , _nabla_velocity(matrix_rans::Zero())
    , _strain_rate(matrix_rans::Zero())
    , _rotation_rate(matrix_rans::Zero())
{
    DataFormatHandler DFH;
    DFH.insert("turb_visc", &_turb_visc, 0.0);
    FileScanner FS(DFH, pf, "Equation");
}

template <int DIM>
void rans_const<DIM>::point(double h, const FemFunction& U, const Vertex2d& v) const
{
    _h = h;
    Multiplex::init_NV<DIM>(_nabla_velocity, *_old);
    Multiplex::init_strain_rot<DIM>(_strain_rate, _rotation_rate, *_old);
    get_rst();
}

template <int DIM>
void rans_const<DIM>::get_rst() const
{
    auto sc = sqrt(_nabla_velocity.cwiseProduct(_nabla_velocity).trace());

    if (abs(sc) > std::numeric_limits<double>::epsilon())
    {
        _rst = _turb_visc * _strain_rate / sc;
    }
    else
    {
        _rst = _turb_visc * _strain_rate;
    }
}
}  // namespace Gascoigne

#endif
