#ifndef RANS_LIN_H
#define RANS_LIN_H

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
template <int DIM, int n_features, int degree>
class rans_lin : public rans_base<DIM>
{
public:
    rans_lin();
    rans_lin(const ParamFile* pf);

    std::string GetName() const
    {
        return "Rans ML RST with Elastic Net";
    }

    void point(double h, const FemFunction& U, const Vertex2d& v) const;
    void lpspoint(double h, const FemFunction& U, const Vertex2d& v) const;

    ~rans_lin();

protected:
    using matrix_rans = typename rans_base<DIM>::matrix_rans;
    using vector_rans = typename rans_base<DIM>::vector_rans;
    using vector_feat = typename Eigen::Matrix<double, n_features, 1>;

    void get_rst() const;

private:
    using poly_feat_coll = poly_feat_collect<n_features, degree, double>;
    const poly_feat_coll _poly_feat{};
    static constexpr auto feature_vec_length = poly_feat_coll::feature_vec_length;
    using rans_base<DIM>::rst_size;
    using rans_base<DIM>::_rst;
    using rans_base<DIM>::_h;
    using rans_base<DIM>::_old;
    using rans_base<DIM>::_visc;
    using rans_base<DIM>::_visc0;
    using rans_base<DIM>::_alpha0;
    using rans_base<DIM>::_alpha;
    using rans_base<DIM>::_delta0;
    using rans_base<DIM>::_delta;

    using vector_poly_feat = Eigen::Matrix<double, feature_vec_length, 1>;
    mutable matrix_rans _nabla_velocity;
    mutable matrix_rans _strain_rate;
    mutable matrix_rans _rotation_rate;
    mutable vector_rans _nabla_pressure;
    mutable vector_feat _feature_vector;
    mutable vector_poly_feat _feature_poly_vector;
    mutable double _turb_visc;

    // Maybe use sparse matrices here, because of L1 regularization we usually get sparse models
    using matrix_linear_model            = Eigen::Matrix<double, 1, n_features>;
    const matrix_linear_model _lin_model = load_matrix<double, 1, n_features>("coefs_lin.csv");

    static inline double v_f1(const matrix_rans& strain_rate, const matrix_rans& rotation_rate,
                              const vector_rans& nabla_pressure, const matrix_rans& nabla_velocity);

    static inline double v_f2(const matrix_rans& strain_rate, const matrix_rans& rotation_rate,
                              const vector_rans& nabla_pressure, const matrix_rans& nabla_velocity);

    static inline double v_f3(const matrix_rans& strain_rate, const matrix_rans& rotation_rate,
                              const vector_rans& nabla_pressure, const matrix_rans& nabla_velocity);

    static inline double v_f4(const matrix_rans& strain_rate, const matrix_rans& rotation_rate,
                              const vector_rans& nabla_pressure, const matrix_rans& nabla_velocity);

    static inline double v_f5(const matrix_rans& strain_rate, const matrix_rans& rotation_rate,
                              const vector_rans& nabla_pressure, const matrix_rans& nabla_velocity);

    static inline double v_q_criterion(const matrix_rans& strain_rate,
                                       const matrix_rans& rotation_rate,
                                       const vector_rans& nabla_pressure,
                                       const matrix_rans& nabla_velocity);

    //! \brief Computing the features, implement feature computations as static member functions
    //! with arguments: strain_rate, rotation_rate, nabla_pressure, nabla_velocity
    template <typename... f>
    void feature_computation(const f&... features) const;
};

/*-----------------------------------------*/

template <int DIM, int n_features, int degree>
rans_lin<DIM, n_features, degree>::~rans_lin()
{
}

template <int DIM, int n_features, int degree>
rans_lin<DIM, n_features, degree>::rans_lin()
    : rans_base<DIM>()
    , _nabla_velocity(matrix_rans::Zero())
    , _strain_rate(matrix_rans::Zero())
    , _rotation_rate(matrix_rans::Zero())
    , _nabla_pressure(vector_rans::Zero())
    , _feature_vector(vector_feat::Zero())
{
}

template <int DIM, int n_features, int degree>
rans_lin<DIM, n_features, degree>::rans_lin(const ParamFile* pf)
    : rans_base<DIM>(pf)
    , _nabla_velocity(matrix_rans::Zero())
    , _strain_rate(matrix_rans::Zero())
    , _rotation_rate(matrix_rans::Zero())
    , _nabla_pressure(vector_rans::Zero())
    , _feature_vector(vector_feat::Zero())
{
}

template <int DIM, int n_features, int degree>
void rans_lin<DIM, n_features, degree>::point(double h, const FemFunction& U,
                                              const Vertex2d& v) const
{
    _h = h;
    Multiplex::init_NV<DIM>(_nabla_velocity, *_old);
    Multiplex::init_strain_rot<DIM>(_strain_rate, _rotation_rate, *_old);
    Multiplex::init_nP<DIM>(_nabla_pressure, *_old);
    get_rst();
    _visc = _visc0 + _turb_visc;
}

template <int DIM, int n_features, int degree>
void rans_lin<DIM, n_features, degree>::lpspoint(double h, const FemFunction& U,
                                                 const Vertex2d& v) const
{
    _visc        = _visc0 + _turb_visc;
    double _norm = std::sqrt(U[1][0] * U[1][0] + U[2][0] * U[2][0]);
    double val   = 6. * _visc / (h * h) + _norm / h;
    _delta       = _delta0 / val;
    _alpha       = _alpha0 / val;
}

template <int DIM, int n_features, int degree>
template <typename... f>
void rans_lin<DIM, n_features, degree>::feature_computation(const f&... features) const
{
    auto feature_ptr = _feature_vector.data();
    for (auto feature : {features...})
    {
        *feature_ptr = feature(_strain_rate, _rotation_rate, _nabla_pressure, _nabla_velocity);
        feature_ptr++;
    }
    _feature_poly_vector = _poly_feat(_feature_vector);
    auto iterator        = _feature_poly_vector.data();
    for (size_t i = 0; i < feature_vec_length; i++)
    {
        if (isnan(*iterator))
        {
            *iterator = 0;
        }
        iterator++;
    }
}

template <int DIM, int n_features, int degree>
void rans_lin<DIM, n_features, degree>::get_rst() const
{
    feature_computation(v_f1, v_f2, v_f3, v_f4, v_f5, v_q_criterion);
    auto sc    = sqrt(_nabla_velocity.cwiseProduct(_nabla_velocity).trace());
    _turb_visc = _lin_model.dot(_feature_poly_vector);

    if (abs(sc) > std::numeric_limits<double>::epsilon())
    {
        _rst = _turb_visc * _strain_rate / sc;
    }
    else
    {
        _rst = _turb_visc * _strain_rate;
        //_rst = matrix_rans::Zero();
    }
}

template <int DIM, int n_features, int degree>
double rans_lin<DIM, n_features, degree>::v_f1(const matrix_rans& strain_rate,
                                               const matrix_rans& rotation_rate,
                                               const vector_rans& nabla_pressure,
                                               const matrix_rans& nabla_velocity)
{
    return (strain_rate * strain_rate).trace() / nabla_velocity.norm();
}

template <int DIM, int n_features, int degree>
double rans_lin<DIM, n_features, degree>::v_f2(const matrix_rans& strain_rate,
                                               const matrix_rans& rotation_rate,
                                               const vector_rans& nabla_pressure,
                                               const matrix_rans& nabla_velocity)
{
    auto sc = nabla_velocity.norm();
    if (abs(sc) < std::numeric_limits<double>::epsilon())
    {
        return 0;
    }
    return (rotation_rate * rotation_rate).trace() / sc;
}

template <int DIM, int n_features, int degree>
double rans_lin<DIM, n_features, degree>::v_f3(const matrix_rans& strain_rate,
                                               const matrix_rans& rotation_rate,
                                               const vector_rans& nabla_pressure,
                                               const matrix_rans& nabla_velocity)
{
    auto sc = nabla_velocity.norm();
    if (abs(sc) < std::numeric_limits<double>::epsilon())
    {
        return 0;
    }
    return (strain_rate * strain_rate * strain_rate).trace() / sc;
}

template <int DIM, int n_features, int degree>
double rans_lin<DIM, n_features, degree>::v_f4(const matrix_rans& strain_rate,
                                               const matrix_rans& rotation_rate,
                                               const vector_rans& nabla_pressure,
                                               const matrix_rans& nabla_velocity)
{
    auto sc = nabla_velocity.norm();
    if (abs(sc) < std::numeric_limits<double>::epsilon())
    {
        return 0;
    }
    return (rotation_rate * rotation_rate * strain_rate).trace() / sc;
}

template <int DIM, int n_features, int degree>
double rans_lin<DIM, n_features, degree>::v_f5(const matrix_rans& strain_rate,
                                               const matrix_rans& rotation_rate,
                                               const vector_rans& nabla_pressure,
                                               const matrix_rans& nabla_velocity)
{
    auto sc = nabla_velocity.norm();
    if (abs(sc) < std::numeric_limits<double>::epsilon())
    {
        return 0;
    }
    return (strain_rate * strain_rate * rotation_rate * rotation_rate).trace() / sc;
}

template <int DIM, int n_features, int degree>
double rans_lin<DIM, n_features, degree>::v_q_criterion(const matrix_rans& strain_rate,
                                                        const matrix_rans& rotation_rate,
                                                        const vector_rans& nabla_pressure,
                                                        const matrix_rans& nabla_velocity)
{
    auto rr = abs((rotation_rate * rotation_rate.transpose()).trace());
    auto ss = abs((strain_rate * strain_rate.transpose()).trace());
    return (rr - ss) / (rr + ss);
}

}  // namespace Gascoigne

#endif
