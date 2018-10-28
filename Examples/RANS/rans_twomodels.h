#ifndef RANS_TWOMODELS_H
#define RANS_TWOMODELS_H

#include <cmath>
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
class rans_twomodels : public rans_base<DIM>
{
public:
    rans_twomodels();
    rans_twomodels(const ParamFile* pf);

    std::string GetName() const
    {
        return "Rans ML RST with Elastic Net";
    }

    void point(double h, const FemFunction& U, const Vertex2d& v) const;

    ~rans_twomodels();

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

    using vector_poly_feat = Eigen::Matrix<double, feature_vec_length, 1>;
    mutable matrix_rans _nabla_velocity;
    mutable matrix_rans _strain_rate;
    mutable matrix_rans _rotation_rate;
    mutable vector_rans _nabla_pressure;
    mutable vector_feat _feature_vector;
    mutable vector_poly_feat _feature_poly_vector;

    // Maybe use sparse matrices here, because of L1 regularization we usually get sparse models
    using matrix_linear_model =
      Eigen::Matrix<double, matrix_rans::RowsAtCompileTime, matrix_rans::ColsAtCompileTime>;
    const matrix_linear_model _lin_model =
      load_matrix<double, matrix_rans::RowsAtCompileTime, matrix_rans::ColsAtCompileTime>(
        "coefs.csv");

    using matrix_model        = Eigen::Matrix<double, rst_size, feature_vec_length>;
    const matrix_model _model = load_matrix<double, rst_size, feature_vec_length>("coefs.csv");
    const Eigen::Matrix<double, n_features, 1> mean_var_normalizations =
      load_matrix<double, 1, n_features>("means.csv").transpose();

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

    static inline double p_f8(const matrix_rans& strain_rate, const matrix_rans& rotation_rate,
                              const vector_rans& nabla_pressure, const matrix_rans& nabla_velocity);

    static inline double v_q_criterion(const matrix_rans& strain_rate,
                                       const matrix_rans& rotation_rate,
                                       const vector_rans& nabla_pressure,
                                       const matrix_rans& nabla_velocity);

    static inline double v_rot_rate_mag(const matrix_rans& strain_rate,
                                        const matrix_rans& rotation_rate,
                                        const vector_rans& nabla_pressure,
                                        const matrix_rans& nabla_velocity);

    static inline double v_strain_rate_mag(const matrix_rans& strain_rate,
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
rans_twomodels<DIM, n_features, degree>::~rans_twomodels()
{
}

template <int DIM, int n_features, int degree>
rans_twomodels<DIM, n_features, degree>::rans_twomodels()
    : rans_base<DIM>()
    , _nabla_velocity(matrix_rans::Zero())
    , _strain_rate(matrix_rans::Zero())
    , _rotation_rate(matrix_rans::Zero())
    , _nabla_pressure(vector_rans::Zero())
    , _feature_vector(vector_feat::Zero())
{
}

template <int DIM, int n_features, int degree>
rans_twomodels<DIM, n_features, degree>::rans_twomodels(const ParamFile* pf)
    : rans_base<DIM>(pf)
    , _nabla_velocity(matrix_rans::Zero())
    , _strain_rate(matrix_rans::Zero())
    , _rotation_rate(matrix_rans::Zero())
    , _nabla_pressure(vector_rans::Zero())
    , _feature_vector(vector_feat::Zero())
{
}

template <int DIM, int n_features, int degree>
void rans_twomodels<DIM, n_features, degree>::point(double h, const FemFunction& U,
                                                    const Vertex2d& v) const
{
    _h = h;
    Multiplex::init_NV<DIM>(_nabla_velocity, *_old);
    Multiplex::init_strain_rot<DIM>(_strain_rate, _rotation_rate, *_old);
    Multiplex::init_nP<DIM>(_nabla_pressure, *_old);

    feature_computation(v_f1, v_f2, v_f3, v_f4, v_f5, p_f8, v_q_criterion, v_rot_rate_mag,
                        v_strain_rate_mag);
    // scale_features();
    get_rst();
}

template <int DIM, int n_features, int degree>
template <typename... f>
void rans_twomodels<DIM, n_features, degree>::feature_computation(const f&... features) const
{
    auto feature_ptr = _feature_vector.data();
    for (auto feature : {features...})
    {
        *feature_ptr = feature(_strain_rate, _rotation_rate, _nabla_pressure, _nabla_velocity);
        feature_ptr++;
    }
    _feature_vector -= mean_var_normalizations.col(0);
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
void rans_twomodels<DIM, n_features, degree>::get_rst() const
{
    std::cout << _feature_vector << '\n';
    auto lin_rst = (_lin_model.cwiseProduct(_strain_rate)).sum();
    auto tmp_rst = _model * _feature_poly_vector;
    for (size_t j = 0; j < DIM; j++)
    {
        for (size_t i = 0; i < j + 1; i++)
        {
            _rst(i, j) = tmp_rst[j + (2 * DIM - 1 - i) * i / 2];
            _rst(j, i) = tmp_rst[j + (2 * DIM - 1 - i) * i / 2];
        }
    }
    _rst = lin_rst * _strain_rate;
    // _rst =
    //   4e-4 * _strain_rate - (1 / 3) * _nabla_velocity.cwiseProduct(matrix_rans::Identity());
    // std::cout << _rst;
}

template <int DIM, int n_features, int degree>
double rans_twomodels<DIM, n_features, degree>::v_f1(const matrix_rans& strain_rate,
                                                     const matrix_rans& rotation_rate,
                                                     const vector_rans& nabla_pressure,
                                                     const matrix_rans& nabla_velocity)
{
    return (strain_rate * strain_rate).trace();
}

template <int DIM, int n_features, int degree>
double rans_twomodels<DIM, n_features, degree>::v_f2(const matrix_rans& strain_rate,
                                                     const matrix_rans& rotation_rate,
                                                     const vector_rans& nabla_pressure,
                                                     const matrix_rans& nabla_velocity)
{
    return (rotation_rate * rotation_rate).trace();
}

template <int DIM, int n_features, int degree>
double rans_twomodels<DIM, n_features, degree>::v_f3(const matrix_rans& strain_rate,
                                                     const matrix_rans& rotation_rate,
                                                     const vector_rans& nabla_pressure,
                                                     const matrix_rans& nabla_velocity)
{
    return (strain_rate * strain_rate * strain_rate).trace();
}

template <int DIM, int n_features, int degree>
double rans_twomodels<DIM, n_features, degree>::v_f4(const matrix_rans& strain_rate,
                                                     const matrix_rans& rotation_rate,
                                                     const vector_rans& nabla_pressure,
                                                     const matrix_rans& nabla_velocity)
{
    return (rotation_rate * rotation_rate * strain_rate).trace();
}

template <int DIM, int n_features, int degree>
double rans_twomodels<DIM, n_features, degree>::v_f5(const matrix_rans& strain_rate,
                                                     const matrix_rans& rotation_rate,
                                                     const vector_rans& nabla_pressure,
                                                     const matrix_rans& nabla_velocity)
{
    return (strain_rate * strain_rate * rotation_rate * rotation_rate).trace();
}

template <int DIM, int n_features, int degree>
double rans_twomodels<DIM, n_features, degree>::p_f8(const matrix_rans& strain_rate,
                                                     const matrix_rans& rotation_rate,
                                                     const vector_rans& nabla_pressure,
                                                     const matrix_rans& nabla_velocity)
{
    return nabla_pressure.norm() / (nabla_pressure.norm() + nabla_velocity.norm());
}

template <int DIM, int n_features, int degree>
double rans_twomodels<DIM, n_features, degree>::v_q_criterion(const matrix_rans& strain_rate,
                                                              const matrix_rans& rotation_rate,
                                                              const vector_rans& nabla_pressure,
                                                              const matrix_rans& nabla_velocity)
{
    auto rr = abs((rotation_rate * rotation_rate.transpose()).trace());
    auto ss = abs((strain_rate * strain_rate.transpose()).trace());
    return (rr - ss) / (rr + ss);
}

template <int DIM, int n_features, int degree>
double rans_twomodels<DIM, n_features, degree>::v_strain_rate_mag(const matrix_rans& strain_rate,
                                                                  const matrix_rans& rotation_rate,
                                                                  const vector_rans& nabla_pressure,
                                                                  const matrix_rans& nabla_velocity)
{
    return std::sqrt(2 * (strain_rate * strain_rate).trace());
}

template <int DIM, int n_features, int degree>
double rans_twomodels<DIM, n_features, degree>::v_rot_rate_mag(const matrix_rans& strain_rate,
                                                               const matrix_rans& rotation_rate,
                                                               const vector_rans& nabla_pressure,
                                                               const matrix_rans& nabla_velocity)
{
    return std::sqrt(2 * (rotation_rate * rotation_rate).trace());
}

}  // namespace Gascoigne

#endif
