#ifndef RANS_BASE_H
#define RANS_BASE_H

#include <cmath>
#include <eigen3/Eigen/Dense>
#include "equation.h"
#include "boundaryequation.h"
#include "lpsequation.h"
#include "paramfile.h"
#include "filescanner.h"
#include "dataformathandler.h"
#include "multiplex.h"

/*-----------------------------------------*/

namespace Gascoigne
{
//! \brief Basic RANS implementation with LPS stabilization.
//!
//! No calculation of RST, it is zero by default. So basically this is Navier Stokes with
//! LPS-stabilization when RST does not get changed
template <int DIM>
class rans_base : public LpsEquation  // , public BoundaryEquation
{
public:
    rans_base();
    rans_base(const ParamFile* pf);

    std::string GetName() const
    {
        return "rans";
    }

    int GetNcomp() const
    {
        return DIM + 1;
    }

    void SetFemData(FemData& q) const override final
    {
        assert(q.find("old") != q.end());
        _old = &q["old"];
    }

    void point(double h, const FemFunction& U, const Vertex2d& v) const;

    void lpspoint(double h, const FemFunction& U, const Vertex2d& v) const;

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const override final;

    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M,
                const TestFunction& N) const override final;

    void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP,
                  const TestFunction& N) const override final;

    void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np,
                    const TestFunction& Mp) const override final;
    // BOUNDARY
    // void Form(VectorIterator b, const FemFunction& U, const TestFunction& N,
    //           int col) const override final;
    // void Matrix(EntryMatrix& E, const FemFunction& U, const TestFunction& M, const TestFunction&
    // N,
    //             int col) const override final;

    void pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v,
                       const Vertex<DIM>& n) const
    {
        //_n     = n;
        _vdotn = 0;
        for (auto i = 1; i <= DIM; i++)
        {
            _vdotn += U[i][0] * n[i];
        }
        _vdotn = (_vdotn < 0) ? _vdotn : 0;
    }

    ~rans_base();

protected:
    using matrix_rans              = Eigen::Matrix<double, DIM, DIM>;
    using vector_rans              = Eigen::Matrix<double, DIM, 1>;
    static constexpr auto rst_size = static_cast<size_t>(DIM * (DIM + 1) / 2);
    mutable FemFunction* _old;
    mutable vector_rans _phi;
    mutable double _h;
    mutable matrix_rans _rst;
    mutable double _vdotn;

private:
    mutable double _visc;
    mutable double _visc0;
    mutable double _theta;

    mutable double _alpha0;
    mutable double _alpha;
    mutable double _delta0;
    mutable double _delta;
};

// Implementation

template <int DIM>
rans_base<DIM>::~rans_base()
{
}

/*-----------------------------------------*/

template <int DIM>
rans_base<DIM>::rans_base()
    : LpsEquation()
    , _h(0.0)
    , _rst(matrix_rans::Zero())
    , _visc(0.01)
    , _theta(0.5)
    , _alpha0(0.5)
    , _alpha(0.5)
    , _delta0(0.5)
    , _delta(0.0)
{
}

/*-----------------------------------------*/

template <int DIM>
rans_base<DIM>::rans_base(const ParamFile* pf)
    : LpsEquation(), _h(0.0), _rst(matrix_rans::Zero()), _alpha(0.0), _delta(0.0)
{
    DataFormatHandler DFH;
    DFH.insert("visc", &_visc, 1.);
    DFH.insert("theta", &_theta, 0.51);
    DFH.insert("alpha0", &_alpha0, 0.5);
    DFH.insert("delta0", &_delta0, 0.5);
    FileScanner FS(DFH, pf, "Equation");
    _visc0 = _visc;
}

/*-----------------------------------------*/

template <int DIM>
void rans_base<DIM>::point(double h, const FemFunction& U, const Vertex2d& v) const
{
    _h = h;
}

/*-----------------------------------------*/

template <int DIM>
void rans_base<DIM>::lpspoint(double h, const FemFunction& U, const Vertex2d& v) const
{
    double _norm = std::sqrt(U[1][0] * U[1][0] + U[2][0] * U[2][0]);
    double val   = 6. * _visc / (h * h) + _norm / h;
    _delta       = _delta0 / val;
    _alpha       = _alpha0 / val;
}

/*-----------------------------------------*/

template <int DIM>
void rans_base<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
    // [0] -> m()
    // [1] -> x()
    // [2] -> y()
    // [3] -> z()
    // Divergence
    for (auto i = 1; i < DIM + 1; ++i)
    {
        b[0] += U[i][i] * N[0];
    }

    Multiplex::init_test<DIM>(_phi, N);
    auto x = _rst * _phi;
    for (auto i = 1; i <= DIM; ++i)
    {
        // Time derivative
        b[i] += (U[i][0] - (*_old)[i][0]) * N[0] / (GetTimeStep());
        // Pressure
        b[i] -= U[0][0] * N[i];
        // Laplace
        for (auto j = 1; j <= DIM + 1; ++j)
        {
            b[i] += _theta * _visc * U[i][j] * N[j];
            b[i] += (1 - _theta) * _visc * (*_old)[i][j] * N[j];
        }
        // Convection
        for (auto j = 1; j <= DIM + 1; ++j)
        {
            b[i] += _theta * U[j][0] * U[i][j] * N[0];
            b[i] += (1 - _theta) * (*_old)[j][0] * (*_old)[i][j] * N[0];
        }
        b[i] += x[i - 1];
    }
}

/*-----------------------------------------*/

template <int DIM>
void rans_base<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M,
                            const TestFunction& N) const
{
    // Laplace
    double lp_conv = 0;
    for (auto i = 1; i <= DIM; ++i)
    {
        lp_conv += _visc * M[i] * N[i];
    }
    // Convection
    for (auto i = 1; i <= DIM; ++i)
    {
        lp_conv += U[i][0] * M[i];
    }
    lp_conv *= (_theta * N[0]);
    for (auto i = 1; i <= DIM; ++i)
    {
        // Time derivative
        A(i, i) += M[0] * N[0] / GetTimeStep();
        // Laplace and Convection
        A(i, i) += lp_conv;
        // Divergence
        A(0, i) += M[i] * N[0];
        // Pressure
        A(i, 0) -= N[i] * M[0];
    }
}

/*-----------------------------------------*/

template <int DIM>
void rans_base<DIM>::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP,
                              const TestFunction& N) const
{
    auto beta = 0.0;
    for (auto i = 1; i <= DIM; ++i)
    {
        beta += U[i][0] * N[i];
    }
    for (auto i = 1; i <= DIM; ++i)
    {
        b[0] += _alpha * UP[0][i] * N[i];
        auto conv = 0.0;
        for (auto j = 1; j <= DIM; ++j)
        {
            conv += U[j][0] * UP[i][j];
        }
        b[i] += _delta * beta * conv;
    }
}

/*-----------------------------------------*/

template <int DIM>
void rans_base<DIM>::StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np,
                                const TestFunction& Mp) const
{
    auto beta_stab  = 0.0;
    auto beta_stab2 = 0.0;
    for (auto j = 1; j <= DIM; ++j)
    {
        beta_stab += U[j][0] * Np[j];
        beta_stab2 += U[j][0] * Mp[j];
    }
    for (auto i = 1; i <= DIM; i++)
    {
        A(0, 0) += _alpha * Mp[i] * Np[i];
        A(i, i) += _delta * beta_stab * beta_stab2;
    }
}
}  // namespace Gascoigne

#endif
