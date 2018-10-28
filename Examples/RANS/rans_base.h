#ifndef RANS_BASE_H
#define RANS_BASE_H

#include <cmath>
#include <eigen3/Eigen/Dense>
#include "equation.h"
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
class rans_base : public LpsEquation
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

    ~rans_base();

protected:
    using matrix_rans              = Eigen::Matrix<double, DIM, DIM>;
    using vector_rans              = Eigen::Matrix<double, DIM, 1>;
    static constexpr auto rst_size = static_cast<size_t>(DIM * (DIM + 1) / 2);
    mutable FemFunction* _old;
    mutable vector_rans _phi;
    mutable double _h;
    mutable matrix_rans _rst;

private:
    mutable double _visc;
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
    DFH.insert("theta", &_theta, 0.5);
    DFH.insert("alpha0", &_alpha0, 0.5);
    DFH.insert("delta0", &_delta0, 0.5);
    FileScanner FS(DFH, pf, "Equation");
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
    double _norm = std::sqrt(U[1].m() * U[1].m() + U[2].m() * U[2].m());
    double val   = 6. * _visc / (h * h) + _norm / h;
    _delta       = _delta0 / val;
    _alpha       = _alpha0 / val;
}

/*-----------------------------------------*/

template <int DIM>
void rans_base<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
    // Divergence
    b[0] += (U[1].x() + U[2].y()) * N.m();

    Multiplex::init_test<DIM>(_phi, N);
    auto x = _rst * _phi;
    for (auto i = 1; i <= DIM; i++)
    {
        // Time derivative
        b[i] += (U[i].m() - (*_old)[i].m()) * N.m() / (GetTimeStep());
        // Pressure
        b[i] -= U[0].m() * N[i];
        // Laplace
        b[i] += _theta * _visc * (U[i].x() * N.x() + U[i].y() * N.y());
        b[i] += (1 - _theta) * _visc * ((*_old)[i].x() * N.x() + (*_old)[i].y() * N.y());
        // Convection
        b[i] += _theta * (U[1].m() * U[i].x() + U[2].m() * U[i].y()) * N.m();
        b[i] += (1 - _theta) * ((*_old)[1].m() * (*_old)[i].x() + (*_old)[2].m() * (*_old)[i].y())
                * N.m();
        b[i] += x[i - 1];
    }
}

/*-----------------------------------------*/

template <int DIM>
void rans_base<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M,
                            const TestFunction& N) const
{
    // Laplace
    double lp_conv = _visc * (M.x() * N.x() + M.y() * N.y());
    // Convection
    for (auto i = 1; i <= DIM; i++)
    {
        lp_conv += U[i].m() * M[i];
    }
    lp_conv *= (_theta * N.m());
    for (auto i = 1; i <= DIM; i++)
    {
        // Time derivative
        A(i, i) += M.m() * N.m() / GetTimeStep();
        // Laplace and Convection
        A(i, i) += lp_conv;
        // Divergence
        A(0, i) += M[i] * N.m();
        // Pressure
        A(i, 0) -= N[i] * M.m();
    }
}

/*-----------------------------------------*/

template <int DIM>
void rans_base<DIM>::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP,
                              const TestFunction& N) const
{
    auto beta = U[1].m() * N.x() + U[2].m() * N.y();
    for (auto i = 1; i <= DIM; i++)
    {
        b[0] += _alpha * UP[0][i] * N[i];
        b[1] += _delta * beta * (U[1].m() * UP[i].x() + U[2].m() * UP[i].y());
        b[2] += _delta * beta * (U[1].m() * UP[i].x() + U[2].m() * UP[i].y());
    }
}

/*-----------------------------------------*/

template <int DIM>
void rans_base<DIM>::StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np,
                                const TestFunction& Mp) const
{
    A(0, 0) += _alpha * (Mp.x() * Np.x() + Mp.y() * Np.y());
    auto beta_stab =
      _delta * (U[1].m() * Mp.x() + U[2].m() * Mp.y()) * (U[1].m() * Np.x() + U[2].m() * Np.y());
    for (auto i = 1; i <= DIM; i++)
        A(i, i) += beta_stab;
}
}  // namespace Gascoigne

#endif
