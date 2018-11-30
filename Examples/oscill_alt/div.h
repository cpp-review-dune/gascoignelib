/*----------------------------   DIV_proj.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __DIV_proj_H
#define __DIV_proj_H
/*----------------------------   DIV_proj.h     ---------------------------*/

#include "equation.h"
#include "boundaryequation.h"
#include "paramfile.h"
#include "lpsequation.h"
#include <eigen3/Eigen/Dense>
#include "chi.h"

/*-----------------------------------------*/

namespace Gascoigne
{
template <int DIM>
class DIV_proj : public LpsEquation  // , public BoundaryEquation
{
protected:
    typedef Eigen::Matrix<double, DIM, DIM> MATRIX;
    typedef Eigen::Matrix<double, DIM, 1> VECTOR;

    // problem parameter
    double rho_f, nu_f, rho_s, lambda_s, mu_s, nu_e;
    double extend0;
    double pp_e;

    mutable double __h;
    mutable Vertex<DIM> __v;

    // stuff from point_M
    mutable double divergence;
    mutable MATRIX PRESSURE_P;
    mutable std::array<double, DIM> Jj, DIVERGENCE_U, DIVERGENCE_V;
    mutable MATRIX Fij[DIM];
    mutable MATRIX TENSOR_dU[DIM];
    mutable MATRIX TENSOR_dV[DIM];
    // mutable MATRIX TENSOR_dU_old[DIM];
    // mutable MATRIX TENSOR_dV_old[DIM];
    mutable MATRIX PRESSURE_U[DIM];

    /* // boundary */
    /* mutable VECTOR normal; */
    /* mutable VECTOR BOUNDARY, BOUNDARY_old; */

    // stuff from point
    mutable double lps, lps0, J, J_old;
    mutable MATRIX NV, NV_old, F, F_old;
    mutable MATRIX SIGMAf, SIGMAf_old;

    Chi chi;
    mutable int domain;

    mutable FemFunction* OLD;

    void SetFemData(FemData& q) const override final
    {
        assert(q.find("old") != q.end());
        OLD = &q["old"];
    }

public:
    ~DIV_proj()
    {
    }
    DIV_proj()
    {
        abort();
    }

    DIV_proj(const ParamFile* pf);

    std::string GetName() const override final
    {
        return "div";
    }

    int GetNcomp() const override final
    {
        return 2 * DIM + 1;
    }

    void point(double h, const FemFunction& U, const Vertex<DIM>& v) const override final;
    void point_M(int j, const FemFunction& U, const TestFunction& M) const;

    void Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const override final;

    void Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M,
                const TestFunction& N) const override final;
    void MatrixBlock(EntryMatrix& A, const FemFunction& U,
                     const FemFunction& NNN) const override final;

    /* ////////////////////////////////////////////////// Boundary */

    /* void Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const; */
    /* void Matrix(EntryMatrix& E, const FemFunction& U, const TestFunction& M, const TestFunction&
     * N, int col) const; */

    /* void pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v, const Vertex<DIM>&
     * n) const; */

    ////////////////////////////////////////////////// LPS

    void lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const override final;

    void StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP,
                  const TestFunction& N) const override final;

    void StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np,
                    const TestFunction& Mp) const override final;
};

}  // namespace Gascoigne

/*----------------------------   DIV_proj.h     ---------------------------*/
/* end of #ifndef __DIV_proj_H */
#endif
/*----------------------------   DIV_proj.h     ---------------------------*/
