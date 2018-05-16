#include "fsi.h"
#include "filescanner.h"

extern double __DT;
extern double __THETA;

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
template <int DIM>
FSI<DIM>::FSI(const ParamFile* pf)
{
    DataFormatHandler DFH;
    DFH.insert("rho_s", &rho_s, 0.0);
    DFH.insert("pp_e", &pp_e, 0.0);
    DFH.insert("lambda_s", &lambda_s, 0.0);
    DFH.insert("mu_s", &mu_s, 0.0);
    DFH.insert("rho_f", &rho_f, 0.0);
    DFH.insert("nu_f", &nu_f, 0.0);
    DFH.insert("extend", &extend0, 0.0);
    DFH.insert("lps", &lps0, 0.0);
    DFH.insert("theta", &theta, 0.51);
    DFH.insert("nu_e", &nu_e, 0.0);

    FileScanner FS(DFH, pf, "Equation");
    assert(rho_f > 0);
    assert(nu_f > 0);
    assert(extend0 > 0);
    assert(lps0 > 0);
    assert(rho_s > 0);
    assert(lambda_s > 0);
    assert(mu_s > 0);
    assert(theta > 0);

    cout << "%%%%%%%%%% Problem %%%%%%%%%%" << endl
         << "  rho_f / nu_f: " << rho_f << " / " << nu_f << endl
         << "  rho_s / mu_s / lambda_s: " << rho_s << " / " << mu_s << " / " << lambda_s << endl
         << "  extend / lps: " << extend0 << " / " << lps0 << endl;

    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl;
}

//////////////////////////////////////////////////

#include "multiplex.xx"

template <int DIM>
void FSI<DIM>::point(double h, const FemFunction& U, const Vertex<DIM>& v) const
{
    __h    = h;
    __v    = v;
    domain = chi(v);

    double dx   = std::max(0.0, fabs(v.x() - 0.4) - 0.4);
    double dy   = std::max(0.0, fabs(v.y() - 0.2) - 0.01);
    double dist = sqrt(dx * dx + dy * dy);

    extend = extend0 / (1.e-2 + dist);

    multiplex_init_NU<DIM>(NU, U);
    multiplex_init_NU<DIM>(NU_old, *OLD);
    // set F, F_old to Identity in fluid
    if (s)
    {
        multiplex_init_F<DIM>(F, U);
        multiplex_init_F<DIM>(F_old, *OLD);

        J     = F.determinant();
        J_old = F_old.determinant();
    }
    multiplex_init_dtU<DIM>(dtU, U, *OLD);
    // dtU = VECTOR::Zero();
    if (domain < 0)
    {
        if (!s)
        {
            F     = MATRIX::Identity();
            F_old = MATRIX::Identity();

            J     = 1.;
            J_old = 1.;
        }
        multiplex_init_NV<DIM>(NV, U);

        multiplex_init_V<DIM>(V, U);

        multiplex_init_NV<DIM>(NV_old, *OLD);
        multiplex_init_V<DIM>(V_old, *OLD);

        // TENSOR
        SIGMAf = rho_f * nu_f * (F.inverse() * NV + NV.transpose() * F.inverse().transpose());
        SIGMAf_old =
          rho_f * nu_f
          * (F_old.inverse() * NV_old + NV_old.transpose() * F_old.inverse().transpose());
    }
    if (domain > 0)
    {
        if (!s)
        {
            multiplex_init_F<DIM>(F, U);
            multiplex_init_F<DIM>(F_old, *OLD);

            J     = F.determinant();
            J_old = F_old.determinant();
        }
        E            = 0.5 * (F.transpose() * F - MATRIX::Identity());
        MATRIX E_old = 0.5 * (F_old.transpose() * F_old - MATRIX::Identity());

        SIGMAs     = (2.0 * mu_s * E + lambda_s * E.trace() * MATRIX::Identity());
        SIGMAs_old = (2.0 * mu_s * E_old + lambda_s * E_old.trace() * MATRIX::Identity());
    }
}

template <int DIM>
void FSI<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N) const
{
    VECTOR phi;
    multiplex_init_test<DIM>(phi, N);

    if (domain < 0)  // fluid
    {
        // divergence
        b[0] += rho_f * J * (F.inverse().transpose().array() * NV.array()).sum() * N.m();

        // time-derivative
        for (int i = 0; i < DIM; ++i)
            b[i + 1] += J * rho_f * (U[i + 1].m() - (*OLD)[i + 1].m()) / GetTimeStep() * N.m();

        // tensor
        VECTOR X = J * SIGMAf * F.inverse().transpose() * phi;
        for (int i = 0; i < DIM; ++i)
            b[i + 1] += theta * X(i, 0);

        X = J_old * SIGMAf_old * F_old.inverse().transpose() * phi;
        for (int i = 0; i < DIM; ++i)
            b[i + 1] += (1.0 - theta) * X(i, 0);

        // convection
        X = theta * J * rho_f * NV * F.inverse() * V * N.m();
        for (int i = 0; i < DIM; ++i)
            b[i + 1] += X(i, 0);

        X = (1.0 - theta) * J_old * rho_f * NV_old * F_old.inverse() * V_old * N.m();
        for (int i = 0; i < DIM; ++i)
            b[i + 1] += X(i, 0);

        // DOMAIN Convection
        X = -rho_f
            * (theta * J * NV * F.inverse() + (1.0 - theta) * J_old * NV_old * F_old.inverse())
            * dtU / GetTimeStep() * N.m();
        for (int i = 0; i < DIM; ++i)
            b[i + 1] += X(i, 0);

        // pressure
        X = -U[0].m() * J * F.inverse().transpose() * phi;
        for (int i = 0; i < DIM; ++i)
            b[i + 1] += X(i, 0);

        // // extension of deformation u with pseudo-elasticity
        // for (int i=0;i<DIM;++i)
        //   {
        //     double mu_e     = extend / (2.0 * (1.0+nu_e));
        //     double lambda_e = nu_e * extend / ( (1.0+nu_e)*(1.0-2*nu_e) );
        //     for (int j=0;j<DIM;++j)
        //       b[i+1+DIM] += mu_e * (U[i+1+DIM][j+1]+U[j+1+DIM][i+1])*N[j+1];
        //     for (int j=0;j<DIM;++j)
        //       b[i+1+DIM] += lambda_e * U[j+1+DIM][j+1]*N[j+1];
        //   }
        // laplace
        for (int i = 0; i < DIM; ++i)
            b[i + 1 + DIM] += extend * (NU * phi)(i, 0);

        // ////////// p-laplace
        // double dp=0;
        // for (int i=0;i<DIM;++i)
        //   for (int j=0;j<DIM;++j)
        //     dp += pow(U[i+1+DIM][j+1],2.0);

        // for (int i=0;i<DIM;++i)
        //   for (int j=0;j<DIM;++j)
        //     b[i+1+DIM] += pow(1.e-1 + dp,pp_e) * U[i+1+DIM][j+1]*N[j+1];
    }
    else
    {
        // time-derivative
        for (int i = 0; i < DIM; ++i)
            b[i + 1] += rho_s * (U[i + 1].m() - (*OLD)[i + 1].m()) / GetTimeStep() * N.m();

        // FULL tensor F Sigma

        for (int i = 0; i < DIM; ++i)
            b[i + 1] += theta * (F * SIGMAs * phi)(i, 0);
        for (int i = 0; i < DIM; ++i)
            b[i + 1] += (1.0 - theta) * (F_old * SIGMAs_old * phi)(i, 0);

        // dt u-v=0
        double scaling = rho_s / GetTimeStep() / theta;
        scaling        = 1.0;
        for (int i = 0; i < DIM; ++i)
        {
            b[i + 1 + DIM] +=
              scaling * (U[i + 1 + DIM].m() - (*OLD)[i + 1 + DIM].m()) / GetTimeStep() * N.m();
            b[i + 1 + DIM] -= scaling * theta * U[i + 1].m() * N.m();
            b[i + 1 + DIM] -= scaling * (1.0 - theta) * (*OLD)[i + 1].m() * N.m();
        }
    }
}

template <int DIM>
void FSI<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M,
                      const TestFunction& N) const
{
    VECTOR phi;
    multiplex_init_test<DIM>(phi, N);
    VECTOR psi;
    multiplex_init_test<DIM>(psi, M);

    if (domain < 0)  // fluid
    {
        //////////////// divergence
        //	b[0] += rho_f * J * (F.inverse().transpose().array() * NV.array()).sum() * N.m();
        //  wrt v
        for (int j = 0; j < DIM; ++j)
        {
            // wrt v
            A(0, j + 1) += DIVERGENCE_V[j] * N.m();

            // wrt u
            A(0, j + 1 + DIM) += rho_f * Jj[j] * divergence * N.m();

            A(0, j + 1 + DIM) += DIVERGENCE_U[j] * N.m();
        }

        ///////////// time-derivative
        // wrt v
        for (int i = 0; i < DIM; ++i)
            A(i + 1, i + 1) += rho_f * J * M.m() * N.m() / GetTimeStep();
        // wrt u
        for (int j = 0; j < DIM; ++j)
            for (int i = 0; i < DIM; ++i)
                A(i + 1, j + 1 + DIM) +=
                  Jj[j] * rho_f * (U[i + 1].m() - (*OLD)[i + 1].m()) / GetTimeStep() * N.m();

        ///////// tensor
        // wrt V
        for (int j = 0; j < DIM; ++j)
        {
            VECTOR X = TENSOR_dV[j] * phi;
            for (int i = 0; i < DIM; ++i)
                A(i + 1, j + 1) += theta * X(i, 0);
        }
        // wrt U
        for (int j = 0; j < DIM; ++j)
        {
            VECTOR X = TENSOR_dU[j] * phi;
            for (int i = 0; i < DIM; ++i)
                A(i + 1, j + 1 + DIM) += theta * X(i, 0);
        }

        // //////////////// // convection
        // X = theta * rho_f * NV * F.inverse()*V*N.m();
        // for (int i=0;i<DIM;++i)
        //   b[i+1] += X(i,0);

        // wrt v
        for (int i = 0; i < DIM; ++i)
            for (int j = 0; j < DIM; ++j)
                A(i + 1, j + 1) += theta * CONV_dV1(i, j) * M.m() * N.m();
        for (int j = 0; j < DIM; ++j)
            A(j + 1, j + 1) += theta * CONV_dV2[j] * N.m();

        // wrt u
        for (int j = 0; j < DIM; ++j)
        {
            for (int i = 0; i < DIM; ++i)
                A(i + 1, j + 1 + DIM) += theta * CONV_dU[j](i, 0) * N.m();
        }

        // DOMAIN Convection
        //	X = -rho_f * 0.5 * (J*NV*F.inverse() + J_old * NV_old * F_old.inverse())
        //	  * dtU/__DT * N.m();

        // wrt dtU
        for (int i = 0; i < DIM; ++i)
            for (int j = 0; j < DIM; ++j)
                A(i + 1, j + 1 + DIM) += DOMAIN_U1(i, j) * M.m() * N.m();

        // wrt V
        for (int j = 0; j < DIM; ++j)
            A(j + 1, j + 1) += DOMAIN_V[j] * N.m();
        // wrt U
        for (int j = 0; j < DIM; ++j)
        {
            for (int i = 0; i < DIM; ++i)
                A(i + 1, j + 1 + DIM) += Jj[j] * DOMAIN_U2(i, 0) * N.m();

            VECTOR Y = -theta * rho_f * J * NV * Fij[j] * dtU / GetTimeStep() * N.m();
            for (int i = 0; i < DIM; ++i)
                A(i + 1, j + 1 + DIM) += Y(i, 0);
        }

        /////////////// pressure
        // X = -U[0].m() * J * F.inverse().transpose()*phi;
        // for (int i=0;i<DIM;++i)
        //   b[i+1] += X(j,0);
        // wrt P
        VECTOR X = PRESSURE_P * phi;
        for (int i = 0; i < DIM; ++i)
            A(i + 1, 0) += X(i, 0);
        // wrt U
        for (int j = 0; j < DIM; ++j)
        {
            X = PRESSURE_U[j] * phi;
            for (int i = 0; i < DIM; ++i)
                A(i + 1, j + 1 + DIM) += X(i, 0);
        }

        ////// EXTEND

        // // extension of deformation u with pseudo-elasticity
        // for (int i=0;i<DIM;++i)
        //   {
        //     double mu_e     = extend / (2.0 * (1.0+nu_e));
        //     double lambda_e = nu_e * extend / ( (1.0+nu_e)*(1.0-2*nu_e) );
        //     for (int j=0;j<DIM;++j)
        //       {
        // 	A(i+1+DIM,i+1+DIM) += mu_e * (M[j+1]+U[j+1+DIM][i+1])*N[j+1];
        // 	A(i+1+DIM,j+1+DIM) += mu_e * (M[i+1])*N[j+1];
        //       }

        //     for (int j=0;j<DIM;++j)
        //       A(i+1+DIM,j+1+DIM) += lambda_e * M[j+1]*N[j+1];
        //   }
        // // extension of deformation u laplace
        for (int i = 0; i < DIM; ++i)
            A(i + 1 + DIM, i + 1 + DIM) += extend * phi.dot(psi);

        // // p-laplace
        // double dp=0;
        // for (int i=0;i<DIM;++i)
        //   for (int j=0;j<DIM;++j)
        //     dp += pow(U[i+1+DIM][j+1],2.0);
        // fixarray<DIM,double> dpD;
        // for (int i=0;i<DIM;++i)
        //   {
        //     dpD[i]=0;
        //     for (int j=0;j<DIM;++j)
        //       dpD[i] += 2.0 * U[i+1+DIM][j+1]*M[j+1];
        //   }
        // for (int i=0;i<DIM;++i)
        //   for (int j=0;j<DIM;++j)
        //     {
        //       for (int k=0;k<DIM;++k)
        // 	A(i+1+DIM,k+1+DIM) += pp_e * pow(1.e-1 + dp,pp_e-1.0) * dpD[k] * U[i+1+DIM][j+1]*N[j+1];
        //       A(i+1+DIM,i+1+DIM)  += pow(1.e-1 + dp,pp_e) * M[j+1]*N[j+1];
        //     }
    }
    else
    {
        // time-derivative
        for (int i = 0; i < DIM; ++i)
            A(i + 1, i + 1) += rho_s * M.m() / GetTimeStep() * N.m();

        // F Sigma
        // wrt F,
        for (int j = 0; j < DIM; ++j)  // Fj Sigma \nabla phi
            A(j + 1, j + 1 + DIM) += theta * (SIGMA_dF[j].transpose() * phi)(0, 0);

        for (int j = 0; j < DIM; ++j)  // F Sigma_j \nabla phi
        {
            VECTOR X = SIGMA_dU[j] * phi;
            for (int i = 0; i < DIM; ++i)
                A(i + 1, j + 1 + DIM) += theta * X(i, 0);
        }

        // dt u-v=0
        double scaling = rho_s / GetTimeStep() / theta;
        scaling        = 1.0;

        for (int i = 0; i < DIM; ++i)
        {
            A(i + 1 + DIM, i + 1 + DIM) += scaling * M.m() * N.m() / GetTimeStep();
            A(i + 1 + DIM, i + 1) -= scaling * theta * M.m() * N.m();
        }
    }
}

template <int DIM>
void FSI<DIM>::point_M(int j, const FemFunction& U, const TestFunction& M) const
{
    VECTOR psi;
    multiplex_init_test<DIM>(psi, M);
    if (s)
    {
        for (int j = 0; j < DIM; ++j)
        {
            // set to zero
            Jj[j]  = (psi.transpose() * J * F.inverse().block(0, j, DIM, 1))(0, 0);
            Fij[j] = -F.inverse().block(0, j, DIM, 1) * psi.transpose() * F.inverse();
        }
    }
    if (domain < 0)
    {
        if (!s)
        {
            for (int j = 0; j < DIM; ++j)
            {
                Jj[j]  = 0.;
                Fij[j] = MATRIX::Zero();
            }
        }
        divergence = (F.inverse().transpose().array() * NV.array()).sum();
        CONV_dV1   = rho_f * J * NV * F.inverse();

        DOMAIN_U1 =
          -rho_f * (theta * J * NV * F.inverse() + (1.0 - theta) * J_old * NV_old * F_old.inverse())
          / GetTimeStep();
        DOMAIN_U2 = -theta * rho_f * NV * F.inverse() * dtU / GetTimeStep();

        for (int j = 0; j < DIM; ++j)
        {
            TENSOR_dV[j] = rho_f * nu_f * J
                           * (F.inverse().block(0, j, DIM, 1) * psi.transpose()
                              + (F.inverse().block(0, j, DIM, 1) * psi.transpose()).transpose())
                           * F.inverse().transpose();
            TENSOR_dU[j] = rho_f * nu_f * J * (Fij[j] * NV + (Fij[j] * NV).transpose())
                             * F.inverse().transpose()                  // p1
                           + J * SIGMAf * Fij[j].transpose()            // p2
                           + Jj[j] * SIGMAf * F.inverse().transpose();  // p3

            DIVERGENCE_U[j] = rho_f * J * (Fij[j].transpose().array() * NV.array()).sum();
            DIVERGENCE_V[j] = rho_f * J * (psi.transpose() * F.inverse().block(0, j, DIM, 1))(0, 0);

            CONV_dV2[j] = J * rho_f * (psi.transpose() * F.inverse() * V)(0, 0);
            CONV_dU[j]  = rho_f * NV * F.inverse() * V * Jj[j] + rho_f * J * NV * Fij[j] * V;

            DOMAIN_V[j] =
              -theta * rho_f * J * (psi.transpose() * F.inverse() * dtU / GetTimeStep())(0, 0);

            PRESSURE_P = -M.m() * J * F.inverse().transpose();

            PRESSURE_U[j] =
              -U[0].m() * Jj[j] * F.inverse().transpose() - U[0].m() * J * Fij[j].transpose();
        }
    }

    if (domain > 0)
    {
        if (!s)
        {
            for (int j = 0; j < DIM; ++j)
            {
                Jj[j]  = (psi.transpose() * J * F.inverse().block(0, j, DIM, 1))(0, 0);
                Fij[j] = -F.inverse().block(0, j, DIM, 1) * psi.transpose() * F.inverse();
            }
        }
        for (int j = 0; j < DIM; ++j)
            SIGMA_dF[j] = (psi.transpose() * SIGMAs).transpose();
        for (int j = 0; j < DIM; ++j)
        {
            MATRIX Ej = 0.5
                        * (psi * (F.block(j, 0, 1, DIM))
                           + F.block(j, 0, 1, DIM).transpose() * psi.transpose());
            SIGMA_dU[j] = F * (2.0 * mu_s * Ej + lambda_s * Ej.trace() * MATRIX::Identity());
        }
    }
}

template <int DIM>
void FSI<DIM>::MatrixBlock(EntryMatrix& A, const FemFunction& U, const FemFunction& N) const
{
    for (int j = 0; j < N.size(); ++j)  // trial
    {
#define M N[j]
        point_M(j, U, M);

        for (int i = 0; i < N.size(); ++i)  // test
        {
            A.SetDofIndex(i, j);
            Matrix(A, U, M, N[i]);
        }
#undef M
    }
}

/*
template<int DIM>
void FSI<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const
TestFunction& N) const
{
  DoubleVector B(GetNcomp());

  B.zero();


  point(__h,U,__v);

  Form(B.begin(), U, N);


  double EPS = 1.e-9;

  assert(U.size()==GetNcomp());


  for (int j=0;j<GetNcomp();++j)
    {

  DoubleVector Bj(GetNcomp());

  Bj.zero();


  FemFunction Uj(GetNcomp());

  for (int i=0; i<GetNcomp(); ++i)
    Uj[i] = U[i];

  Uj[j].add(EPS, M);




  point(__h,Uj,__v);

  Form(Bj.begin(), Uj, N);


  for (int i=0;i<GetNcomp();++i)
    A(i,j)+= (Bj[i]-B[i])/EPS;

    }

  return;


}
*/

////////////////////////////////////////////////// BOUNDARY

// template<int DIM>
// void FSI<DIM>::Form(VectorIterator b, const FemFunction& U, const TestFunction& N, int col) const
// {
//   if (domain>0) return;

//   if (col==1)
//     {
// 	for (int i=0;i<DIM;++i)
// 	  b[i+1] -= theta       * BOUNDARY(i,0)*N.m();
// 	for (int i=0;i<DIM;++i)
// 	  b[i+1] -= (1.0-theta) * BOUNDARY_old(i,0)*N.m();
//     }

// }

// template<int DIM>
// void FSI<DIM>::Matrix(EntryMatrix& A, const FemFunction& U, const TestFunction& M, const
// TestFunction& N, int col) const
// {
//   if (domain>0) return;

//   if (col==1)
//     {
//   	VECTOR psi; multiplex_init_test<DIM>(psi,M);
//   	for (int j=0;j<DIM;++j)
//   	  {
// 	    MATRIX FIJ = -F.inverse().block(0,j,DIM,1)*psi.transpose()*F.inverse();
// 	    double JJ =  (psi.transpose() * J  * F.inverse().block(0,j,DIM,1))(0,0);

// 	    VECTOR BOUNDARY_U = rho_f * nu_f * N.m() *
// 	      (  JJ * NV.transpose()*F.inverse().transpose()*F.inverse().transpose() +
// 		 JJ * NV.transpose() * (FIJ.transpose()*F.inverse().transpose() +
// F.inverse().transpose()*FIJ.transpose()) )  * normal; 	    double BOUNDARY_V = rho_f * nu_f *
// N.m() * J * 	      (psi.transpose() *
// (F.inverse().transpose()*F.inverse().transpose())*normal)(0,0); 	    for (int i=0;i<DIM;++i)
// 	      {
// 		A(i+1,j+1+DIM) -= theta * BOUNDARY_U(i,0);
// 		A(i+1,j+1)     -= theta * BOUNDARY_V;
// 	      }
//   	  }

//     }

// }

// template<int DIM>
// void FSI<DIM>::pointboundary(double h, const FemFunction& U, const Vertex<DIM>& v, const
// Vertex<DIM>& n) const
// {
//   domain = chi(v);
//   if (domain>0) return; // no boundary eq in solid

//   multiplex_init_F<DIM>(F,U);
//   multiplex_init_F<DIM>(F_old,*OLD);

//   multiplex_init_NV<DIM>(NV,U);
//   multiplex_init_NV<DIM>(NV_old,*OLD);

//   J = F.determinant();
//   J_old = F_old.determinant();

//   multiplex_init_normal<DIM>(normal,n);

//   BOUNDARY = rho_f * nu_f * J *
//   NV.transpose()*F.inverse().transpose()*F.inverse().transpose()*normal; BOUNDARY_old = rho_f *
//   nu_f * J_old *
//   NV_old.transpose()*F_old.inverse().transpose()*F_old.inverse().transpose()*normal;

// }

////////////////////////////////////////////////// LPS

template <int DIM>
void FSI<DIM>::lpspoint(double h, const FemFunction& U, const Vertex<DIM>& v) const
{
    double vel = 1.0;

    lps    = lps0 / (vel / h + nu_f / h / h);
    domain = chi(v);
}
template <int DIM>
void FSI<DIM>::StabForm(VectorIterator b, const FemFunction& U, const FemFunction& UP,
                        const TestFunction& N) const
{
    if (domain < 0)  /// fludi
        for (int i = 0; i < DIM; ++i)
            b[0] += lps * UP[0][i + 1] * N[i + 1];
}

template <int DIM>
void FSI<DIM>::StabMatrix(EntryMatrix& A, const FemFunction& U, const TestFunction& Np,
                          const TestFunction& Mp) const
{
    if (domain < 0)
        for (int i = 0; i < DIM; ++i)
            A(0, 0) += lps * Mp[i + 1] * Np[i + 1];
}

template class FSI<2>;
template class FSI<3>;

}  // namespace Gascoigne

/*-----------------------------------------*/
