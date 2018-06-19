namespace Gascoigne::Multiplex
{
template <int DIM>
void inline init_NV(Eigen::Matrix<double, DIM, DIM>& NV, const FemFunction& U)
{
    assert(0);
}

template <>
void inline init_NV<2>(Eigen::Matrix<double, 2, 2>& NV, const FemFunction& U)
{
    NV << U[1].x(), U[1].y(), U[2].x(), U[2].y();
}

template <>
void inline init_NV<3>(Eigen::Matrix<double, 3, 3>& NV, const FemFunction& U)
{
    NV << U[1].x(), U[1].y(), U[1].z(), U[2].x(), U[2].y(), U[2].z(), U[3].x(), U[3].y(), U[3].z();
}

template <int DIM>
void inline init_normal(Eigen::Matrix<double, DIM, 1>& normal, const Vertex<DIM>& n)
{
    assert(0);
}

template <>
void inline init_normal<2>(Eigen::Matrix<double, 2, 1>& normal, const Vertex<2>& n)
{
    normal << n.x(), n.y();
}
template <>
void inline init_normal<3>(Eigen::Matrix<double, 3, 1>& normal, const Vertex<3>& n)
{
    normal << n.x(), n.y(), n.z();
}

template <int DIM>
void inline init_NU(Eigen::Matrix<double, DIM, DIM>& NU, const FemFunction& U)
{
    assert(0);
}

template <>
void inline init_NU<2>(Eigen::Matrix<double, 2, 2>& NU, const FemFunction& U)
{
    NU << U[3].x(), U[3].y(), U[4].x(), U[4].y();
}

template <>
void inline init_NU<3>(Eigen::Matrix<double, 3, 3>& NU, const FemFunction& U)
{
    NU << U[4].x(), U[4].y(), U[4].z(), U[5].x(), U[5].y(), U[5].z(), U[6].x(), U[6].y(), U[6].z();
}

template <int DIM>
void inline init_F(Eigen::Matrix<double, DIM, DIM>& F, const FemFunction& U)
{
    assert(0);
}

template <>
void inline init_F<2>(Eigen::Matrix<double, 2, 2>& F, const FemFunction& U)
{
    F << 1.0 + U[3].x(), U[3].y(), U[4].x(), 1.0 + U[4].y();
}

template <>
void inline init_F<3>(Eigen::Matrix<double, 3, 3>& F, const FemFunction& U)
{
    F << 1.0 + U[4].x(), U[4].y(), U[4].z(), U[5].x(), 1.0 + U[5].y(), U[5].z(), U[6].x(), U[6].y(),
      1.0 + U[6].z();
}

template <int DIM>
void inline init_V(Eigen::Matrix<double, DIM, 1>& V, const FemFunction& U)
{
    assert(0);
}

template <>
void inline init_V<2>(Eigen::Matrix<double, 2, 1>& V, const FemFunction& U)
{
    V << U[1].m(), U[2].m();
}

template <>
void inline init_V<3>(Eigen::Matrix<double, 3, 1>& V, const FemFunction& U)
{
    V << U[1].m(), U[2].m(), U[3].m();
}

template <int DIM>
void inline init_dtU(Eigen::Matrix<double, DIM, 1>& dtU, const FemFunction& U,
                     const FemFunction& OLD)
{
    assert(0);
}

template <>
void inline init_dtU<2>(Eigen::Matrix<double, 2, 1>& dtU, const FemFunction& U,
                        const FemFunction& OLD)
{
    dtU << U[3].m() - OLD[3].m(), U[4].m() - OLD[4].m();
}
template <>
void inline init_dtU<3>(Eigen::Matrix<double, 3, 1>& dtU, const FemFunction& U,
                        const FemFunction& OLD)
{
    dtU << U[4].m() - OLD[4].m(), U[5].m() - OLD[5].m(), U[6].m() - OLD[6].m();
}

template <int DIM>
void inline init_test(Eigen::Matrix<double, DIM, 1>& phi, const TestFunction& N)
{
    assert(0);
}
template <>
void inline init_test<2>(Eigen::Matrix<double, 2, 1>& phi, const TestFunction& N)
{
    phi << N.x(), N.y();
}
template <>
void inline init_test<3>(Eigen::Matrix<double, 3, 1>& phi, const TestFunction& N)
{
    phi << N.x(), N.y(), N.z();
}
}  // namespace Gascoigne::Multiplex
