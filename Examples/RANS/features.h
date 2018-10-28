namespace Gascoigne::feat
{
//! \brief Functions which calculate flow features from velocity gradient or pressure
//! gradient in a flow.
//!
//! Prefix velocity or v for a feature calculated from the velocity gradient
//! Prefix pressure or p for a feature calculated from the pressure gradient
'f1', 'f2', 'f3', 'f4', 'f5', 'f7', 'f8', 'f9', 'vorticity',
                'Q-criterion'
double v_f1
double v_f2
double v_f3
double v_f4
double v_f5
double v_f7
double p_f8
double p_f9
double v_vorticity
double v_q_criterion
}  // namespace Gascoigne::feat

    _feature_vector[0] = (std::abs(std::pow(r12, 2) + std::pow(r21, 2))
                          - std::abs(std::pow(s11, 2) + 2 * std::pow(s12, 2) + std::pow(s22, 2)))
                         / std::abs(std::pow(r12, 2) + std::pow(r21, 2) + std::pow(s11, 2)
                                    + 2 * std::pow(s12, 2) + std::pow(s22, 2));
    // Rrtm
    _feature_vector[1] = std::sqrt(_feature_vector[1]);
    // Srtm
    _feature_vector[2] = std::sqrt(_feature_vector[0]);
    // f1
    _feature_vector[3] = std::pow(s11, 2) + 2 * std::pow(s12, 2) + std::pow(s22, 2);
    // f2
    _feature_vector[4] = 2 * r12 * r21;
    // f3
    _feature_vector[5] =
      std::pow(s11, 3) + 3 * s11 * std::pow(s12, 2) + 3 * std::pow(s12, 2) * s22 + std::pow(s22, 3);
    // f4
    _feature_vector[6] = r12 * r21 * (s11 + s22);
    // f5
    _feature_vector[7] = r12 * r21 * (std::pow(s11, 2) + 2 * std::pow(s12, 2) + std::pow(s22, 2));
    // f6
    _feature_vector[8] =
      r12 * r21 * s12
      * (r12 * s11 * (s11 + s22) + r12 * (std::pow(s12, 2) + std::pow(s22, 2))
         + r21 * s22 * (s11 + s22) + r21 * (std::pow(s11, 2) + std::pow(s12, 2)));
    // f7
    _feature_vector[9] = r12 * r21
                         * (std::pow(s11, 3) + 3 * s11 * std::pow(s12, 2)
                            + 3 * std::pow(s12, 2) * s22 + std::pow(s22, 3));
    // f8 pressure magnitude
    _feature_vector[10] = (*rans_base<DIM>::_old)[0].x() * (*rans_base<DIM>::_old)[0].x()
                          + (*rans_base<DIM>::_old)[0].y() * (*rans_base<DIM>::_old)[0].y();
    // f9 pressure
    _feature_vector[11] = _feature_vector[7] - (*_old)[0].x() * (*_old)[0].y();
    // vorticity
    _feature_vector[12] = r12 / std::abs(s12);
