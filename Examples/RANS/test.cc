#include <cmath>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>

template <typename m_type, size_t rows, size_t cols>
Eigen::Matrix<m_type, rows, cols> load_matrix(const std::string& file, const char delim = ',')
{
    std::ifstream input;
    input.open(file);
    std::string row_data;
    std::array<m_type, cols * rows> val;
    for (size_t i = 0; i < rows; i++)
    {
        std::getline(input, row_data);
        std::stringstream row_data_stream(row_data);
        std::string element;
        for (size_t j = 0; j < cols; j++)
        {
            std::getline(row_data_stream, element, delim);
            val[i * cols + j] = static_cast<m_type>(std::stod(element));
        }
    }
    return Eigen::Map<const Eigen::Matrix<m_type, rows, cols, Eigen::RowMajor>>(val.data(), rows,
                                                                                cols);
}

int main()
{
    Eigen::Matrix<double, 10, 1> _feature_vector;
    _feature_vector << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9;
    Eigen::Matrix<double, 10, 2> mean_var_normalizations =
      load_matrix<double, 2, 10>("var_mean_scales.csv").transpose();
    _feature_vector -= mean_var_normalizations.col(0);
    _feature_vector.cwiseProduct(mean_var_normalizations.col(1).cwiseInverse());
    return 0;
}
