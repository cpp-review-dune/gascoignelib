#ifndef POLYNOMIAL_FEATURES
#define POLYNOMIAL_FEATURES

#include <utility>
#include <iostream>
#include <array>
#include <vector>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <type_traits>

namespace Gascoigne
{
//////////////////////////////////////////////////////////////////////////
// Definitions of free functions                                        //
//////////////////////////////////////////////////////////////////////////
template <size_t n_features>
constexpr inline size_t binom(size_t k);

template <size_t n_features>
constexpr inline size_t summed_binom(size_t k);

//
//! \brief Returns compile-time array initialized with index dependent values. Function is given
//! as reference (e.g. binom, summed_binom).
template <size_t... i, typename fct>
constexpr inline auto f_array(const fct& f, std::index_sequence<i...>);

template <size_t size, typename fct>
constexpr inline auto f_array(const fct& f);

//
//! \brief summing up entries in an array. The STL does not provide constexpr versions of sum, so we
//! have to implement this ourselves
template <size_t size, typename T>
constexpr inline size_t sum(const std::array<T, size>& arr);

//
//! \brief Summing up entries in an array weighted by their index. Simple utility function for
//! calculation of array sizes.
template <size_t size>
constexpr inline size_t sum_weighted(const std::array<size_t, size>& arr);

//
// adapted from https://www.codeproject.com/Articles/21335/Combinations-in-C-Part
//! \brief Generates the next combination of length degree (with repetition) chosen from 0,...,
//! n_features-1 in lexicograpgical order from given array (in-place!).
template <size_t n_features, size_t degree>
constexpr inline bool next_combination_w_r(std::array<size_t, degree>& inp);

//
//! \brief Generates every possible combination of length degree (with repetition) of integers from
//! 0,..., n_features-1 and writes them into one array.
template <size_t n_features, size_t degree, size_t n_poly_features>
constexpr inline std::array<size_t, n_poly_features * degree> generate_combs();

//
//! \brief Interface for the user to generate a static loop. Function for call to static loop
//! implementation.
template <size_t start, size_t end, typename fctor, typename... fctor_types>
inline constexpr void static_loop(fctor_types&&... fctor_args);

//
//! \brief Function for reading eigen matrix from textfile. Delimiter can be specified and
//! dimensions have to be known beforehand.
template <typename m_type, size_t rows, size_t cols>
Eigen::Matrix<m_type, rows, cols> load_matrix(const std::string& file, const char delim = ',');

//
// adapted from
// https://www.codeproject.com/Articles/857354/Compile-Time-Loops-with-Cplusplus-Creating-a-Gener
//! \brief Implements a loop which is evaluated at compile time. Used here for instantiation of
//! multiple templates depending on the loop index.
template <size_t for_start, size_t for_end, typename fctor, typename... fctor_types>
struct static_loop_impl
{
    static inline constexpr void loop(fctor_types&&... fctor_args)
    {
        iteration(std::integral_constant<size_t, for_start>(),
                  std::forward<fctor_types>(fctor_args)...);
    }

private:
    //! \brief End of loop: end index reached, doing nothing causes termination.
    static inline constexpr void iteration(std::integral_constant<size_t, for_end + 1>,
                                           fctor_types&&...)
    {
    }

    //! \brief One loop iteration: execute func and next iteration with advanced index.
    template <size_t index>
    static inline constexpr void iteration(std::integral_constant<size_t, index>,
                                           fctor_types&&... fctor_args)
    {
        fctor::template func<index>(std::forward<fctor_types>(fctor_args)...);
        iteration(std::integral_constant<size_t, index + 1>(),
                  std::forward<fctor_types>(fctor_args)...);
    }
};

//! \brief Class for generating all polynomial interactions from `n_features` variables of degree
//! `degree`
template <size_t n_features, size_t degree>
struct polynomial_features
{
private:
    static constexpr auto n_poly_features = binom<n_features>(degree);
    static constexpr auto comb_size       = []() constexpr
    {
        if constexpr (degree == 0)
            return 1;
        else
            return degree * n_poly_features;
    }
    ();

    const std::array<size_t, comb_size> combinations = []() constexpr
    {
        if constexpr (degree == 0)
            return generate_combs<n_features, 1, n_poly_features>();
        else
            return generate_combs<n_features, degree, n_poly_features>();
    }
    ();

public:
    const std::array<size_t, comb_size>& operator()() const
    {
        return combinations;
    }

    constexpr polynomial_features()
    {
    }

    //! \brief Write combinations and some metadata into a stream in a human readable way
    friend std::ostream& operator<<(std::ostream& stream, const polynomial_features& poly)
    {
        stream << "Degree: " << degree << " number of features: " << n_features
               << " number of poly features: " << n_poly_features << '\n';

        auto combs_iterator = poly().begin();
        for (size_t j = 0; j < n_poly_features; ++j)
        {
            stream << j + 1 << ".  ";
            for (size_t i = 0; i < degree; ++i)
            {
                stream << *combs_iterator << ' ';
                combs_iterator++;
            }
            stream << '\n';
        }
        return stream;
    }
};

// we need an additional wrapper here because we can't do partial function specialization
// otherwise we could maybe implement this as free function
// \brief Class for collecting all the polynomial interactions from `n_features` variables of degree
// 0 to `degree`
template <size_t n_features, size_t degree, typename realtype, bool eigen_feat_vector = true>
struct poly_feat_collect
{
private:
    static constexpr auto number_polys = f_array<degree + 1>(binom<n_features>);
    static constexpr auto idx_loc      = f_array<degree + 1>(summed_binom<n_features>);

    using comb_vector = typename std::array<size_t, 1 + sum_weighted(number_polys)>;
    const comb_vector all_combinations = accumulate_combinations();

    constexpr inline comb_vector accumulate_combinations()
    {
        comb_vector ret{};
        static_loop<0, degree, poly_feat_collect<n_features, degree, realtype>>(ret);
        return ret;
    }

public:
    static constexpr auto feature_vec_length = sum(number_polys);

    template <size_t sz>
    using vector_feat = typename std::conditional<eigen_feat_vector, Eigen::Matrix<realtype, sz, 1>,
                                                  std::array<realtype, sz>>::type;

    // helper function for static loop
    template <size_t index>
    static inline constexpr void func(comb_vector& comb_collection)
    {
        polynomial_features<n_features, index> poly;
        if constexpr (index == 0)
        {
            comb_collection[index] = 0;
        }
        else
        {
            auto comb_collection_it = comb_collection.begin() + idx_loc[index - 1];
            std::copy(poly().cbegin(), poly().cend(), comb_collection_it);
        }
    }

    //! Most important method: Taking vector of features and returning the vector with polynomial
    //! features up to degree `degree`
    vector_feat<feature_vec_length> operator()(vector_feat<n_features>& features) const
    {
        vector_feat<feature_vec_length> res;
        res.fill(1.0);
        size_t k = 0;

        auto comb_iterator = all_combinations.begin();
        for (size_t i = 0; i <= degree; i++)
        {
            for (size_t j = 0; j < number_polys[i]; j++)
            {
                for (size_t l = 0; l < ((i < 1) ? 1 : i); l++)
                {
                    if (k == 0)
                    {
                        res[k] = 1;
                    }
                    else
                    {
                        res[k] *= features[*comb_iterator];
                    }
                    comb_iterator++;
                }
                k++;
            }
        }
        return res;
    }

    constexpr poly_feat_collect()
    {
    }

    //! \brief Write combinations and some metadata into a stream in a human readable way
    friend std::ostream& operator<<(std::ostream& stream, const poly_feat_collect& poly)
    {
        stream << "Degree: " << degree << " number of features: " << n_features
               << "\n[degree] number of combinations: ";
        for (size_t i = 0; i < poly.number_polys.size(); i++)
        {
            stream << "    [" << i << "] " << poly.number_polys[i];
        }
        stream << "\nFeature vector length: " << poly.feature_vec_length << '\n';

        size_t k = 0, m = 0;

        for (size_t i = 0; i <= degree; i++)
        {
            for (size_t j = 0; j < number_polys[i]; j++)
            {
                stream << k + 1 << ".  ";
                for (size_t l = 0; l < ((i < 1) ? 1 : i); l++)
                {
                    stream << poly.all_combinations[m] << ' ';
                    m++;
                }
                stream << '\n';
                k++;
            }
        }
        return stream;
    }
};

//////////////////////////////////////////////////////////////////////////
// Implementation of free functions                                     //
//////////////////////////////////////////////////////////////////////////

template <size_t start, size_t end, typename fctor, typename... fctor_types>
inline constexpr void static_loop(fctor_types&&... fctor_args)
{
    static_loop_impl<start, end, fctor, fctor_types...>::loop(
      std::forward<fctor_types>(fctor_args)...);
}

template <size_t n_features>
constexpr inline size_t binom(size_t k)
{
    auto n = n_features + k - 1;
    return (k > n) ? 0 :                                                         // special cases
             (k == 0 || k == n) ? 1 :                                            //
               (k == 1 || k == n - 1) ? n : (binom<n_features>(k - 1) * n) / k;  // recursive call
}

//
template <size_t n_features>
constexpr inline size_t summed_binom(size_t k)
{
    auto res = k * binom<n_features>(k);
    return (k == 0) ? 1 : res + summed_binom<n_features>(k - 1);
}

//
template <size_t... i, typename fct>
constexpr inline auto f_array(const fct& f, std::index_sequence<i...>)
{
    return std::array<size_t, sizeof...(i)>{{f(i)...}};
}

//
template <size_t size, typename fct>
constexpr inline auto f_array(const fct& f)
{
    return f_array(f, std::make_index_sequence<size>{});
}

//
template <size_t size, typename T>
constexpr inline size_t sum(const std::array<T, size>& arr)
{
    size_t res = 0;
    for (size_t i = 0; i < size; ++i)
        res += arr[i];

    return res;
}

//
template <size_t size>
constexpr inline size_t sum_weighted(const std::array<size_t, size>& arr)
{
    size_t res = 0;
    for (size_t i = 0; i < size; ++i)
        res += i * arr[i];

    return res;
}

//
template <size_t n_features, size_t degree>
constexpr inline bool next_combination_w_r(std::array<size_t, degree>& inp)
{
    if constexpr (n_features == 0 || degree == 0)
    {
        return false;
    }

    bool reached_set_end = false;
    for (auto i = degree - 1; i >= 0; --i)
    {
        if (i == 0 && inp[i] == n_features - 1)
        {
            return false;
        }

        if (reached_set_end)
        {
            if (inp[i] != n_features - 1)
            {
                auto level = inp[i] + 1;
                for (auto j = i; j < degree; ++j)
                    inp[j] = level;

                return true;
            }
        }

        if (inp[i] == n_features - 1)
        {
            reached_set_end = true;
        }
        else if (inp[i] < n_features - 1)
        {
            (inp[i])++;
            return true;
        }
    }
    return true;
}

//
template <size_t n_features, size_t degree, size_t n_poly_features>
constexpr inline std::array<size_t, n_poly_features * degree> generate_combs()
{
    std::array<size_t, n_poly_features * degree> combs{};
    auto iter = combs.begin();

    using comb_array = typename std::array<size_t, degree>;
    comb_array single_comb{};
    bool has_next = true;
    for (size_t j = 0; j < n_poly_features; ++j)
    {
        std::copy(single_comb.begin(), single_comb.end(), iter);
        iter += degree;
        assert(has_next && "No next combination");
        has_next = next_combination_w_r<n_features>(single_comb);
    }
    return combs;
}

template <typename m_type, size_t rows, size_t cols>
Eigen::Matrix<m_type, rows, cols> load_matrix(const std::string& file, const char delim)
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
}  // namespace Gascoigne

#endif
