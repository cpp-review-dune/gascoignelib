#ifndef POLYNOMIAL_FEATURES
#define POLYNOMIAL_FEATURES

#include <utility>
#include <array>
#include <vector>
#include <eigen3/Eigen/Dense>

namespace Gascoigne
{
template <size_t n_features, size_t degree>
struct polynomial_features
{
    static constexpr inline size_t binom(size_t k)
    {
        auto n = n_features + k - 1;
        return (k > n) ? 0 :                                             // special cases
                 (k == 0 || k == n) ? 1 :                                //
                   (k == 1 || k == n - 1) ? n : (binom(k - 1) * n) / k;  // recursive call
    }
    //! \brief Returns compile-time array initialized with index dependent values. Function is given
    //! as reference (e.g. binom).
    template <size_t... i, typename fct>
    static constexpr inline auto f_array(const fct& f, std::index_sequence<i...>)
    {
        return std::array<size_t, sizeof...(i)>{{f(i)...}};
    }

    template <size_t size, typename fct>
    static constexpr inline auto f_array(const fct& f)
    {
        return f_array(f, std::make_index_sequence<size>{});
    }

    template <size_t size>
    static constexpr inline size_t sum_weighted(const std::array<size_t, size>& arr)
    {
        size_t res = 1;
        for (auto i = 1; i < size; ++i)
            res += i * arr[i];

        return res;
    }

    template <size_t size>
    static constexpr inline size_t sum(const std::array<size_t, size>& arr)
    {
        size_t res = 0;
        for (auto i = 0; i < size; ++i)
            res += arr[i];

        return res;
    }

    static constexpr auto n_poly_features = f_array<degree + 1>(binom);
    static constexpr size_t featvec_size  = sum(n_poly_features);
    static constexpr size_t comb_size     = sum_weighted(n_poly_features);

    // adapted from https://www.codeproject.com/Articles/21335/Combinations-in-C-Part
    template <size_t sz>
    static inline bool next_combination(std::array<size_t, sz>& inp)
    {
        auto k = inp.size();
        if (n_features == 0 || k == 0)
        {
            return false;
        }

        bool reached_set_end = false;
        for (auto i = k - 1; i >= 0; --i)
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
                    for (auto j = i; j < k; ++j)
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

    // static inline std::array<size_t, comb_size> generate_combs()
    // {
    //     std::array<size_t, comb_size> combs;
    //     combs.fill(0);
    //     auto iter = combs.begin();
    //     for (auto i = 0; i <= degree; i++)
    //     {
    //         using comb_array       = size_t[i];
    //         comb_array single_comb;
    //         for(auto j = 0;j< i;j++)
    //             single_comb[j] = 0;
    //         for (auto j = 0; j < n_poly_features[i]; ++j)
    //         {
    //             std::copy(single_comb.begin(), single_comb.end(), iter);
    //             iter += (single_comb.end() - single_comb.begin());
    //             next_combination(single_comb);
    //         }
    //     }
    //     return combs;
    // }

    static inline void
      // TODO implememt as free functions, member functions need full specialization
      template <size_t curr_deg>
      static inline void generate_comb_impl(std::array<size_t, comb_size>& combs,
                                            typename std::array<size_t, comb_size>::iterator iter)
    {
        generate_comb_impl<curr_deg - 1>(combs, iter - n_poly_features[curr_deg - 1]);
        using comb_array = typename std::array<size_t, curr_deg>;
        comb_array single_comb;
        single_comb.fill(0);
        for (auto j = 0; j < n_poly_features[curr_deg]; ++j)
        {
            std::copy(single_comb.begin(), single_comb.end(), iter);
            iter += (single_comb.end() - single_comb.begin());
            next_combination(single_comb);
        }
    }

    template <>
    static inline void generate_comb_impl<0>(std::array<size_t, comb_size>& combs,
                                             typename std::array<size_t, comb_size>::iterator iter)
    {
        auto is_valid = (iter == combs.begin());
        // if valid we want entry to be zero
        combs[0] = static_cast<size_t>(!is_valid);
    }

    const std::array<size_t, comb_size> combinations = generate_combs();

    using vector_feat = Eigen::Matrix<double, n_features, 1>;
    constexpr polynomial_features()
    {
    }

    friend std::ostream& operator<<(std::ostream& stream, const polynomial_features& poly)
    {
        stream << "Number of polynomial features:";
        for (auto i = 0; i < n_poly_features.size(); ++i)
        {
            stream << "Degree: " << i << " number: " << n_poly_features[i] << '\n';
        }
        stream << "Sum: " << featvec_size << '\n'
               << "Size of combinations vector: " << comb_size << '\n';

        auto combs_iterator = poly.combinations.begin();
        for (auto i = 0; i <= degree; ++i)
        {
            for (auto j = 0; j < n_poly_features[i]; ++j)
            {
                for (auto i = 0; i <= degree; ++i)
                {
                    stream << *combs_iterator << ' ';
                    combs_iterator++;
                }
                stream << '\n';
            }
        }
        return stream;
    }

    // Eigen::Matrix<double, n_poly_features, 1>& operator()(const vector_feat& feature_vector)
    // const
    // {
    // }
};

}  // namespace Gascoigne

#endif
