#include <string>
#include <sstream>
#include <iostream>
#include <utility>
#include <array>
#include "polynomial_features.h"

int main()
{
    // Gascoigne::polynomial_features<10, 3> polyy;
    // std::cout << polyy << '\n';
    // auto p = polyy::

    Gascoigne::poly_feat_collect<20, 5, double, false> polyx;
    // std::cout << polyx << '\n';
    std::array<double, 20> x = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                                11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    auto feat                = polyx(x);

    for (auto x : feat)
    {
        std::cout << x << ' ';
    }
    // auto x = polyx.all_combinations;
    // for (auto y : x)
    //     std::cout << y << ' ';
    // std::cout << x.size() << '\n';
}
