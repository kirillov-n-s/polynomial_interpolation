#pragma once
#include "valarray"
#include "matplot/matplot.h"

namespace num
{
    template<typename T>
    std::valarray<T> uniform(T a, T b, std::size_t n);

    template<typename T>
    std::valarray<T> chebyshev(T a, T b, std::size_t n);

    template<typename T>
    std::valarray<T> subsplit(const std::valarray<T>& x, std::size_t k);
}

namespace num
{
    template<typename T>
    std::valarray<T> uniform(T a, T b, std::size_t n)
    {
        std::valarray<T> i(n + 1);
        std::iota(begin(i), end(i), 0);
        return a + (b - a) * i / n;
    }

    template<typename T>
    std::valarray<T> chebyshev(T a, T b, std::size_t n)
    {
        std::valarray<T> i(n + 1);
        std::iota(begin(i), end(i), 0);
        return (a + b) / 2 + (b - a) / 2 * std::cos((2 * i + 1) / (2 * n + 2) * matplot::pi);
    }

    template<typename T>
    std::valarray<T> subsplit(const std::valarray<T>& x, std::size_t k)
    {
        auto n = x.size();
        auto m = (n - 1) * k + 1;
        std::valarray<T> res(m);

        for (std::size_t i = 0; i < n - 1; ++i)
        {
            auto tmp = uniform(x[i], x[i + 1], k);
            for (std::size_t j = 0; j < k; ++j)
                res[i * k + j] = tmp[j];
        }

        res[m - 1] = x[n - 1];
        return res;
    }
}
