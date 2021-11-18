#pragma once
#include "valarray"
#include "functional"
#include "algorithm"

namespace num
{
    template<typename T>
    auto newton_polynomial(const std::valarray<T>& x, std::valarray<T> f);

    template<typename T>
    std::valarray<T> diff(const std::valarray<T>& f, const std::valarray<T>& Pn);

    template<typename T>
    T error(const std::valarray<T>& f, const std::valarray<T>& Pn);

    template<typename T, typename F>
    std::valarray<T> apply(const std::valarray<T>& x, const F& f);
}

namespace num
{
    template<typename T>
    auto newton_polynomial(const std::valarray<T>& x, std::valarray<T> f)
    {
        auto n = x.size() - 1;
        std::valarray<T> d(n + 1);
        d[0] = f[0];

        for (std::size_t i = 1; i <= n; ++i)
        {
            for (std::size_t k = 0; k <= n - i; ++k)
                f[k] = (f[k + 1] - f[k]) / (x[i + k] - x[k]);
            d[i] = f[0];
        }

        return [d, x](T X) -> T
        {
            auto n = d.size();
            T prod = 1, sum = 0;

            for (std::size_t i = 0; i < n; ++i)
            {
                sum += d[i] * prod;
                prod *= (X - x[i]);
            }

            return sum;
        };
    }

    template<typename T>
    std::valarray<T> diff(const std::valarray<T>& f, const std::valarray<T>& Pn)
    {
        return std::abs(f - Pn);
    }

    template<typename T>
    T error(const std::valarray<T>& f, const std::valarray<T>& Pn)
    {
        return diff(f, Pn).max();
    }

    template<typename T, typename F>
    std::valarray<T> apply(const std::valarray<T>& x, const F& f)
    {
        auto n = x.size();
        std::valarray<T> y(n);
        for (std::size_t i = 0; i < n; ++i)
            y[i] = f(x[i]);
        return y;
    }
}
