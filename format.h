#pragma once

#include <iostream>
#include <iomanip>
#include <limits>

//func decl
namespace num
{
    template<typename T>
    std::size_t format(std::ostream& out, bool left = true);
}

//func def
namespace num
{
    template<typename T>
    std::size_t format(std::ostream& out, bool left)
    {
        int digits = std::numeric_limits<T>::max_digits10;
        int places = digits + 7 + std::log10(std::numeric_limits<T>::max_exponent10);
        out << std::scientific << std::setprecision(digits) << (left ? std::left : std::right) << std::showpos << std::setfill(' ');
        return places;
    }
}
