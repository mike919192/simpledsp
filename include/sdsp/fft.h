#pragma once
#include <array>
#include <complex>
#include <cmath>
#include "math.h"
#include <algorithm>

namespace sdsp
{
//determine highest bit and return its position
//std method is not constexpr so thats why this was created
constexpr uint log2(uint num)
{
    uint ret_value{ 0 };
    while ((num = num >> 1) > 0u) {
        ret_value++;
    }
    return ret_value;
}

constexpr uint log4(uint num)
{
    uint ret_value{ 0 };
    while ((num = num >> 2) > 0u) {
        ret_value++;
    }
    return ret_value;
}

//determine if num is a power of 2
constexpr bool is_power_of_2(uint num)
{
    if (num == 0)
        return false;
    else
        return (num & (num - 1)) == 0;
}

//determine if num is a power of 4
constexpr bool is_power_of_4(uint num)
{
    return is_power_of_2(num) && (log2(num) % 2 == 0);
}

template <size_t n_t>
using trig_array_t = std::array<std::array<double, n_t>, log2(n_t)>;

template <size_t n_t>
using coeff_array_t = std::array<std::array<std::complex<double>, n_t>, log2(n_t)>;

template <size_t n_t>
using complex_array_t = std::array<std::complex<double>, n_t>;

template <size_t n_t, class trig_t>
constexpr trig_array_t<n_t> calc_trigs_naive()
{
    trig_array_t<n_t> trigs{ 0 };
    for (size_t i{ 0 }; i < trigs.size(); i++) {
        uint pow2{ 1u << (i + 1u) };
        for (size_t j{ 0 }; j < trigs.at(i).size(); j++) {
            trigs.at(i).at(j) = trig_t::value(2 * M_PI * j / pow2);
        }
    }
    return trigs;
}

class sine_calculator {
public:
    constexpr static double value_0()
    {
        return 0.0;
    }
    constexpr static double value_90()
    {
        return 1.0;
    }
    constexpr static double value(double rad)
    {
        return std::sin(rad);
    };

    //symmetry sign changes
    //return either +1 or -1 to indicate sign change
    constexpr static double sym_0()
    {
        return -1.0;
    }
    constexpr static double sym_90()
    {
        return 1.0;
    }
};

class cosine_calculator {
public:
    constexpr static double value_0()
    {
        return 1.0;
    }
    constexpr static double value_90()
    {
        return 0.0;
    }
    constexpr static double value(double rad)
    {
        return std::cos(rad);
    };

    //symmetry sign changes
    //return either +1 or -1 to indicate sign change
    constexpr static double sym_0()
    {
        return 1.0;
    }
    constexpr static double sym_90()
    {
        return -1.0;
    }
};

class reverse_fft {
public:
    constexpr static double sign()
    {
        return -1.0;
    }

    template <size_t n_t>
    constexpr static void scale_values(complex_array_t<n_t> &data)
    {
        std::for_each(data.begin(), data.end(), [](auto &n) { n *= (1.0 / n_t); });
    }
};

class forward_fft {
public:
    constexpr static double sign()
    {
        return 1.0;
    }

    template <size_t n_t>
    constexpr static void scale_values(complex_array_t<n_t> &)
    {
    }
};

template <size_t n_t, class trig_t>
constexpr trig_array_t<n_t> calc_trigs()
{
    trig_array_t<n_t> values{ 0 };
    size_t i{ 0 };
    for (auto &value : values) {
        uint pow2{ 1u << (i + 1u) };

        //first element at 0 deg
        value.at(0) = trig_t::value_0();

        //i = 0 is special
        //just alternate +1 and -1
        if (i == 0) {
            for (size_t j{ 1 }; j < value.size(); j++) {
                value.at(j) = value.at(j - 1) * -1.0;
            }
        } else {
            //calculate up to but not including 90 deg
            uint num{ 1u << (i - 1u) };
            for (uint j{ 1 }; j < num; j++) {
                value.at(j) = trig_t::value(2 * M_PI * j / pow2);
            }
            //next element is 90 deg
            value.at(num) = trig_t::value_90();

            //now use symmetry to get the rest of the values
            int dir{ -1 };
            double sign{ trig_t::sym_90() };
            uint bouncy_index{ num };

            for (size_t j{ num + 1 }; j < value.size(); j++) {
                bouncy_index += dir;
                value.at(j) = value.at(bouncy_index) * sign;
                if (bouncy_index == 0) {
                    dir = 1;
                    sign *= trig_t::sym_0();
                } else if (bouncy_index == num) {
                    dir = -1;
                    sign *= trig_t::sym_90();
                }
            }
        }
        i++;
    }
    return values;
}

//form the W coefficients from the cosines and sines arrays
template <size_t n_t, class trig_t>
constexpr coeff_array_t<n_t> calc_w_coeffs()
{
    //auto cosines {calc_trigs_naive<n_t, cosine_calculator>()};
    auto cosines{ calc_trigs<n_t, cosine_calculator>() };

    //auto sines = calc_trigs_naive<n_t, sine_calculator>();
    auto sines{ calc_trigs<n_t, sine_calculator>() };

    coeff_array_t<n_t> w_coeffs{ 0 };

    for (size_t i{ 0 }; i < cosines.size(); i++) {
        for (size_t j{ 0 }; j < cosines.at(i).size(); j++) {
            w_coeffs.at(i).at(j) = std::complex(cosines.at(i).at(j), trig_t::sign() * -1.0 * sines.at(i).at(j));
        }
    }
    return w_coeffs;
}

//reverse the order of bits/digits
template <size_t n_t, uint base_t>
constexpr uint digit_reverse(uint n)
{
    constexpr uint num_bits{ log2(base_t) };
    uint ret_value{ 0 };
    uint shift{ log2(n_t) - num_bits };
    uint upper_bits{ (base_t - 1) << shift };
    uint lower_bits{ base_t - 1 };
    while (upper_bits > lower_bits) {
        ret_value |= (n & upper_bits) >> shift;
        ret_value |= (n & lower_bits) << shift;
        upper_bits = upper_bits >> num_bits;
        lower_bits = lower_bits << num_bits;
        shift -= num_bits * 2;
    }
    if (upper_bits == lower_bits) {
        ret_value |= n & upper_bits;
    }
    return ret_value;
}

template <size_t n_t, uint base_t>
constexpr std::array<uint, n_t> calc_swap_lookup()
{
    //first create an array where every element value is reversed bits of the index
    std::array<uint, n_t> swap_lookup{ 0 };

    for (size_t i{ 0 }; i < n_t; i++) {
        swap_lookup.at(i) = digit_reverse<n_t, base_t>(static_cast<uint>(i));
    }
    //then go through one more time and for every pair, unreverse the one with the higher index
    //this prevents the swap from occuring twice, which would undo the swap
    for (size_t i{ 1 }; i < n_t - 1; i++) {
        uint i2{ swap_lookup.at(i) };
        if (i2 != i) {
            swap_lookup.at(i2) = i2;
        }
    }
    return swap_lookup;
}

template <class trig_t = forward_fft, size_t n_t>
void fft_radix2(complex_array_t<n_t> &data)
{
    static_assert(is_power_of_2(n_t), "FFT size must be a power of 2!");

    //compile time calculations
    constexpr static auto w_coeffs{ calc_w_coeffs<n_t, trig_t>() };
    constexpr static auto swap_lookup{ calc_swap_lookup<n_t, 2>() };

    //decimation in time
    //perform swap on inputs
    for (uint i{ 1 }; i < n_t - 1; i++) {
        uint i2{ swap_lookup.at(i) };
        if (i2 != i)
            std::swap(data.at(i), data.at(i2));
    }

    //outer most loop is the FFT stages
    for (uint i{ 0 }; i < log2(n_t); i++) {
        uint two_raised_to_i{ 1u << i };

        //the inside 2 loops perform the butterfly pattern
        for (uint j = 0; j < n_t; j += (two_raised_to_i << 1)) {
            for (uint k{ 0 }; k < two_raised_to_i; k++) {
                uint index1{ j + k };
                uint index2{ j + k + two_raised_to_i };

                //optimization from page 145
                std::complex<double> temp{ data.at(index2) * w_coeffs.at(i).at(index1) };
                std::complex<double> val1{ data.at(index1) + temp };
                std::complex<double> val2{ data.at(index1) - temp };

                data.at(index1) = val1;
                data.at(index2) = val2;
            }
        }
    }

    //if reverse, scale the data by 1/n_t
    //if forward, does nothing
    trig_t::scale_values(data);
}

template <class trig_t = forward_fft, size_t n_t>
void fft_radix4(complex_array_t<n_t> &data)
{
    static_assert(is_power_of_4(n_t), "FFT radix 4 size must be a power of 4!");

    //compile time calculations
    constexpr static auto w_coeffs{ calc_w_coeffs<n_t, trig_t>() };
    constexpr static auto swap_lookup{ calc_swap_lookup<n_t, 4>() };
    constexpr static uint coeff_subscript{ log2(n_t) - 1 };

    for (uint i{ 0 }; i < log4(n_t); i++) {
        uint i_group_size{ static_cast<uint>(n_t) / (4u << (2u * i)) };
        uint four_raised_to_n{ i > 0 ? 1u << ((i - 1u) * 2u) : 0 };

        uint j2{ 0 };
        for (uint j{ 0 }; j < n_t; j += (i_group_size << 2)) {
            for (uint k{ 0 }; k < i_group_size; k++) {
                uint index1{ k + j };
                uint index2{ k + i_group_size + j };
                uint index3{ k + 2 * i_group_size + j };
                uint index4{ k + 3 * i_group_size + j };
                uint coeff_index1{ (index1 - j) * four_raised_to_n * (j2 % 4) };
                uint coeff_index2{ (index2 - j) * four_raised_to_n * (j2 % 4) };
                uint coeff_index3{ (index3 - j) * four_raised_to_n * (j2 % 4) };
                uint coeff_index4{ (index4 - j) * four_raised_to_n * (j2 % 4) };

                std::complex<double> temp1{ coeff_index1 > 0 ?
                                                data.at(index1) * w_coeffs.at(coeff_subscript).at(coeff_index1) :
                                                data.at(index1) };
                std::complex<double> temp2{ coeff_index2 > 0 ?
                                                data.at(index2) * w_coeffs.at(coeff_subscript).at(coeff_index2) :
                                                data.at(index2) };
                std::complex<double> temp3{ coeff_index3 > 0 ?
                                                data.at(index3) * w_coeffs.at(coeff_subscript).at(coeff_index3) :
                                                data.at(index3) };
                std::complex<double> temp4{ coeff_index4 > 0 ?
                                                data.at(index4) * w_coeffs.at(coeff_subscript).at(coeff_index4) :
                                                data.at(index4) };
                std::complex<double> temp2_timesi{ trig_t::sign() * std::complex<double>(-temp2.imag(), temp2.real()) };
                std::complex<double> temp4_timesi{ trig_t::sign() * std::complex<double>(-temp4.imag(), temp4.real()) };

                data.at(index1) = temp1 + temp2 + temp3 + temp4;
                data.at(index2) = temp1 - temp2_timesi - temp3 + temp4_timesi;
                data.at(index3) = temp1 - temp2 + temp3 - temp4;
                data.at(index4) = temp1 + temp2_timesi - temp3 - temp4_timesi;
            }
            j2++;
        }
    }

    for (uint i{ 1 }; i < n_t - 1; i++) {
        uint i2{ swap_lookup.at(i) };
        if (i2 != i)
            std::swap(data.at(i), data.at(i2));
    }

    //if reverse, scale the data by 1/n_t
    //if forward, does nothing
    trig_t::scale_values(data);
}
}
