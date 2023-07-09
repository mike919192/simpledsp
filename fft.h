#pragma once
#include <array>
#include <complex>
#include <cmath>
#include "math.h"
#include <algorithm>

//determine highest bit and return its position
//std method is not constexpr so thats why this was created
constexpr uint fft_log2(uint num)
{
    uint retValue {0};
    while( (num = num >> 1) > 0) {
        retValue++;
    }
    return retValue;
}

template <size_t N>
using trig_array = std::array<std::array<double, N>, fft_log2(N)>;

template <size_t N>
using coeff_array = std::array<std::array<std::complex<double>, N>, fft_log2(N)>;

template <class T, size_t N>
using complex_array = std::array<std::complex<T>, N>;

template <size_t N>
constexpr trig_array<N> calc_cosines_naive()
{
    trig_array<N> cosines {0};
    for (size_t i = 0; i < cosines.size(); i++) {
        uint pow2 = 1 << i + 1;
        for (size_t j = 0; j < cosines.at(i).size(); j++) {
            cosines.at(i).at(j) = std::cos(2 * M_PI * j / pow2);
        }
    }
    return cosines;
}

template <size_t N>
constexpr trig_array<N> calc_sines_naive()
{
    trig_array<N> sines {0};
    for (size_t i = 0; i < sines.size(); i++) {
        uint pow2 = 1 << i + 1;
        for (size_t j = 0; j < sines.at(i).size(); j++) {
            sines.at(i).at(j) = -std::sin(2 * M_PI * j / pow2);
        }
    }
    return sines;
}

class sine_calculator
{
    public:
    constexpr static double Value0() { return 0.0; }
    constexpr static double Value90() { return 1.0; }
    constexpr static double Value(double rad) { return std::sin(rad); };

    //symmetry sign changes
    //return either +1 or -1 to indicate sign change
    constexpr static double Sym0() { return -1.0; }
    constexpr static double Sym90() { return 1.0; }
};

class cosine_calculator
{
    public:
    constexpr static double Value0() { return 1.0; }
    constexpr static double Value90() { return 0.0; }
    constexpr static double Value(double rad) { return std::cos(rad); };

    //symmetry sign changes
    //return either +1 or -1 to indicate sign change
    constexpr static double Sym0() { return 1.0; }
    constexpr static double Sym90() { return -1.0; }
};

class reverse_fft
{
    public:
    constexpr static double Sign() { return 1.0; }

    template<size_t N>
    constexpr static void ScaleValues(complex_array<double, N> & data)
    {
        std::for_each(data.begin(), data.end(), [](auto &n) { n *= (1.0 / N); });
    }
};

class forward_fft
{
    public:
    constexpr static double Sign() { return -1.0; }

    template<size_t N>
    constexpr static void ScaleValues(complex_array<double, N> &) { }
};

template <size_t N, class T>
constexpr trig_array<N> calc_trigs()
{
    trig_array<N> values {0};
    for (size_t i = 0; i < values.size(); i++) {
        uint pow2 = 1 << (i + 1);
        
        //first element at 0 deg
        values.at(i).at(0) = T::Value0();

        //i = 0 is special
        //just alternate +1 and -1
        if (i == 0) {
            for (size_t j = 1; j < values.at(i).size(); j++) {
                values.at(i).at(j) = values.at(i).at(j-1) * -1.0;
            }
        } else {
            //calculate up to but not including 90 deg
            uint num = 1 << (i - 1);
            for (size_t j = 1; j < num; j++) {
                values.at(i).at(j) = T::Value(2 * M_PI * j / pow2);
            }
            //next element is 90 deg
            values.at(i).at(num) = T::Value90();    

            //now use symmetry to get the rest of the values
            int dir = -1;
            double sign = T::Sym90();
            uint bouncyIndex = num;

            for (size_t j = num + 1; j < values.at(i).size(); j++) {
                bouncyIndex += dir;
                values.at(i).at(j) = values.at(i).at(bouncyIndex) * sign;                
                if (bouncyIndex == 0) {
                    dir = 1;
                    sign *= T::Sym0();
                } else if (bouncyIndex == num) {
                    dir = -1;
                    sign *= T::Sym90();
                }
            }
        }
    }
    return values;
}

//form the W coefficients from the cosines and sines arrays
template<size_t N, class T>
constexpr coeff_array<N> calc_wCoeffs()
{
    //auto cosines = calc_cosines_naive<N>();
    auto cosines = calc_trigs<N, cosine_calculator>();
    
    //auto sines = calc_sines_naive<N>();
    auto sines = calc_trigs<N, sine_calculator>();

    coeff_array<N> wCoeffs {0};

    for (size_t i = 0; i < cosines.size(); i++) {

        for (size_t j = 0; j < cosines.at(i).size(); j++) {
            wCoeffs.at(i).at(j) = std::complex(cosines.at(i).at(j), T::Sign() * sines.at(i).at(j));
        }
    }
    return wCoeffs;
}

//determine if num is a power of 2
constexpr bool isPowerOf2(uint num)
{
    if (num == 0)
        return false;
    else
        return (num & (num - 1)) == 0;
}

//determine if num is a power of 4
constexpr bool isPowerOf4(uint num)
{
    return isPowerOf2(num) && (fft_log2(num) % 2 == 0);
}

//reverse the order of bits in the mask of N - 1
template <size_t N> 
constexpr uint reverse(uint n)
{
    uint retValue {0};
    uint shift = fft_log2(N) - 1;
    uint upperBit = N >> 1;
    uint lowerBit = 1;
    while (upperBit > lowerBit) {
        retValue |= (n & upperBit) >> shift;
        retValue |= (n & lowerBit) << shift;
        upperBit = upperBit >> 1;
        lowerBit = lowerBit << 1;
        shift -= 2;
    }
    if (upperBit == lowerBit) {
        retValue |= n & upperBit;
    }
    return retValue;
}

template <size_t N> 
constexpr std::array<uint, N> calc_swap_lookup()
{
    //first create an array where every element value is reversed bits of the index
    std::array<uint, N> swapLookup {0};
    for (uint i = 0; i < N; i++) {
        swapLookup.at(i) = reverse<N>(i);
    }
    //then go through one more time and for every pair, unreverse the one with the higher index
    //this prevents the swap from occuring twice, which would undo the swap
    for (uint i = 1; i < N - 1; i++) {
        uint i2 = swapLookup.at(i);
        if (i2 != i) {
            swapLookup.at(i2) = i2;
        }
    }
    return swapLookup;
}

template <size_t N>
constexpr std::array<uint, N> calc_radix4_swap_lookup()
{
    //first create an array with the output pattern of the radix 4
    std::array<uint, N> swapLookup {0};
    uint index {0};
    for (uint i = 0; i < 4; i++) {

        for (uint j = 0; j < N / 4; j++) {
            swapLookup.at(index) = j * 4 + i;
            index++;
        }        
    }
    //then go through one more time and for every pair, unreverse the one with the higher index
    //this prevents the swap from occuring twice, which would undo the swap
    for (uint i = 1; i < N - 1; i++) {
        uint i2 = swapLookup.at(i);
        if (i2 != i) {
            swapLookup.at(i2) = i2;
        }
    }
    return swapLookup;
}

template <class T = forward_fft, size_t N>
void fft(complex_array<double, N> & data) {

    static_assert(isPowerOf2(N), "FFT size must be a power of 2!");

    //compile time calculations
    constexpr auto wCoeffs = calc_wCoeffs<N, T>();
    constexpr auto swapLookup = calc_swap_lookup<N>();

    //decimation in time
    //perform swap on inputs
    for (uint i = 1; i < N - 1; i++) {
        uint i2 = swapLookup.at(i);
        if (i2 != i)
            std::swap(data.at(i), data.at(i2));
    }

    //outer most loop is the FFT stages
    //2pt FFT -> 4pt FFT -> 8pt FFT -> etc
    uint i2 = 0;
    for (uint i = 1; i < N; i = i << 1) {

        //the inside 2 loops perform the butterfly pattern
        for (uint j = 0; j < N; j += (i << 1)) {

            for (uint k = 0; k < i; k++) {
                uint index1 = j + k;
                uint index2 = j + k + i;

                //original formula page 137
                //std::complex<double> val1 = data.at(index1) + data.at(index2) * wCoeffs.at(i2).at(index1);
                //std::complex<double> val2 = data.at(index1) + data.at(index2) * wCoeffs.at(i2).at(index2);

                //optimization from page 145
                std::complex<double> temp = data.at(index2) * wCoeffs.at(i2).at(index1);
                std::complex<double> val1 = data.at(index1) + temp;
                std::complex<double> val2 = data.at(index1) - temp;

                data.at(index1) = val1;
                data.at(index2) = val2;
            }
        }
        i2++;
    }

    //if reverse, scale the data by 1/N
    //if forward, does nothing
    T::ScaleValues(data);

}

template <class T = forward_fft, size_t N>
void fft_radix4(complex_array<double, N> & data) {

    static_assert(isPowerOf4(N), "FFT radix 4 size must be a power of 4!");

    //compile time calculations
    constexpr auto wCoeffs = calc_wCoeffs<N, T>();
    constexpr auto swapLookup = calc_radix4_swap_lookup<N>();

    uint i2 = 0;
    for (uint i = N / 4; i >= 1; i = i >> 2) {

        uint j2 = 0;
        for (uint j = 0; j < N; j += (i << 2)) {

            for (uint k = 0; k < i; k++) {
                uint index1 = k + j;
                uint index2 = k + i + j;
                uint index3 = k + 2 * i + j;
                uint index4 = k + 3 * i + j;
                uint coeff_subscript = 2 * i2 + 1;
                uint coeff_index1 = (index1 - j) * i2 * i2 * (j2 % 4);
                uint coeff_index2 = (index2 - j) * i2 * i2 * (j2 % 4);
                uint coeff_index3 = (index3 - j) * i2 * i2 * (j2 % 4);
                uint coeff_index4 = (index4 - j) * i2 * i2 * (j2 % 4);

                std::complex<double> temp1 = data.at(index1) * wCoeffs.at(coeff_subscript).at(coeff_index1);
                std::complex<double> temp2 = data.at(index2) * wCoeffs.at(coeff_subscript).at(coeff_index2);
                std::complex<double> temp3 = data.at(index3) * wCoeffs.at(coeff_subscript).at(coeff_index3);
                std::complex<double> temp4 = data.at(index4) * wCoeffs.at(coeff_subscript).at(coeff_index4);
                std::complex<double> temp2_timesi = std::complex<double>(-temp2.imag(), temp2.real());
                std::complex<double> temp4_timesi = std::complex<double>(-temp4.imag(), temp4.real());

                data.at(index1) = temp1 + temp2 + temp3 + temp4;
                data.at(index2) = temp1 - temp2_timesi - temp3 + temp4_timesi;
                data.at(index3) = temp1 - temp2 + temp3 - temp4;
                data.at(index4) = temp1 + temp2_timesi - temp3 - temp4_timesi;
            }
            j2++;
        }
        i2++;
    }

    for (uint i = 1; i < N - 1; i++) {
        uint i2 = swapLookup.at(i);
        if (i2 != i)
            std::swap(data.at(i), data.at(i2));
    }

    //if reverse, scale the data by 1/N
    //if forward, does nothing
    T::ScaleValues(data);

}
