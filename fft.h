#pragma once
#include <array>
#include <complex>
#include <cmath>
#include "math.h"

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
constexpr std::array<std::array<double, N>, fft_log2(N)> calc_cosines_bad()
{
    std::array<std::array<double, N>, fft_log2(N)> cosines {0};
    for (size_t i = 0; i < cosines.size(); i++) {
        uint pow2 = 1 << i + 1;
        for (size_t j = 0; j < cosines.at(i).size(); j++) {
            cosines.at(i).at(j) = std::cos(2 * M_PI * j / pow2);
        }
    }
    return cosines;
}

template <size_t N>
constexpr std::array<std::array<double, N>, fft_log2(N)> calc_cosines()
{
    std::array<std::array<double, N>, fft_log2(N)> cosines {0};
    for (size_t i = 0; i < cosines.size(); i++) {
        uint pow2 = 1 << (i + 1);
        
        //first element is always 1
        cosines.at(i).at(0) = 1.0;
        uint values90 = 1; //number of values from 0 up to and including 90deg

        //i = 0 is special
        //just alternate +1 and -1
        if (i == 0) {
            for (size_t j = 1; j < cosines.at(i).size(); j++) {
                cosines.at(i).at(j) = cosines.at(i).at(j-1) * -1.0;
            }
        }

        else {
            //calculate up to but not including 90 deg
            uint num = 1 << (i - 1);
            for (size_t j = 1; j < num; j++) {
                cosines.at(i).at(j) = std::cos(2 * M_PI * j / pow2);
            }
            //next element is 0
            cosines.at(i).at(num) = 0.0;
            values90 = num + 1;        

            //now use symmetry to get the rest of the values
            int dir = -1;
            double sign = -1.0;
            uint bouncyIndex = values90 - 1;

            for (size_t j = values90; j < cosines.at(i).size(); j++) {
                bouncyIndex += dir;
                cosines.at(i).at(j) = cosines.at(i).at(bouncyIndex) * sign;                
                if (bouncyIndex == 0) {
                    dir = 1;
                } else if (bouncyIndex == values90 - 1) {
                    dir = -1;
                    sign *= -1.0;
                }
            }
        }
    }
    return cosines;
}

template <size_t N>
constexpr std::array<std::array<double, N>, fft_log2(N)> calc_sines_bad()
{
    std::array<std::array<double, N>, fft_log2(N)> sines {0};
    for (size_t i = 0; i < sines.size(); i++) {
        uint pow2 = 1 << i + 1;
        for (size_t j = 0; j < sines.at(i).size(); j++) {
            sines.at(i).at(j) = -std::sin(2 * M_PI * j / pow2);
        }
    }
    return sines;
}

template <size_t N>
constexpr std::array<std::array<double, N>, fft_log2(N)> calc_sines()
{
    std::array<std::array<double, N>, fft_log2(N)> sines {0};
    for (size_t i = 0; i < sines.size(); i++) {
        uint pow2 = 1 << (i + 1);
        
        //first element is always 0
        sines.at(i).at(0) = 0.0;
        uint values90 = 1; //number of values from 0 up to and including 90deg

        //i = 0 is special
        //all values are 0

        if (i > 0) {
            //calculate up to but not include 90 deg
            uint num = 1 << (i - 1);
            for (size_t j = 1; j < num; j++) {
                sines.at(i).at(j) = -std::sin(2 * M_PI * j / pow2);
            }
            //next element is -1
            sines.at(i).at(num) = -1.0;
            values90 = num + 1;        

            //now use symmetry to get the rest of the values
            int dir = -1;
            double sign = 1.0;
            uint bouncyIndex = values90 - 1;

            for (size_t j = values90; j < sines.at(i).size(); j++) {
                bouncyIndex += dir;
                sines.at(i).at(j) = sines.at(i).at(bouncyIndex) * sign;                
                if (bouncyIndex == 0) {
                    dir = 1;
                    sign *= -1.0;
                } else if (bouncyIndex == values90 - 1) {
                    dir = -1;                    
                }
            }
        }
    }
    return sines;
}

//form the W coefficients from the cosines and sines arrays
template<size_t N>
constexpr std::array<std::array<std::complex<double>, N>, fft_log2(N)> calc_wCoeffs()
{
    auto cosines = calc_cosines<N>();
    auto sines = calc_sines<N>();
    std::array<std::array<std::complex<double>, N>, fft_log2(N)> wCoeffs {0};

    for (size_t i = 0; i < cosines.size(); i++) {

        for (size_t j = 0; j < cosines.at(i).size(); j++) {
            wCoeffs.at(i).at(j) = std::complex(cosines.at(i).at(j), sines.at(i).at(j));
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
void fft(std::array<std::complex<double>, N>& data) {

    static_assert(isPowerOf2(N), "FFT size must be a power of 2!");

    //compile time calculations
    constexpr auto wCoeffs = calc_wCoeffs<N>();
    constexpr auto swapLookup = calc_swap_lookup<N>();

    //decimation in time
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
                std::complex<double> val1 = data.at(index1) + data.at(index2) * wCoeffs.at(i2).at(index1);
                std::complex<double> val2 = data.at(index1) + data.at(index2) * wCoeffs.at(i2).at(index2);
                data.at(index1) = val1;
                data.at(index2) = val2;
            }
        }
        i2++;
    }

}