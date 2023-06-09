#pragma once
#include <array>
#include <complex>
#include <cmath>

#define PI 3.14159265358979323846

//8, 3
template <size_t N, size_t N2>
constexpr std::array<std::array<double, N>, N2> calc_cosines()
{
    std::array<std::array<double, N>, N2> cosines {0};
    for (size_t i = 0; i < cosines.size(); i++) {
        unsigned int pow2 = 1 << i + 1;
        for (size_t j = 0; j < cosines.at(i).size(); j++) {
            cosines.at(i).at(j) = std::cos(2 * PI * j / pow2);
        }
    }
    return cosines;
}

//8, 3
template <size_t N, size_t N2>
constexpr std::array<std::array<double, N>, N2> calc_sines()
{
    std::array<std::array<double, N>, N2> sines {0};
    for (size_t i = 0; i < sines.size(); i++) {
        unsigned int pow2 = 1 << i + 1;
        for (size_t j = 0; j < sines.at(i).size(); j++) {
            sines.at(i).at(j) = -std::sin(2 * PI * j / pow2);
        }
    }
    return sines;
}

//determine highest bit and return its position
//std method is not constexpr so thats why this was created
constexpr unsigned int fft_log2(unsigned int num)
{
    unsigned int retValue {0};
    while(num = num >> 1) {
        retValue++;
    }
    return retValue;
}

//determine if num is a power of 2
constexpr bool isPowerOf2(unsigned int num)
{
    if (num == 0)
        return false;
    else
        return (num & (num - 1)) == 0;
}

template <size_t N> uint reverse(uint n)
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
void fft(std::array<std::complex<double>, N>& data) {

    static_assert(isPowerOf2(N), "FFT size must be a power of 2!");
    using namespace std::complex_literals;

    //compile time calculations
    auto cosines = calc_cosines<N, fft_log2(N)>();
    auto sines = calc_sines<N, fft_log2(N)>();

    for (uint i = 1; i < N / 2; i++) {
        uint i2 = reverse<N>(i);
        if (i2 != i)
            std::swap(data.at(i), data.at(i2));
    }

    int i2 = 0;
    for (int i = 1; i < N; i = i << 1) {

        for (int j = 0; j < N; j+=(i<<1)) {

            for (int k = 0; k < i; k++) {
                int index1 = j+k;
                int index2 = j+k+i;
                std::complex<double> W = cosines.at(i2).at(index1) + sines.at(i2).at(index1) * 1.0i;
                std::complex<double> W2 = cosines.at(i2).at(index2) + sines.at(i2).at(index2) * 1.0i;
                std::complex<double> val1 = data.at(index1) + data.at(index2) * W;
                std::complex<double> val2 = data.at(index1) + data.at(index2) * W2;
                data.at(index1) = val1;
                data.at(index2) = val2;
            }
        }
        i2++;
    }

}