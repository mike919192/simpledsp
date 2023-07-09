
#include "catch2/catch_amalgamated.hpp"
#include "fft.h"

TEST_CASE("FFT test", "[single-file]")
{
    complex_array<double, 1024> complexValues {0.3535, 0.3535, 0.6464, 1.0607, 0.3535, -1.0607, -1.3535, -0.3535};

    BENCHMARK("Radix2 FFT")
    {
        fft(complexValues);
        return complexValues;
    };

    BENCHMARK("Radix4 FFT")
    {
        fft_radix4(complexValues);
        return complexValues;
    };
}
