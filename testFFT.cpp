
#include "catch2/catch_amalgamated.hpp"
#include "fft.h"

TEST_CASE("FFT Radix 2 test", "[single-file]")
{
    constexpr uint N {1024};
    constexpr uint n {200};
    complex_array<double, N> s {0};
    
    for (size_t i = 0; i < s.size(); i++) {
        s.at(i) = std::cos(n * 2 * M_PI * i / N);
    }

    complex_array<double, N> S {0};
    S.at(n) = N / 2;
    S.at(N - n) = N / 2;

    SECTION("FFT Forward Test")
    {
        fft(s);

        std::array<double, N> error {0};

        for (size_t i = 0; i < error.size(); i++) {
            error.at(i) = std::abs(S.at(i) - s.at(i));
        }

        double max_error = *std::max_element(error.begin(), error.end());

        REQUIRE(max_error < 4 * N * std::numeric_limits<double>::epsilon());
    }

    SECTION("FFT Reverse Test")
    {        
        fft<reverse_fft>(S);

        std::array<double, N> error {0};

        for (size_t i = 0; i < error.size(); i++) {
            error.at(i) = std::abs(S.at(i) - s.at(i));
        }

        double max_error = *std::max_element(error.begin(), error.end());

        REQUIRE(max_error < 4 * N * std::numeric_limits<double>::epsilon());
    }
}

TEST_CASE("FFT benchmark", "[single-file]")
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
