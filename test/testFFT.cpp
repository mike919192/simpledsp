
#define CATCH_CONFIG_MAIN
#include "catch2/catch_all.hpp"
#include "sdsp/fft.h"

template <size_t N>
double calc_max_error(const sdsp::complex_array<N> &observed, const sdsp::complex_array<N> &expected)
{
    std::array<double, N> error {0};

    for (size_t i {0}; i < error.size(); i++) {
        error.at(i) = std::abs(observed.at(i) - expected.at(i));
    }

    return *std::max_element(error.begin(), error.end());
}

TEST_CASE("FFT Radix 2 test", "[single-file]")
{
    constexpr uint N {64};
    constexpr uint n {7};
    sdsp::complex_array<N> s {0};
    
    for (size_t i = 0; i < s.size(); i++) {
        s.at(i) = std::cos(n * 2 * M_PI * i / N);
    }

    sdsp::complex_array<N> S {0};
    S.at(n) = N / 2;
    S.at(N - n) = N / 2;

    SECTION("FFT Forward Test")
    {
        sdsp::fft_radix2(s);

        double max_error {calc_max_error(S, s)};

        REQUIRE(max_error < 4 * N * std::numeric_limits<double>::epsilon());
    }

    SECTION("FFT Reverse Test")
    {        
        sdsp::fft_radix2<sdsp::reverse_fft>(S);

        double max_error {calc_max_error(S, s)};

        REQUIRE(max_error < 4 * N * std::numeric_limits<double>::epsilon());
    }

    SECTION("FFT Time shift test")
    {
        //s2 shifted 90 degrees from s
        sdsp::complex_array<N> s2 {0};
    
        for (size_t i {0}; i < s2.size(); i++) {
            s2.at(i) = std::cos(n * 2 * M_PI * i / N + (M_PI / 2.0));
        }

        sdsp::fft_radix2(s);
        sdsp::fft_radix2(s2);

        sdsp::complex_array<N> S2 {0};
        S2.at(n) = std::complex<double>(0, N / 2);
        S2.at(N - n) = std::complex<double>(0, -(N / 2.0));

        double max_error {calc_max_error(S2, s2)};
        REQUIRE(max_error < 4 * N * std::numeric_limits<double>::epsilon());
    }
}

TEST_CASE("FFT Radix 2 linearity", "[single-file]")
{
    constexpr double freq1 {1000.0};
    constexpr double freq2 {500.0};
    constexpr double fs {8000.0};
    constexpr double a1 {1.5};
    constexpr double a2 {2.5};
    constexpr uint N {256};

    //test linearity
    //FFT(a1x1[n]+a2x2[n])=a1FFT(x1[n])+a2FFT(x2[n])

    sdsp::complex_array<N> x1 {0};
    sdsp::complex_array<N> x2 {0};
    for (size_t i {0}; i < x1.size(); i++) {
        double value {std::sin(2.0 * M_PI * freq1 * (1.0 / fs) * i)};
        x1.at(i) = std::complex<double>(value, 0);
        
        double value2 {std::sin(2.0 * M_PI * freq2 * (1.0 / fs) * i)};
        x2.at(i) = std::complex<double>(value2, 0);
    }

    sdsp::complex_array<N> a1x1 {x1};
    std::for_each(a1x1.begin(), a1x1.end(), [a1](auto &n) { n *= a1; });

    sdsp::complex_array<N> a2x2 {x2};
    std::for_each(a2x2.begin(), a2x2.end(), [a2](auto &n) { n *= a2; });

    sdsp::complex_array<N> a1x1_plus_a2x2 {0};
    for (size_t i {0}; i < a1x1_plus_a2x2.size(); i++) {
        a1x1_plus_a2x2.at(i) = a1x1.at(i) + a2x2.at(i);
    }

    sdsp::complex_array<N> fft_data1 {a1x1_plus_a2x2};
    sdsp::fft_radix2(fft_data1);

    sdsp::complex_array<N> fft_data2 {x1};
    sdsp::complex_array<N> fft_data3 {x2};
    sdsp::fft_radix2(fft_data2);
    sdsp::fft_radix2(fft_data3);

    sdsp::complex_array<N> a1_fft_data2 {fft_data2};
    std::for_each(a1_fft_data2.begin(), a1_fft_data2.end(), [a1](auto &n) { n *= a1; });

    sdsp::complex_array<N> a2_fft_data3 {fft_data3};
    std::for_each(a2_fft_data3.begin(), a2_fft_data3.end(), [a2](auto &n) { n *= a2; });

    sdsp::complex_array<N> a1_fft_data2_plus_a2_fft_data3 {0};
    for (size_t i {0}; i < a1_fft_data2_plus_a2_fft_data3.size(); i++) {
        a1_fft_data2_plus_a2_fft_data3.at(i) = a1_fft_data2.at(i) + a2_fft_data3.at(i);
    }

    //difference should be very close to 0
    double max_error {calc_max_error(fft_data1, a1_fft_data2_plus_a2_fft_data3)};
    REQUIRE(max_error < 4 * N * std::numeric_limits<double>::epsilon());
}

TEST_CASE("FFT Radix 4 test", "[single-file]")
{
    constexpr uint N {64};
    constexpr uint n {7};
    sdsp::complex_array<N> s {0};
    
    for (size_t i {0}; i < s.size(); i++) {
        s.at(i) = std::cos(n * 2 * M_PI * i / N);
    }

    sdsp::complex_array<N> S {0};
    S.at(n) = N / 2;
    S.at(N - n) = N / 2;

    SECTION("FFT Forward Test")
    {
        sdsp::fft_radix4(s);

        double max_error {calc_max_error(S, s)};

        REQUIRE(max_error < 4 * N * std::numeric_limits<double>::epsilon());
    }

    SECTION("FFT Reverse Test")
    {        
        sdsp::fft_radix4<sdsp::reverse_fft>(S);

        double max_error {calc_max_error(S, s)};

        REQUIRE(max_error < 4 * N * std::numeric_limits<double>::epsilon());
    }

    SECTION("FFT Time shift test")
    {
        //s2 shifted 90 degrees from s
        sdsp::complex_array<N> s2 {0};
    
        for (size_t i {0}; i < s2.size(); i++) {
            s2.at(i) = std::cos(n * 2 * M_PI * i / N + (M_PI / 2.0));
        }

        sdsp::fft_radix4(s);
        sdsp::fft_radix4(s2);

        sdsp::complex_array<N> S2 {0};
        S2.at(n) = std::complex<double>(0, N / 2);
        S2.at(N - n) = std::complex<double>(0, -(N / 2.0));

        double max_error {calc_max_error(S2, s2)};
        REQUIRE(max_error < 4 * N * std::numeric_limits<double>::epsilon());
    }
}

TEST_CASE("FFT Radix 4 linearity", "[single-file]")
{
    constexpr double freq1 {1000.0};
    constexpr double freq2 {500.0};
    constexpr double fs {8000.0};
    constexpr double a1 {1.5};
    constexpr double a2 {2.5};
    constexpr uint N {256};

    //test linearity
    //FFT(a1x1[n]+a2x2[n])=a1FFT(x1[n])+a2FFT(x2[n])

    sdsp::complex_array<N> x1 {0};
    sdsp::complex_array<N> x2 {0};
    for (size_t i {0}; i < x1.size(); i++) {
        double value {std::sin(2.0 * M_PI * freq1 * (1.0 / fs) * i)};
        x1.at(i) = std::complex<double>(value, 0);
        
        double value2 {std::sin(2.0 * M_PI * freq2 * (1.0 / fs) * i)};
        x2.at(i) = std::complex<double>(value2, 0);
    }

    sdsp::complex_array<N> a1x1 {x1};
    std::for_each(a1x1.begin(), a1x1.end(), [a1](auto &n) { n *= a1; });

    sdsp::complex_array<N> a2x2 {x2};
    std::for_each(a2x2.begin(), a2x2.end(), [a2](auto &n) { n *= a2; });

    sdsp::complex_array<N> a1x1_plus_a2x2 {0};
    for (size_t i {0}; i < a1x1_plus_a2x2.size(); i++) {
        a1x1_plus_a2x2.at(i) = a1x1.at(i) + a2x2.at(i);
    }

    sdsp::complex_array<N> fft_data1 {a1x1_plus_a2x2};
    sdsp::fft_radix4(fft_data1);

    sdsp::complex_array<N> fft_data2 {x1};
    sdsp::complex_array<N> fft_data3 {x2};
    sdsp::fft_radix4(fft_data2);
    sdsp::fft_radix4(fft_data3);

    sdsp::complex_array<N> a1_fft_data2 {fft_data2};
    std::for_each(a1_fft_data2.begin(), a1_fft_data2.end(), [a1](auto &n) { n *= a1; });

    sdsp::complex_array<N> a2_fft_data3 {fft_data3};
    std::for_each(a2_fft_data3.begin(), a2_fft_data3.end(), [a2](auto &n) { n *= a2; });

    sdsp::complex_array<N> a1_fft_data2_plus_a2_fft_data3 {0};
    for (size_t i = 0; i < a1_fft_data2_plus_a2_fft_data3.size(); i++) {
        a1_fft_data2_plus_a2_fft_data3.at(i) = a1_fft_data2.at(i) + a2_fft_data3.at(i);
    }

    //difference should be very close to 0
    double max_error {calc_max_error(fft_data1, a1_fft_data2_plus_a2_fft_data3)};
    REQUIRE(max_error < 4 * N * std::numeric_limits<double>::epsilon());
}

TEST_CASE("FFT benchmark", "[single-file]")
{
    sdsp::complex_array<1024> complexValues {0.3535, 0.3535, 0.6464, 1.0607, 0.3535, -1.0607, -1.3535, -0.3535};

    BENCHMARK("Radix2 FFT")
    {
        sdsp::fft_radix2(complexValues);
        return complexValues;
    };

    BENCHMARK("Radix4 FFT")
    {
        sdsp::fft_radix4(complexValues);
        return complexValues;
    };
}
