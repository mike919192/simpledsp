
#include "fourier.h"
#include "fft.h"
#include <complex>
#include <algorithm>

int main()
{
    std::array<double, 32> values {0.3535, 0, 0.3535, 0, 0.6464, 0, 1.0607, 0, 0.3535, 0, -1.0607, 0, -1.3535, 0, -0.3535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    four1(values, 1);

    std::array<std::complex<double>, 8> complexValues {0.3535, 0.3535, 0.6464, 1.0607, 0.3535, -1.0607, -1.3535, -0.3535};
    std::array<std::complex<double>, 16> complexValues2 {0.3535, 0.3535, 0.6464, 1.0607, 0.3535, -1.0607, -1.3535, -0.3535, 0, 0, 0, 0, 0, 0, 0, 0};

    std::array<std::complex<double>, 256> complexValues3 {0};
    std::array<double, 256> outputMag {0};
    std::array<double, 256> outputPhase {0};
    const double freq1 = 1000.0;
    const double freq2 = 500.0;
    const double fs = 8000.0;
    const double a1 = 1.5;
    const double a2 = 2.5;

    for (size_t i = 0; i < complexValues3.size(); i++)
	{
		double value1 = std::sin(2.0 * M_PI * freq1 * (1.0/fs) * i);
		double value2 = std::sin(2.0 * M_PI * freq2 * (1.0/fs) * i);

        complexValues3.at(i) = std::complex<double>(value1 + value2, 0);
	}

    fft(complexValues3);

    for (size_t i = 0; i < complexValues3.size(); i++)
	{
		outputMag.at(i) = std::abs(complexValues3.at(i));
        outputPhase.at(i) = std::arg(complexValues3.at(i));
	}

    fft<reverse_fft>(complexValues3);

    //test linearity
    std::array<std::complex<double>, 256> x1 {0};
    std::array<std::complex<double>, 256> x2 {0};
    for (size_t i = 0; i < x1.size(); i++) {
        double value = std::sin(2.0 * M_PI * freq1 * (1.0 / fs) * i);
        x1.at(i) = std::complex<double>(value, 0);
        
        double value2 = std::sin(2.0 * M_PI * freq2 * (1.0 / fs) * i);
        x2.at(i) = std::complex<double>(value2, 0);
    }

    std::array<std::complex<double>, 256> a1x1 = x1;
    std::for_each(a1x1.begin(), a1x1.end(), [a1](auto &n) { n *= a1; });

    std::array<std::complex<double>, 256> a2x2 = x2;
    std::for_each(a2x2.begin(), a2x2.end(), [a2](auto &n) { n *= a2; });

    std::array<std::complex<double>, 256> a1x1_plus_a2x2 {0};
    for (size_t i = 0; i < a1x1_plus_a2x2.size(); i++) {
        a1x1_plus_a2x2.at(i) = a1x1.at(i) + a2x2.at(i);
    }

    std::array<std::complex<double>, 256> fft_data1 = a1x1_plus_a2x2;
    fft(fft_data1);

    std::array<std::complex<double>, 256> fft_data2 = x1;
    std::array<std::complex<double>, 256> fft_data3 = x2;
    fft(fft_data2);
    fft(fft_data3);

    std::array<std::complex<double>, 256> a1_fft_data2 = fft_data2;
    std::for_each(a1_fft_data2.begin(), a1_fft_data2.end(), [a1](auto &n) { n *= a1; });

    std::array<std::complex<double>, 256> a2_fft_data3 = fft_data3;
    std::for_each(a2_fft_data3.begin(), a2_fft_data3.end(), [a2](auto &n) { n *= a2; });

    std::array<std::complex<double>, 256> a1_fft_data2_plus_a2_fft_data3 {0};
    for (size_t i = 0; i < a1_fft_data2_plus_a2_fft_data3.size(); i++) {
        a1_fft_data2_plus_a2_fft_data3.at(i) = a1_fft_data2.at(i) + a2_fft_data3.at(i);
    }

    //difference should be very close to 0
    std::array<double, 256> difference {0};
    for (size_t i = 0; i < difference.size(); i++) {
        difference.at(i) = std::abs(fft_data1.at(i) - a1_fft_data2_plus_a2_fft_data3.at(i));
    }
    double max = *std::max_element(difference.begin(), difference.end());
}
