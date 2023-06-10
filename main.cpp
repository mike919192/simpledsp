
#include "fourier.h"
#include "fft.h"
#include <complex>

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
}
