
#include "fourier.h"
#include "fft.h"
#include <complex>

int main()
{
    std::array<double, 32> values {0.3535, 0, 0.3535, 0, 0.6464, 0, 1.0607, 0, 0.3535, 0, -1.0607, 0, -1.3535, 0, -0.3535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    four1(values, 1);

    std::array<std::complex<double>, 8> complexValues {0.3535, 0.3535, 0.6464, 1.0607, 0.3535, -1.0607, -1.3535, -0.3535};
    std::array<std::complex<double>, 16> complexValues2 {0.3535, 0.3535, 0.6464, 1.0607, 0.3535, -1.0607, -1.3535, -0.3535, 0, 0, 0, 0, 0, 0, 0, 0};

    fft(complexValues2);
}
