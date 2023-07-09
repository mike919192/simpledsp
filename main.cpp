
#include "fft.h"
#include <complex>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>

std::string exec(const char* cmd) {
    std::array<char, 1024> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

void binwrite(const std::vector<std::complex<double>>& data)
{
    std::remove("example.bin");
    std::ofstream myfile("example.bin", std::ios::binary);
    for (size_t i = 0; i < data.size(); i++) {
        double real = data.at(i).real();
        double imag = data.at(i).imag();
        myfile.write(reinterpret_cast<char*>(&real), sizeof real);
        myfile.write(reinterpret_cast<char*>(&imag), sizeof imag);
    }
}

std::vector<std::complex<double>> binread()
{
    auto filesize = std::filesystem::file_size("example2.bin");
    size_t numElements = filesize / (sizeof(double) * 2);
    std::vector<std::complex<double>> outVec(numElements);
    std::ifstream myfile("example2.bin", std::ios::binary);
    
    for (size_t i = 0; i < outVec.size(); i++) {
        double real {0};
        double imag {0};        
        myfile.read(reinterpret_cast<char*>(&real), sizeof real);
        myfile.read(reinterpret_cast<char*>(&imag), sizeof imag);
        outVec.at(i) = std::complex(real, imag);
    }
    return outVec;
}

int main()
{
    complex_array<double, 1024> complexValues {0.3535, 0.3535, 0.6464, 1.0607, 0.3535, -1.0607, -1.3535, -0.3535};
    //complex_array<double, 16> complexValues2 {0.3535, 0.3535, 0.6464, 1.0607, 0.3535, -1.0607, -1.3535, -0.3535, 0, 0, 0, 0, 0, 0, 0, 0};
    complex_array<double, 1024> complexValues2 {0.3535, 0.3535, 0.6464, 1.0607, 0.3535, -1.0607, -1.3535, -0.3535};
    fft(complexValues);
    fft_radix4(complexValues2);

    std::vector<std::complex<double>> vec(complexValues.begin(), complexValues.end());
    //csvwrite(vec);
    //std::cout << exec("octave-cli octavefft.m");
    

    complex_array<double, 256> complexValues3 {0};
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

    std::vector<std::complex<double>> vec2(complexValues3.begin(), complexValues3.end());
    binwrite(vec2);
    std::cout << exec("octave-cli octavefft.m");
    std::vector<std::complex<double>> readVec = binread();

    fft(complexValues3);

    //difference should be very close to 0
    std::array<double, 256> difference2 {0};
    for (size_t i = 0; i < difference2.size(); i++) {
        difference2.at(i) = std::abs(complexValues3.at(i) - readVec.at(i));
    }
    double max2 = *std::max_element(difference2.begin(), difference2.end());

    for (size_t i = 0; i < complexValues3.size(); i++)
	{
		outputMag.at(i) = std::abs(complexValues3.at(i));
        outputPhase.at(i) = std::arg(complexValues3.at(i));
	}

    fft<reverse_fft>(complexValues3);

    //test linearity
    complex_array<double, 256> x1 {0};
    complex_array<double, 256> x2 {0};
    for (size_t i = 0; i < x1.size(); i++) {
        double value = std::sin(2.0 * M_PI * freq1 * (1.0 / fs) * i);
        x1.at(i) = std::complex<double>(value, 0);
        
        double value2 = std::sin(2.0 * M_PI * freq2 * (1.0 / fs) * i);
        x2.at(i) = std::complex<double>(value2, 0);
    }

    complex_array<double, 256> a1x1 = x1;
    std::for_each(a1x1.begin(), a1x1.end(), [a1](auto &n) { n *= a1; });

    complex_array<double, 256> a2x2 = x2;
    std::for_each(a2x2.begin(), a2x2.end(), [a2](auto &n) { n *= a2; });

    complex_array<double, 256> a1x1_plus_a2x2 {0};
    for (size_t i = 0; i < a1x1_plus_a2x2.size(); i++) {
        a1x1_plus_a2x2.at(i) = a1x1.at(i) + a2x2.at(i);
    }

    complex_array<double, 256> fft_data1 = a1x1_plus_a2x2;
    fft(fft_data1);

    complex_array<double, 256> fft_data2 = x1;
    complex_array<double, 256> fft_data3 = x2;
    fft(fft_data2);
    fft(fft_data3);

    complex_array<double, 256> a1_fft_data2 = fft_data2;
    std::for_each(a1_fft_data2.begin(), a1_fft_data2.end(), [a1](auto &n) { n *= a1; });

    complex_array<double, 256> a2_fft_data3 = fft_data3;
    std::for_each(a2_fft_data3.begin(), a2_fft_data3.end(), [a2](auto &n) { n *= a2; });

    complex_array<double, 256> a1_fft_data2_plus_a2_fft_data3 {0};
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
