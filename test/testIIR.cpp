
#include "catch2/catch_all.hpp"
#include "sdsp/casc2orderIIR.h"
#include <filesystem>
#include <fstream>

std::tuple<std::vector<double>, sdsp::FilterType, double, double, double>
csvreadImpulse2(const std::string& filename)
{
    std::ifstream myfile;
    myfile.open(filename);
    unsigned int fType{ 0 };
    char comma{ 0 };
    unsigned int n{ 0 };
    sdsp::FilterType fType_out{ sdsp::FilterType::None };
    double fs_out{ 0.0 };
    double f0_out{ 0.0 };
    double Q_out{ 0.0 };
    myfile >> fType >> comma >> fs_out >> comma >> f0_out >> comma >> Q_out >>
      comma >> n >> comma;
    fType_out = static_cast<sdsp::FilterType>(fType);

    std::vector<double> impulse(n);

    for (unsigned int i = 0; i < n - 1; i++) {
        myfile >> impulse.at(i) >> comma;
    }
    myfile >> impulse.back();
    return { impulse, fType_out, fs_out, f0_out, Q_out };
}

TEST_CASE("Filter test")
{
    SECTION("Test impulse response")
    {
        std::string path = "../../../test_data/impulse_response";
        for (const auto& entry : std::filesystem::directory_iterator(path)) {

            auto [readImpulse, fType, fs, f0, Q] =
              csvreadImpulse2(entry.path().string());
            sdsp::casc2orderIIR<4> df;

            if (fType == sdsp::FilterType::LowPass) {
                df.SetLPCoeff(f0, fs);
            } else if (fType == sdsp::FilterType::HighPass) {
                df.SetHPCoeff(f0, fs);
            } else if (fType == sdsp::FilterType::BandPass) {
                df.SetBPCoeff(f0, fs, Q);
            } else {
                throw std::runtime_error("Unknown filter type");
            }
            sdsp::casc2orderIIR<4> df2 = df;
            std::vector<double> data(readImpulse.size());
            data.at(0) = 1.0;
            df.Process(data.begin(), data.end());

            std::vector<double> error(readImpulse.size());

            for (size_t i = 0; i < error.size(); i++) {
                error.at(i) = std::abs(data.at(i) - readImpulse.at(i));
            }
            double maxError = *std::max_element(error.begin(), error.end());
            REQUIRE(maxError < 1e-12); // typically error 1e-16 or less

            // also check that we can process in blocks and get the same result
            std::vector<double> data2(readImpulse.size());
            data2.at(0) = 1.0;

            constexpr unsigned int blockSize{ 32 };
            unsigned int index{ 0 };
            for (index = 0; index <= data2.size() - blockSize;
                 index += blockSize) {
                df2.Process(std::next(data2.begin(), index),
                            std::next(data2.begin(), index + blockSize));
            }

            if (index < data2.size()) {
                df2.Process(std::next(data2.begin(), index), data2.end());
            }

            REQUIRE(data == data2);
        }
    }

    SECTION("Test filter preload")
    {
        constexpr double fs{ 100e3 };
        constexpr double f0{ 10e3 };
        constexpr double Q{ 1.1 };
        constexpr double steadyValue{ 10.0 };

        {
            std::array<double, 1024> steadyLP;
            std::fill(steadyLP.begin(), steadyLP.end(), steadyValue);
            sdsp::casc2orderIIR<4> lpFilter;
            lpFilter.SetLPCoeff(f0, fs);
            lpFilter.PreloadFilter(steadyValue);
            lpFilter.Process(steadyLP.begin(), steadyLP.end());
            auto calcError = [steadyValue](double& n) {
                n = std::abs(n - steadyValue);
            };
            std::for_each(steadyLP.begin(), steadyLP.end(), calcError);
            double maxErrorLP =
              *std::max_element(steadyLP.begin(), steadyLP.end());
            REQUIRE(maxErrorLP < 1e-12);
        }

        {
            std::array<double, 1024> steadyHP;
            std::fill(steadyHP.begin(), steadyHP.end(), steadyValue);
            sdsp::casc2orderIIR<4> hpFilter;
            hpFilter.SetHPCoeff(f0, fs);
            hpFilter.PreloadFilter(steadyValue);
            hpFilter.Process(steadyHP.begin(), steadyHP.end());
            auto calcError = [](double& n) { n = std::abs(n); };
            std::for_each(steadyHP.begin(), steadyHP.end(), calcError);
            double maxErrorHP =
              *std::max_element(steadyHP.begin(), steadyHP.end());
            REQUIRE(maxErrorHP < 1e-12);
        }

        {
            std::array<double, 1024> steadyBP;
            std::fill(steadyBP.begin(), steadyBP.end(), steadyValue);
            sdsp::casc2orderIIR<4> bpFilter;
            bpFilter.SetBPCoeff(f0, fs, Q);
            bpFilter.PreloadFilter(steadyValue);
            bpFilter.Process(steadyBP.begin(), steadyBP.end());
            auto calcError = [](double& n) { n = std::abs(n); };
            std::for_each(steadyBP.begin(), steadyBP.end(), calcError);
            double maxErrorBP =
              *std::max_element(steadyBP.begin(), steadyBP.end());
            REQUIRE(maxErrorBP < 1e-12);
        }
    }
}

TEST_CASE("LP compile time filter test")
{
    SECTION("Test impulse response of LP specialized")
    {
        std::string path = "../../../test_data/impulse_response";
        for (const auto& entry : std::filesystem::directory_iterator(path)) {

            //check that filename starts with LP
            if (entry.path().filename().string().rfind("LP", 0) != 0)
                continue;

            auto [readImpulse, fType, fs, f0, Q] =
              csvreadImpulse2(entry.path().string());
            sdsp::casc_2o_IIR_lp<4> df;

            if (fType == sdsp::FilterType::LowPass) {
                df.set_coeff(f0, fs);
            } else {
                throw std::runtime_error("Unknown filter type");
            }
            sdsp::casc_2o_IIR_lp<4> df2 = df;
            std::vector<double> data(readImpulse.size());
            data.at(0) = 1.0;
            df.process(data.begin(), data.end());

            std::vector<double> error(readImpulse.size());

            for (size_t i = 0; i < error.size(); i++) {
                error.at(i) = std::abs(data.at(i) - readImpulse.at(i));
            }
            double maxError = *std::max_element(error.begin(), error.end());
            REQUIRE(maxError < 1e-12); // typically error 1e-16 or less

            // also check that we can process in blocks and get the same result
            std::vector<double> data2(readImpulse.size());
            data2.at(0) = 1.0;

            constexpr unsigned int blockSize{ 32 };
            unsigned int index{ 0 };
            for (index = 0; index <= data2.size() - blockSize;
                 index += blockSize) {
                df2.process(std::next(data2.begin(), index),
                            std::next(data2.begin(), index + blockSize));
            }

            if (index < data2.size()) {
                df2.process(std::next(data2.begin(), index), data2.end());
            }

            REQUIRE(data == data2);
        }
    }

    SECTION("Test the LP filter gain")
    {
        constexpr double fs{ 100e3 };
        constexpr double f0{ 10e3 };

        std::array<double, 1024> impulse1{ 0 };
        impulse1.at(0) = 1.0;
        std::array<double, 1024> impulse2{ 0 };
        impulse2.at(0) = 1.0;

        sdsp::casc_2o_IIR_lp<4> df;
        df.set_coeff(f0, fs);

        sdsp::casc_2o_IIR_lp<4> df2;
        df2.set_coeff(f0, fs, 2.0);

        df.process(impulse1.begin(), impulse1.end());
        df2.process(impulse2.begin(), impulse2.end());

        auto times_two = [](double& n) { n = 2.0 * n; };
        std::for_each(impulse1.begin(), impulse1.end(), times_two);

        std::array<double, 1024> error{ 0 };

        for (size_t i = 0; i < error.size(); i++) {
            error.at(i) = std::abs(impulse1.at(i) - impulse2.at(i));
        }
        double maxError = *std::max_element(error.begin(), error.end());
        REQUIRE(maxError < 1e-12);
    }
}

TEST_CASE("HP compile time filter test")
{
    SECTION("Test impulse response of HP specialized")
    {
        std::string path = "../../../test_data/impulse_response";
        for (const auto& entry : std::filesystem::directory_iterator(path)) {

            //check that filename starts with HP
            if (entry.path().filename().string().rfind("HP", 0) != 0)
                continue;

            auto [readImpulse, fType, fs, f0, Q] =
              csvreadImpulse2(entry.path().string());
            sdsp::casc_2o_IIR_hp<4> df;

            if (fType == sdsp::FilterType::HighPass) {
                df.set_coeff(f0, fs);
            } else {
                throw std::runtime_error("Unknown filter type");
            }
            sdsp::casc_2o_IIR_hp<4> df2 = df;
            std::vector<double> data(readImpulse.size());
            data.at(0) = 1.0;
            df.process(data.begin(), data.end());

            std::vector<double> error(readImpulse.size());

            for (size_t i = 0; i < error.size(); i++) {
                error.at(i) = std::abs(data.at(i) - readImpulse.at(i));
            }
            double maxError = *std::max_element(error.begin(), error.end());
            REQUIRE(maxError < 1e-12); // typically error 1e-16 or less

            // also check that we can process in blocks and get the same result
            std::vector<double> data2(readImpulse.size());
            data2.at(0) = 1.0;

            constexpr unsigned int blockSize{ 32 };
            unsigned int index{ 0 };
            for (index = 0; index <= data2.size() - blockSize;
                 index += blockSize) {
                df2.process(std::next(data2.begin(), index),
                            std::next(data2.begin(), index + blockSize));
            }

            if (index < data2.size()) {
                df2.process(std::next(data2.begin(), index), data2.end());
            }

            REQUIRE(data == data2);
        }
    }

    SECTION("Test the HP filter gain")
    {
        constexpr double fs{ 100e3 };
        constexpr double f0{ 10e3 };

        std::array<double, 1024> impulse1{ 0 };
        impulse1.at(0) = 1.0;
        std::array<double, 1024> impulse2{ 0 };
        impulse2.at(0) = 1.0;

        sdsp::casc_2o_IIR_hp<4> df;
        df.set_coeff(f0, fs);

        sdsp::casc_2o_IIR_hp<4> df2;
        df2.set_coeff(f0, fs, 2.0);

        df.process(impulse1.begin(), impulse1.end());
        df2.process(impulse2.begin(), impulse2.end());

        auto times_two = [](double& n) { n = 2.0 * n; };
        std::for_each(impulse1.begin(), impulse1.end(), times_two);

        std::array<double, 1024> error{ 0 };

        for (size_t i = 0; i < error.size(); i++) {
            error.at(i) = std::abs(impulse1.at(i) - impulse2.at(i));
        }
        double maxError = *std::max_element(error.begin(), error.end());
        REQUIRE(maxError < 1e-12);
    }
}

TEST_CASE("BP compile time filter test")
{
    SECTION("Test impulse response of BP specialized")
    {
        std::string path = "../../../test_data/impulse_response";
        for (const auto& entry : std::filesystem::directory_iterator(path)) {

            //check that filename starts with BP
            if (entry.path().filename().string().rfind("BP", 0) != 0)
                continue;

            auto [readImpulse, fType, fs, f0, Q] =
              csvreadImpulse2(entry.path().string());
            sdsp::casc_2o_IIR_bp<4> df;

            if (fType == sdsp::FilterType::BandPass) {
                df.set_coeff(f0, fs, Q);
            } else {
                throw std::runtime_error("Unknown filter type");
            }
            sdsp::casc_2o_IIR_bp<4> df2 = df;
            std::vector<double> data(readImpulse.size());
            data.at(0) = 1.0;
            df.process(data.begin(), data.end());

            std::vector<double> error(readImpulse.size());

            for (size_t i = 0; i < error.size(); i++) {
                error.at(i) = std::abs(data.at(i) - readImpulse.at(i));
            }
            double maxError = *std::max_element(error.begin(), error.end());
            REQUIRE(maxError < 1e-12); // typically error 1e-16 or less

            // also check that we can process in blocks and get the same result
            std::vector<double> data2(readImpulse.size());
            data2.at(0) = 1.0;

            constexpr unsigned int blockSize{ 32 };
            unsigned int index{ 0 };
            for (index = 0; index <= data2.size() - blockSize;
                 index += blockSize) {
                df2.process(std::next(data2.begin(), index),
                            std::next(data2.begin(), index + blockSize));
            }

            if (index < data2.size()) {
                df2.process(std::next(data2.begin(), index), data2.end());
            }

            REQUIRE(data == data2);
        }
    }

    SECTION("Test the BP filter gain")
    {
        constexpr double fs{ 100e3 };
        constexpr double f0{ 10e3 };
        constexpr double Q{ 1.1 };

        std::array<double, 1024> impulse1{ 0 };
        impulse1.at(0) = 1.0;
        std::array<double, 1024> impulse2{ 0 };
        impulse2.at(0) = 1.0;

        sdsp::casc_2o_IIR_bp<4> df;
        df.set_coeff(f0, fs, Q);

        sdsp::casc_2o_IIR_bp<4> df2;
        df2.set_coeff(f0, fs, Q, 2.0);

        df.process(impulse1.begin(), impulse1.end());
        df2.process(impulse2.begin(), impulse2.end());

        auto times_two = [](double& n) { n = 2.0 * n; };
        std::for_each(impulse1.begin(), impulse1.end(), times_two);

        std::array<double, 1024> error{ 0 };

        for (size_t i = 0; i < error.size(); i++) {
            error.at(i) = std::abs(impulse1.at(i) - impulse2.at(i));
        }
        double maxError = *std::max_element(error.begin(), error.end());
        REQUIRE(maxError < 1e-12);
    }
}

TEST_CASE("Filter benchmarks")
{
    SECTION("LP Filter benchmark section")
    {
        constexpr double fs{ 100e3 };
        constexpr double f0{ 10e3 };
        // constexpr double Q {1.1};

        sdsp::casc2orderIIR<4> df;
        df.SetLPCoeff(f0, fs);

        sdsp::casc_2o_IIR_lp<4> df2;
        df2.set_coeff(f0, fs);

        std::array<double, 4096> data{ 0.0 };
        data.at(0) = 1.0;

        BENCHMARK("Runtime configurable LP filter benchmark")
        {
            df.Process(data.begin(), data.end());
            return data;
        };

        BENCHMARK("Specialized LP filter benchmark")
        {
            df2.process(data.begin(), data.end());
            return data;
        };
        SUCCEED();
    }

    SECTION("HP Filter benchmark section")
    {
        constexpr double fs{ 100e3 };
        constexpr double f0{ 10e3 };
        // constexpr double Q {1.1};

        sdsp::casc2orderIIR<4> df;
        df.SetHPCoeff(f0, fs);

        sdsp::casc_2o_IIR_hp<4> df2;
        df2.set_coeff(f0, fs);

        std::array<double, 4096> data{ 0.0 };
        data.at(0) = 1.0;

        BENCHMARK("Runtime configurable HP filter benchmark")
        {
            df.Process(data.begin(), data.end());
            return data;
        };

        BENCHMARK("Specialized HP filter benchmark")
        {
            df2.process(data.begin(), data.end());
            return data;
        };
        SUCCEED();
    }

    SECTION("BP Filter benchmark section")
    {
        constexpr double fs{ 100e3 };
        constexpr double f0{ 10e3 };
        constexpr double Q{ 1.1 };

        sdsp::casc2orderIIR<4> df;
        df.SetBPCoeff(f0, fs, Q);

        sdsp::casc_2o_IIR_bp<4> df2;
        df2.set_coeff(f0, fs, Q);

        std::array<double, 4096> data{ 0.0 };
        data.at(0) = 1.0;

        BENCHMARK("Runtime configurable BP filter benchmark")
        {
            df.Process(data.begin(), data.end());
            return data;
        };

        BENCHMARK("Specialized BP filter benchmark")
        {
            df2.process(data.begin(), data.end());
            return data;
        };
        SUCCEED();
    }
}
