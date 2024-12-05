
#include "catch2/catch_all.hpp"
#include "sdsp/casc_2o_iir.h"
#include <filesystem>
#include <fstream>

std::tuple<std::vector<double>, sdsp::filter_type, double, double, double> csvreadImpulse2(const std::string &filename)
{
    std::ifstream myfile;
    myfile.open(filename);
    unsigned int fType{ 0 };
    char comma{ 0 };
    unsigned int n{ 0 };
    sdsp::filter_type fType_out{ sdsp::filter_type::none };
    double fs_out{ 0.0 };
    double f0_out{ 0.0 };
    double Q_out{ 0.0 };
    myfile >> fType >> comma >> fs_out >> comma >> f0_out >> comma >> Q_out >> comma >> n >> comma;
    fType_out = static_cast<sdsp::filter_type>(fType);

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
        for (const auto &entry : std::filesystem::directory_iterator(path)) {
            auto [readImpulse, fType, fs, f0, Q] = csvreadImpulse2(entry.path().string());
            sdsp::casc_2o_iir<4> df;

            if (fType == sdsp::filter_type::low_pass) {
                df.set_lp_coeff(f0, fs);
            } else if (fType == sdsp::filter_type::high_pass) {
                df.set_hp_coeff(f0, fs);
            } else if (fType == sdsp::filter_type::band_pass) {
                df.set_bp_coeff(f0, fs, Q);
            } else {
                throw std::runtime_error("Unknown filter type");
            }
            sdsp::casc_2o_iir<4> df2 = df;
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
            for (index = 0; index <= data2.size() - blockSize; index += blockSize) {
                df2.process(std::next(data2.begin(), index), std::next(data2.begin(), index + blockSize));
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

        sdsp::casc_2o_iir<4> df;
        df.set_lp_coeff(f0, fs);

        sdsp::casc_2o_iir<4> df2;
        df2.set_lp_coeff(f0, fs, 2.0);

        df.process(impulse1.begin(), impulse1.end());
        df2.process(impulse2.begin(), impulse2.end());

        auto times_two = [](double &n) { n = 2.0 * n; };
        std::for_each(impulse1.begin(), impulse1.end(), times_two);

        std::array<double, 1024> error{ 0 };

        for (size_t i = 0; i < error.size(); i++) {
            error.at(i) = std::abs(impulse1.at(i) - impulse2.at(i));
        }
        double maxError = *std::max_element(error.begin(), error.end());
        REQUIRE(maxError < 1e-12);
    }

    SECTION("Test the HP filter gain")
    {
        constexpr double fs{ 100e3 };
        constexpr double f0{ 10e3 };

        std::array<double, 1024> impulse1{ 0 };
        impulse1.at(0) = 1.0;
        std::array<double, 1024> impulse2{ 0 };
        impulse2.at(0) = 1.0;

        sdsp::casc_2o_iir<4> df;
        df.set_hp_coeff(f0, fs);

        sdsp::casc_2o_iir<4> df2;
        df2.set_hp_coeff(f0, fs, 2.0);

        df.process(impulse1.begin(), impulse1.end());
        df2.process(impulse2.begin(), impulse2.end());

        auto times_two = [](double &n) { n = 2.0 * n; };
        std::for_each(impulse1.begin(), impulse1.end(), times_two);

        std::array<double, 1024> error{ 0 };

        for (size_t i = 0; i < error.size(); i++) {
            error.at(i) = std::abs(impulse1.at(i) - impulse2.at(i));
        }
        double maxError = *std::max_element(error.begin(), error.end());
        REQUIRE(maxError < 1e-12);
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

        sdsp::casc_2o_iir<4> df;
        df.set_bp_coeff(f0, fs, Q);

        sdsp::casc_2o_iir<4> df2;
        df2.set_bp_coeff(f0, fs, Q, 2.0);

        df.process(impulse1.begin(), impulse1.end());
        df2.process(impulse2.begin(), impulse2.end());

        auto times_two = [](double &n) { n = 2.0 * n; };
        std::for_each(impulse1.begin(), impulse1.end(), times_two);

        std::array<double, 1024> error{ 0 };

        for (size_t i = 0; i < error.size(); i++) {
            error.at(i) = std::abs(impulse1.at(i) - impulse2.at(i));
        }
        double maxError = *std::max_element(error.begin(), error.end());
        REQUIRE(maxError < 1e-12);
    }

    SECTION("Test filter preload")
    {
        constexpr double fs{ 100e3 };
        constexpr double f0{ 10e3 };
        constexpr double Q{ 1.1 };
        constexpr double steadyValue{ 10.0 };

        {
            std::array<double, 1024> steadyLP;
            steadyLP.fill(steadyValue);
            sdsp::casc_2o_iir<4> lpFilter;
            lpFilter.set_lp_coeff(f0, fs);
            lpFilter.preload_filter(steadyValue);
            lpFilter.process(steadyLP.begin(), steadyLP.end());
            auto calcError = [steadyValue](double &n) { n = std::abs(n - steadyValue); };
            std::for_each(steadyLP.begin(), steadyLP.end(), calcError);
            double maxErrorLP = *std::max_element(steadyLP.begin(), steadyLP.end());
            REQUIRE(maxErrorLP < 1e-12);
        }

        {
            std::array<double, 1024> steadyHP;
            steadyHP.fill(steadyValue);
            sdsp::casc_2o_iir<4> hpFilter;
            hpFilter.set_hp_coeff(f0, fs);
            hpFilter.preload_filter(steadyValue);
            hpFilter.process(steadyHP.begin(), steadyHP.end());
            auto calcError = [](double &n) { n = std::abs(n); };
            std::for_each(steadyHP.begin(), steadyHP.end(), calcError);
            double maxErrorHP = *std::max_element(steadyHP.begin(), steadyHP.end());
            REQUIRE(maxErrorHP < 1e-12);
        }

        {
            std::array<double, 1024> steadyBP;
            steadyBP.fill(steadyValue);
            sdsp::casc_2o_iir<4> bpFilter;
            bpFilter.set_bp_coeff(f0, fs, Q);
            bpFilter.preload_filter(steadyValue);
            bpFilter.process(steadyBP.begin(), steadyBP.end());
            auto calcError = [](double &n) { n = std::abs(n); };
            std::for_each(steadyBP.begin(), steadyBP.end(), calcError);
            double maxErrorBP = *std::max_element(steadyBP.begin(), steadyBP.end());
            REQUIRE(maxErrorBP < 1e-12);
        }
    }
}

TEST_CASE("LP compile time filter test")
{
    SECTION("Test impulse response of LP specialized")
    {
        std::string path = "../../../test_data/impulse_response";
        for (const auto &entry : std::filesystem::directory_iterator(path)) {
            //check that filename starts with LP
            if (entry.path().filename().string().rfind("LP", 0) != 0)
                continue;

            auto [readImpulse, fType, fs, f0, Q] = csvreadImpulse2(entry.path().string());
            sdsp::casc_2o_iir_lp<4> df;

            if (fType == sdsp::filter_type::low_pass) {
                df.set_lp_coeff(f0, fs);
            } else {
                throw std::runtime_error("Unknown filter type");
            }
            sdsp::casc_2o_iir_lp<4> df2 = df;
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
            for (index = 0; index <= data2.size() - blockSize; index += blockSize) {
                df2.process(std::next(data2.begin(), index), std::next(data2.begin(), index + blockSize));
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

        sdsp::casc_2o_iir_lp<4> df;
        df.set_lp_coeff(f0, fs);

        sdsp::casc_2o_iir_lp<4> df2;
        df2.set_lp_coeff(f0, fs, 2.0);

        df.process(impulse1.begin(), impulse1.end());
        df2.process(impulse2.begin(), impulse2.end());

        auto times_two = [](double &n) { n = 2.0 * n; };
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
        for (const auto &entry : std::filesystem::directory_iterator(path)) {
            //check that filename starts with HP
            if (entry.path().filename().string().rfind("HP", 0) != 0)
                continue;

            auto [readImpulse, fType, fs, f0, Q] = csvreadImpulse2(entry.path().string());
            sdsp::casc_2o_iir_hp<4> df;

            if (fType == sdsp::filter_type::high_pass) {
                df.set_hp_coeff(f0, fs);
            } else {
                throw std::runtime_error("Unknown filter type");
            }
            sdsp::casc_2o_iir_hp<4> df2 = df;
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
            for (index = 0; index <= data2.size() - blockSize; index += blockSize) {
                df2.process(std::next(data2.begin(), index), std::next(data2.begin(), index + blockSize));
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

        sdsp::casc_2o_iir_hp<4> df;
        df.set_hp_coeff(f0, fs);

        sdsp::casc_2o_iir_hp<4> df2;
        df2.set_hp_coeff(f0, fs, 2.0);

        df.process(impulse1.begin(), impulse1.end());
        df2.process(impulse2.begin(), impulse2.end());

        auto times_two = [](double &n) { n = 2.0 * n; };
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
        for (const auto &entry : std::filesystem::directory_iterator(path)) {
            //check that filename starts with BP
            if (entry.path().filename().string().rfind("BP", 0) != 0)
                continue;

            auto [readImpulse, fType, fs, f0, Q] = csvreadImpulse2(entry.path().string());
            sdsp::casc_2o_iir_bp<4> df;

            if (fType == sdsp::filter_type::band_pass) {
                df.set_bp_coeff(f0, fs, Q);
            } else {
                throw std::runtime_error("Unknown filter type");
            }
            sdsp::casc_2o_iir_bp<4> df2 = df;
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
            for (index = 0; index <= data2.size() - blockSize; index += blockSize) {
                df2.process(std::next(data2.begin(), index), std::next(data2.begin(), index + blockSize));
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

        sdsp::casc_2o_iir_bp<4> df;
        df.set_bp_coeff(f0, fs, Q);

        sdsp::casc_2o_iir_bp<4> df2;
        df2.set_bp_coeff(f0, fs, Q, 2.0);

        df.process(impulse1.begin(), impulse1.end());
        df2.process(impulse2.begin(), impulse2.end());

        auto times_two = [](double &n) { n = 2.0 * n; };
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

        sdsp::casc_2o_iir<4> df;
        df.set_lp_coeff(f0, fs);

        sdsp::casc_2o_iir_lp<4> df2;
        df2.set_lp_coeff(f0, fs);

        std::array<double, 4096> data{ 0.0 };
        data.at(0) = 1.0;

        BENCHMARK("Runtime configurable LP filter benchmark")
        {
            std::array<double, 4096> data2 = data;
            df.process(data2.begin(), data2.end());
            return data2;
        };

        BENCHMARK("Specialized LP filter benchmark")
        {
            std::array<double, 4096> data2 = data;
            df2.process(data2.begin(), data2.end());
            return data2;
        };
        SUCCEED();
    }

    SECTION("HP Filter benchmark section")
    {
        constexpr double fs{ 100e3 };
        constexpr double f0{ 10e3 };
        // constexpr double Q {1.1};

        sdsp::casc_2o_iir<4> df;
        df.set_hp_coeff(f0, fs);

        sdsp::casc_2o_iir_hp<4> df2;
        df2.set_hp_coeff(f0, fs);

        std::array<double, 4096> data{ 0.0 };
        data.at(0) = 1.0;

        BENCHMARK("Runtime configurable HP filter benchmark")
        {
            std::array<double, 4096> data2 = data;
            df.process(data2.begin(), data2.end());
            return data2;
        };

        BENCHMARK("Specialized HP filter benchmark")
        {
            std::array<double, 4096> data2 = data;
            df2.process(data2.begin(), data2.end());
            return data2;
        };
        SUCCEED();
    }

    SECTION("BP Filter benchmark section")
    {
        constexpr double fs{ 100e3 };
        constexpr double f0{ 10e3 };
        constexpr double Q{ 1.1 };

        sdsp::casc_2o_iir<4> df;
        df.set_bp_coeff(f0, fs, Q);

        sdsp::casc_2o_iir_bp<4> df2;
        df2.set_bp_coeff(f0, fs, Q);

        std::array<double, 4096> data{ 0.0 };
        data.at(0) = 1.0;

        BENCHMARK("Runtime configurable BP filter benchmark")
        {
            std::array<double, 4096> data2 = data;
            df.process(data2.begin(), data2.end());
            return data2;
        };

        BENCHMARK("Specialized BP filter benchmark")
        {
            std::array<double, 4096> data2 = data;
            df2.process(data2.begin(), data2.end());
            return data2;
        };
        SUCCEED();
    }
}
