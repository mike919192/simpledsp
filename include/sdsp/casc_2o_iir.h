#pragma once
#include "filter_type.h"
#include <array>
#include <cmath>

namespace sdsp
{
template <size_t m_t>
class casc_2o_iir {
private:
    int m_pos{ 0 };

    double m_gain{ 1.0 };

    std::array<std::array<double, 3>, m_t + 1> m_mem{ 0 };

    std::array<std::array<double, 3>, m_t> m_b_coeff{ 0 };
    std::array<std::array<double, 3>, m_t> m_a_coeff{ 0 };

    filter_type m_f_type{ filter_type::none };

public:
    casc_2o_iir()
    {
        static_assert(m_t % 2 == 0, "M must be even!");
    }

    void copy_coeff_from(const casc_2o_iir<m_t> &other_filter)
    {
        m_gain = other_filter.m_gain;
        m_b_coeff = other_filter.m_b_coeff;
        m_a_coeff = other_filter.m_a_coeff;
        m_f_type = other_filter.m_f_type;
    }

    template <typename iter_t>
    void process(iter_t begin, iter_t end)
    {
        constexpr int order{ 2 };
        constexpr int j1{ 1 };
        constexpr int j2{ 2 };

        int p{ m_pos };

        std::array<std::array<double, 3>, m_t + 1> y = m_mem;

        std::array<std::array<double, 3>, m_t> b = m_b_coeff;

        std::array<std::array<double, 3>, m_t> a = m_a_coeff;

        while (begin < end) {
            y.at(0).at(p) = *begin * m_gain;

            int d1{ p - j1 };
            if (d1 < 0)
                d1 += order + 1;

            int d2{ p - j2 };
            if (d2 < 0)
                d2 += order + 1;

            uint j;

            for (j = 0; j < m_t; j++) {
                y.at(j + 1).at(p) = y.at(j).at(p);

                y.at(j + 1).at(p) += y.at(j).at(d1) * b.at(j).at(j1) - y.at(j + 1).at(d1) * a.at(j).at(j1);
                y.at(j + 1).at(p) += y.at(j).at(d2) * b.at(j).at(j2) - y.at(j + 1).at(d2) * a.at(j).at(j2);
            }

            *begin = y.at(j).at(p);

            p++;
            if (p > order)
                p = 0;
            begin++;
        }
        m_pos = p;
        m_mem = y;
    }

    void set_bp_coeff(double f0, double fs, double q, double gain_in = 1.0)
    {
        m_gain = gain_in;
        double q2{ 2 * q };
        m_f_type = filter_type::band_pass;

        double e0{ 2 * M_PI * f0 / fs };
        double dnm{ std::sin(e0) };

        double de{ 2 * std::tan(e0 / q2) / dnm };

        for (unsigned int k = 0; k < m_t / 2; k++) {
            double d{ 2 * std::sin((2 * k + 1) * M_PI / (2.0 * m_t)) };

            double a{ (1 + de * de / 4.0) * 2 / d / de };
            double dk{ std::sqrt(de * d / (a + std::sqrt(a * a - 1))) };

            double b{ d * de / dk / 2.0 };
            double w{ b + std::sqrt(b * b - 1) };

            double t{ std::tan(e0 / 2.0) };

            double e1{ 2.0 * std::atan(t / w) };
            double e2{ 2.0 * std::atan(w * t) };

            t = dk * std::sin(e1) / 2.0;
            dnm = (1 + t);
            double beta1{ (1 - t) / dnm / 2.0 };

            t = dk * std::sin(e2) / 2.0;
            dnm = (1 + t);
            double beta2{ (1 - t) / dnm / 2.0 };

            double gamma1{ (0.5 + beta1) * std::cos(e1) };
            double gamma2{ (0.5 + beta2) * std::cos(e2) };

            t = std::sqrt(1 + (w - 1 / w) / dk * (w - 1 / w) / dk);
            double alpha1{ (0.5 - beta1) * t / 2.0 };
            double alpha2{ (0.5 - beta2) * t / 2.0 };

            m_gain *= 4 * alpha1 * alpha2;

            m_b_coeff.at(2 * k).at(0) = 1.0;
            m_b_coeff.at(2 * k + 1).at(0) = 1.0;
            m_b_coeff.at(2 * k).at(1) = 0;
            m_b_coeff.at(2 * k + 1).at(1) = 0;
            m_b_coeff.at(2 * k).at(2) = -1.0;
            m_b_coeff.at(2 * k + 1).at(2) = -1.0;

            m_a_coeff.at(2 * k).at(0) = 1;
            m_a_coeff.at(2 * k + 1).at(0) = 1;
            m_a_coeff.at(2 * k).at(1) = -2 * gamma1;
            m_a_coeff.at(2 * k + 1).at(1) = -2 * gamma2;
            m_a_coeff.at(2 * k).at(2) = 2 * beta1;
            m_a_coeff.at(2 * k + 1).at(2) = 2 * beta2;
        }
    }

    void set_hp_coeff(double f0, double fs, double gain_in = 1.0)
    {
        m_gain = gain_in;
        m_f_type = filter_type::high_pass;

        double e0{ 2 * M_PI * f0 / fs };
        for (unsigned int k = 0; k < m_t; k++) {
            double dk{ 2 * std::sin((2 * k + 1) * M_PI / (4.0 * m_t)) };

            double t{ dk * std::sin(e0) / 2 };
            double dnm{ 1 + t };

            double beta1{ (1 - t) / dnm / 2 };
            double gamma1{ (0.5 + beta1) * std::cos(e0) };
            double alpha1{ (0.5 + beta1 + gamma1) / 4 };

            m_gain *= 2 * alpha1;

            m_b_coeff.at(k).at(0) = 1.0;
            m_b_coeff.at(k).at(1) = -2.0;
            m_b_coeff.at(k).at(2) = 1.0;

            m_a_coeff.at(k).at(0) = 1;
            m_a_coeff.at(k).at(1) = -2 * gamma1;
            m_a_coeff.at(k).at(2) = 2 * beta1;
        }
    }

    void set_lp_coeff(double f0, double fs, double gain_in = 1.0)
    {
        m_gain = gain_in;
        m_f_type = filter_type::low_pass;

        double e0{ 2 * M_PI * f0 / fs };
        for (unsigned int k = 0; k < m_t; k++) {
            double dk{ 2 * std::sin((2 * k + 1) * M_PI / (4.0 * m_t)) };

            double t{ dk * std::sin(e0) / 2 };
            double dnm{ 1 + t };

            double beta1{ (1 - t) / dnm / 2 };
            double gamma1{ (0.5 + beta1) * std::cos(e0) };
            double alpha1{ (0.5 + beta1 - gamma1) / 4 };

            m_gain *= 2 * alpha1;

            m_b_coeff.at(k).at(0) = 1.0;
            m_b_coeff.at(k).at(1) = 2.0;
            m_b_coeff.at(k).at(2) = 1.0;

            m_a_coeff.at(k).at(0) = 1;
            m_a_coeff.at(k).at(1) = -2 * gamma1;
            m_a_coeff.at(k).at(2) = 2 * beta1;
        }
    }

    // preload the filter memory for steady state input equal to value parameter
    void preload_filter(double value)
    {
        double preload_value = value * m_gain;
        std::array<std::array<double, 3>, m_t + 1> mem_vals{ 0 };
        for (int i = 0; i < 3; i++) {
            mem_vals.at(0).at(i) = preload_value;
        }
        if (m_f_type == filter_type::low_pass) {
            for (uint j = 1; j < m_t + 1; j++) {
                preload_value /= 1 + m_a_coeff.at(j - 1).at(1) + m_a_coeff.at(j - 1).at(2);
                preload_value *= m_b_coeff.at(j - 1).at(0) + m_b_coeff.at(j - 1).at(1) + m_b_coeff.at(j - 1).at(2);
                for (uint i = 0; i < 3; i++) {
                    mem_vals.at(j).at(i) = preload_value;
                }
            }
        }
        m_mem = mem_vals;
    }
};

template <size_t m_t>
class casc_2o_iir_base {
protected:
    int m_pos{ 0 };

    double m_gain{ 1.0 };

    std::array<std::array<double, 3>, m_t + 1> m_mem{ 0 };

    std::array<std::array<double, 3>, m_t> m_a_coeff{ 0 };

    template <typename iir_t, typename iter_t>
    void process_base(iter_t begin, iter_t end)
    {
        constexpr int order{ 2 };
        constexpr int j1{ 1 };
        constexpr int j2{ 2 };

        int p{ m_pos };

        std::array<std::array<double, 3>, m_t + 1> y = m_mem;

        std::array<std::array<double, 3>, m_t> a = m_a_coeff;

        while (begin < end) {
            y.at(0).at(p) = *begin * m_gain;

            int d1{ p - j1 };
            if (d1 < 0)
                d1 += order + 1;

            int d2{ p - j2 };
            if (d2 < 0)
                d2 += order + 1;

            iir_t::process_spec(y, a, p, d1, d2);

            *begin = y.at(m_t).at(p);

            p++;
            if (p > order)
                p = 0;
            begin++;
        }
        m_pos = p;
        m_mem = y;
    }
};

template <size_t m_t>
class casc_2o_iir_lp : casc_2o_iir_base<m_t> {
public:
    casc_2o_iir_lp()
    {
        static_assert(m_t % 2 == 0, "M must be even!");
    }

    void copy_coeff_from(const casc_2o_iir_lp<m_t> &other_filter)
    {
        this->gain = other_filter.gain;
        this->aCoeff = other_filter.aCoeff;
    }

    template <typename iter_t>
    void process(iter_t begin, iter_t end)
    {
        this->template process_base<casc_2o_iir_lp<m_t>>(begin, end);
    }

    static void process_spec(std::array<std::array<double, 3>, m_t + 1> &y,
                             const std::array<std::array<double, 3>, m_t> &a, int p, int d1, int d2)
    {
        for (uint j = 0; j < m_t; j++) {
            y.at(j + 1).at(p) = y.at(j).at(p);

            y.at(j + 1).at(p) += y.at(j).at(d1) + y.at(j).at(d1) - y.at(j + 1).at(d1) * a.at(j).at(1);
            y.at(j + 1).at(p) += y.at(j).at(d2) - y.at(j + 1).at(d2) * a.at(j).at(2);
        }
    }

    void set_lp_coeff(double f0, double fs, double gain_in = 1.0)
    {
        this->m_gain = gain_in;

        double e0{ 2 * M_PI * f0 / fs };
        for (unsigned int k = 0; k < m_t; k++) {
            double dk{ 2 * std::sin((2 * k + 1) * M_PI / (4.0 * m_t)) };

            double t{ dk * std::sin(e0) / 2 };
            double dnm{ 1 + t };

            double beta1{ (1 - t) / dnm / 2 };
            double gamma1{ (0.5 + beta1) * std::cos(e0) };
            double alpha1{ (0.5 + beta1 - gamma1) / 4 };

            this->m_gain *= 2 * alpha1;
            // bCoeff.at(k).at(0) = 1.0;
            // bCoeff.at(k).at(1) = 2.0;
            // bCoeff.at(k).at(2) = 1.0;

            this->m_a_coeff.at(k).at(0) = 1;
            this->m_a_coeff.at(k).at(1) = -2 * gamma1;
            this->m_a_coeff.at(k).at(2) = 2 * beta1;
        }
    }
};

template <size_t m_t>
class casc_2o_iir_hp : casc_2o_iir_base<m_t> {
public:
    casc_2o_iir_hp()
    {
        static_assert(m_t % 2 == 0, "M must be even!");
    }

    void copy_coeff_from(const casc_2o_iir_hp<m_t> &other_filter)
    {
        this->gain = other_filter.gain;
        this->aCoeff = other_filter.aCoeff;
    }

    template <typename iter_t>
    void process(iter_t begin, iter_t end)
    {
        this->template process_base<casc_2o_iir_hp<m_t>>(begin, end);
    }

    static void process_spec(std::array<std::array<double, 3>, m_t + 1> &y,
                             const std::array<std::array<double, 3>, m_t> &a, int p, int d1, int d2)
    {
        for (uint j = 0; j < m_t; j++) {
            y.at(j + 1).at(p) = y.at(j).at(p);

            y.at(j + 1).at(p) += -y.at(j).at(d1) - y.at(j).at(d1) - y.at(j + 1).at(d1) * a.at(j).at(1);
            y.at(j + 1).at(p) += y.at(j).at(d2) - y.at(j + 1).at(d2) * a.at(j).at(2);
        }
    }

    void set_hp_coeff(double f0, double fs, double gain_in = 1.0)
    {
        this->m_gain = gain_in;

        double e0{ 2 * M_PI * f0 / fs };
        for (unsigned int k = 0; k < m_t; k++) {
            double dk{ 2 * std::sin((2 * k + 1) * M_PI / (4.0 * m_t)) };

            double t{ dk * std::sin(e0) / 2 };
            double dnm{ 1 + t };

            double beta1{ (1 - t) / dnm / 2 };
            double gamma1{ (0.5 + beta1) * std::cos(e0) };
            double alpha1{ (0.5 + beta1 + gamma1) / 4 };

            this->m_gain *= 2 * alpha1;
            // bCoeff.at(k).at(0) = 1.0;
            // bCoeff.at(k).at(1) = -2.0;
            // bCoeff.at(k).at(2) = 1.0;

            this->m_a_coeff.at(k).at(0) = 1;
            this->m_a_coeff.at(k).at(1) = -2 * gamma1;
            this->m_a_coeff.at(k).at(2) = 2 * beta1;
        }
    }
};

template <size_t m_t>
class casc_2o_iir_bp : casc_2o_iir_base<m_t> {
public:
    casc_2o_iir_bp()
    {
        static_assert(m_t % 2 == 0, "M must be even!");
    }

    void copy_coeff_from(const casc_2o_iir_bp<m_t> &other_filter)
    {
        this->gain = other_filter.gain;
        this->aCoeff = other_filter.aCoeff;
    }

    template <typename iter_t>
    void process(iter_t begin, iter_t end)
    {
        this->template process_base<casc_2o_iir_bp<m_t>>(begin, end);
    }

    static void process_spec(std::array<std::array<double, 3>, m_t + 1> &y,
                             const std::array<std::array<double, 3>, m_t> &a, int p, int d1, int d2)
    {
        for (uint j = 0; j < m_t; j++) {
            y.at(j + 1).at(p) = y.at(j).at(p);

            y.at(j + 1).at(p) += -y.at(j + 1).at(d1) * a.at(j).at(1);
            y.at(j + 1).at(p) += -y.at(j).at(d2) - y.at(j + 1).at(d2) * a.at(j).at(2);
        }
    }

    void set_bp_coeff(double f0, double fs, double q, double gain_in = 1.0)
    {
        double q2{ 2 * q };
        this->m_gain = gain_in;

        double e0{ 2 * M_PI * f0 / fs };
        double dnm{ std::sin(e0) };

        double de{ 2 * std::tan(e0 / q2) / dnm };

        for (unsigned int k = 0; k < m_t / 2; k++) {
            double d{ 2 * std::sin((2 * k + 1) * M_PI / (2.0 * m_t)) };

            double a{ (1 + de * de / 4.0) * 2 / d / de };
            double dk{ std::sqrt(de * d / (a + std::sqrt(a * a - 1))) };

            double b{ d * de / dk / 2.0 };
            double w{ b + std::sqrt(b * b - 1) };

            double t{ std::tan(e0 / 2.0) };

            double e1{ 2.0 * std::atan(t / w) };
            double e2{ 2.0 * std::atan(w * t) };

            t = dk * std::sin(e1) / 2.0;
            dnm = (1 + t);
            double beta1{ (1 - t) / dnm / 2.0 };

            t = dk * std::sin(e2) / 2.0;
            dnm = (1 + t);
            double beta2{ (1 - t) / dnm / 2.0 };

            double gamma1{ (0.5 + beta1) * std::cos(e1) };
            double gamma2{ (0.5 + beta2) * std::cos(e2) };

            t = std::sqrt(1 + (w - 1 / w) / dk * (w - 1 / w) / dk);
            double alpha1{ (0.5 - beta1) * t / 2.0 };
            double alpha2{ (0.5 - beta2) * t / 2.0 };

            this->m_gain *= 4 * alpha1 * alpha2;
            // bCoeff.at(2 * k).at(0) = 1.0;
            // bCoeff.at(2 * k + 1).at(0) = 1.0;
            // bCoeff.at(2 * k).at(1) = 0;
            // bCoeff.at(2 * k + 1).at(1) = 0;
            // bCoeff.at(2 * k).at(2) = -1.0;
            // bCoeff.at(2 * k + 1).at(2) = -1.0;

            this->m_a_coeff.at(2 * k).at(0) = 1;
            this->m_a_coeff.at(2 * k + 1).at(0) = 1;
            this->m_a_coeff.at(2 * k).at(1) = -2 * gamma1;
            this->m_a_coeff.at(2 * k + 1).at(1) = -2 * gamma2;
            this->m_a_coeff.at(2 * k).at(2) = 2 * beta1;
            this->m_a_coeff.at(2 * k + 1).at(2) = 2 * beta2;
        }
    }
};
}
