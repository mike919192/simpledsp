#pragma once
#include <array>
#include <cmath>
#include "FilterType.h"

namespace sdsp
{
    template <size_t M>
    class casc2orderIIR 
    {
    private:
        
        int pos {0};
        
        std::array<std::array<double, 3>, M + 1> mem {0};

	    std::array<std::array<double, 3>, M> bCoeff {0};
	    std::array<std::array<double, 3>, M> aCoeff {0};

        FilterType fType {FilterType::None};

    public:

        casc2orderIIR()
        {
            static_assert(M % 2 == 0, "M must be even!");
        }

        template <typename Iter>
        void Process(Iter begin, Iter end) 
        {
            constexpr int order {2};
            constexpr int j1 {1};
            constexpr int j2 {2};

            int p {pos};

            std::array<std::array<double, 3>, M + 1> y = mem;

            std::array<std::array<double, 3>, M> b = bCoeff;

            std::array<std::array<double, 3>, M> a = aCoeff;

            while (begin < end) {
                y.at(0).at(p) = *begin;

                int d1 {p - j1};
                if (d1 < 0)
                    d1 += order + 1;

                int d2 {p - j2};
                if (d2 < 0)
                    d2 += order + 1;

                uint j;

                for (j = 0; j < M; j++) {
                    y.at(j + 1).at(p) = y.at(j).at(p) * b.at(j).at(0);

                    y.at(j + 1).at(p) += y.at(j).at(d1) * b.at(j).at(j1) - y.at(j + 1).at(d1) * a.at(j).at(j1);
                    y.at(j + 1).at(p) += y.at(j).at(d2) * b.at(j).at(j2) - y.at(j + 1).at(d2) * a.at(j).at(j2);
                }

                *begin = y.at(j).at(p);

                p++;
                if (p > order)
                    p = 0;
                begin++;
            }
            pos = p;
            mem = y;
        }
        
        void SetBPCoeff(double f0, double fs, double Q) 
        {
            double q2 {2 * Q};
            fType = FilterType::BandPass;

            double e0 {2 * M_PI * f0 / fs};
            double dnm {std::sin(e0)};
            
            double de {2 * std::tan(e0 / q2) / dnm};
            
            for (unsigned int k = 0; k < M / 2; k++) {
                double D {2 * std::sin((2 * k + 1) * M_PI / (2.0 * M))};
                
                double A {(1 + de * de / 4.0) * 2 / D / de};
                double dk {std::sqrt(de * D / (A + std::sqrt(A * A - 1)))};
                
                double B {D * de / dk / 2.0};
                double W {B + std::sqrt(B * B - 1)};
                
                double t {std::tan(e0 / 2.0)};
                
                double e1 {2.0 * std::atan(t / W)};
                double e2 {2.0 * std::atan(W * t)};

                t = dk * std::sin(e1) / 2.0;
                dnm = (1 + t);
                double beta1 {(1 - t) / dnm / 2.0};

                t = dk * std::sin(e2) / 2.0;
                dnm = (1 + t);
                double beta2 {(1 - t) / dnm / 2.0};

                double gamma1 {(0.5 + beta1) * std::cos(e1)};
                double gamma2 {(0.5 + beta2) * std::cos(e2)};

                t = std::sqrt(1 + (W - 1 / W) / dk * (W - 1 / W) / dk);
                double alpha1 {(0.5 - beta1) * t / 2.0};
                double alpha2 {(0.5 - beta2) * t / 2.0};

                bCoeff.at(2 * k).at(0) = 2 * alpha1;
                bCoeff.at(2 * k + 1).at(0) = 2 * alpha2;
                bCoeff.at(2 * k).at(1) = 0;
                bCoeff.at(2 * k + 1).at(1) = 0;
                bCoeff.at(2 * k).at(2) = -bCoeff.at(2 * k).at(0);
                bCoeff.at(2 * k + 1).at(2) = -bCoeff.at(2 * k + 1).at(0);

                aCoeff.at(2 * k).at(0) = 1;
                aCoeff.at(2 * k + 1).at(0) = 1;
                aCoeff.at(2 * k).at(1) = -2 * gamma1;
                aCoeff.at(2 * k + 1).at(1) = -2 * gamma2;
                aCoeff.at(2 * k).at(2) = 2 * beta1;
                aCoeff.at(2 * k + 1).at(2) = 2 * beta2;
            }
        }
        
        void SetHPCoeff(double f0, double fs) 
        {
            fType = FilterType::HighPass;
            
            double e0 {2 * M_PI * f0 / fs};
            for (unsigned int k = 0; k < M; k++) {
                double dk {2 * std::sin((2 * k + 1) * M_PI / (4.0 * M))};

                double t {dk * std::sin(e0) / 2};
                double dnm {1 + t};
                
                double beta1 {(1 - t) / dnm / 2};
                double gamma1 {(0.5 + beta1) * std::cos(e0)};
                double alpha1 {(0.5 + beta1 + gamma1) / 4}; //-

                bCoeff.at(k).at(0) = 2 * alpha1;
                bCoeff.at(k).at(1) = -2 * bCoeff.at(k).at(0); //-
                bCoeff.at(k).at(2) = bCoeff.at(k).at(0);

                aCoeff.at(k).at(0) = 1;
                aCoeff.at(k).at(1) = -2 * gamma1;
                aCoeff.at(k).at(2) = 2 * beta1;
            }
        }
        
        void SetLPCoeff(double f0, double fs) 
        {
            fType = FilterType::LowPass;

            double e0 {2 * M_PI * f0 / fs};
            for (unsigned int k = 0; k < M; k++) {
                double dk {2 * std::sin((2 * k + 1) * M_PI / (4.0 * M))};

                double t {dk * std::sin(e0) / 2};
                double dnm {1 + t};
                
                double beta1 {(1 - t) / dnm / 2};
                double gamma1 {(0.5 + beta1) * std::cos(e0)};
                double alpha1 {(0.5 + beta1 - gamma1) / 4};

                bCoeff.at(k).at(0) = 2 * alpha1;
                bCoeff.at(k).at(1) = 2 * bCoeff.at(k).at(0);
                bCoeff.at(k).at(2) = bCoeff.at(k).at(0);

                aCoeff.at(k).at(0) = 1;
                aCoeff.at(k).at(1) = -2 * gamma1;
                aCoeff.at(k).at(2) = 2 * beta1;
            }
        }
        
        //preload the filter memory for steady state input equal to value parameter
        void PreloadFilter(double value)
        {
            std::array<std::array<double, 3>, M + 1> memVals {0};
            for (int i = 0; i < 3; i++) {
                memVals.at(0).at(i) = value;
            }
            if (fType == FilterType::LowPass) {
                for (uint j = 1; j < M + 1; j++) {
                    for (uint i = 0; i < 3; i++) {
                        memVals.at(j).at(i) = value;
                    }
                }
            }
            mem = memVals;
        }
    };

    template<size_t M>
    class casc_2o_IIR
    {
        protected:
        int pos {0};

        double gain {1.0};
        
        std::array<std::array<double, 3>, M + 1> mem {0};
        
	    std::array<std::array<double, 3>, M> aCoeff {0};
    };

    template<size_t M>
    class casc_2o_IIR_lp : casc_2o_IIR<M>
    {
    public:

        casc_2o_IIR_lp()
        {
            static_assert(M % 2 == 0, "M must be even!");
        }

        template <typename Iter>
        void process(Iter begin, Iter end) 
        {
            constexpr int order {2};
            constexpr int j1 {1};
            constexpr int j2 {2};

            int p {this->pos};

            std::array<std::array<double, 3>, M + 1> y = this->mem;

            std::array<std::array<double, 3>, M> a = this->aCoeff;

            while (begin < end) {
                y.at(0).at(p) = *begin * this->gain;

                int d1 {p - j1};
                if (d1 < 0)
                    d1 += order + 1;

                int d2 {p - j2};
                if (d2 < 0)
                    d2 += order + 1;

                uint j {0};

                for (j = 0; j < M; j++) {
                    y.at(j + 1).at(p) = y.at(j).at(p);

                    y.at(j + 1).at(p) += y.at(j).at(d1) + y.at(j).at(d1) - y.at(j + 1).at(d1) * a.at(j).at(j1);
                    y.at(j + 1).at(p) += y.at(j).at(d2) - y.at(j + 1).at(d2) * a.at(j).at(j2);
                }

                *begin = y.at(j).at(p);

                p++;
                if (p > order)
                    p = 0;
                begin++;
            }
            this->pos = p;
            this->mem = y;
        }

        void set_coeff(double f0, double fs, double gainIn = 1.0) 
        {
            this->gain = gainIn;

            double e0 {2 * M_PI * f0 / fs};
            for (unsigned int k = 0; k < M; k++) {
                double dk {2 * std::sin((2 * k + 1) * M_PI / (4.0 * M))};

                double t {dk * std::sin(e0) / 2};
                double dnm {1 + t};
                
                double beta1 {(1 - t) / dnm / 2};
                double gamma1 {(0.5 + beta1) * std::cos(e0)};
                double alpha1 {(0.5 + beta1 - gamma1) / 4};

                this->gain *= 2 * alpha1;
                //bCoeff.at(k).at(0) = 1.0;
                //bCoeff.at(k).at(1) = 2.0;
                //bCoeff.at(k).at(2) = 1.0;

                this->aCoeff.at(k).at(0) = 1;
                this->aCoeff.at(k).at(1) = -2 * gamma1;
                this->aCoeff.at(k).at(2) = 2 * beta1;
            }
        }
    };

    template<size_t M>
    class casc_2o_IIR_hp : casc_2o_IIR<M>
    {
    public:

        casc_2o_IIR_hp()
        {
            static_assert(M % 2 == 0, "M must be even!");
        }

        template <typename Iter>
        void process(Iter begin, Iter end) 
        {
            constexpr int order {2};
            constexpr int j1 {1};
            constexpr int j2 {2};

            int p {this->pos};

            std::array<std::array<double, 3>, M + 1> y = this->mem;

            std::array<std::array<double, 3>, M> a = this->aCoeff;

            while (begin < end) {
                y.at(0).at(p) = *begin * this->gain;

                int d1 {p - j1};
                if (d1 < 0)
                    d1 += order + 1;

                int d2 {p - j2};
                if (d2 < 0)
                    d2 += order + 1;

                uint j {0};

                for (j = 0; j < M; j++) {
                    y.at(j + 1).at(p) = y.at(j).at(p);

                    y.at(j + 1).at(p) += - y.at(j).at(d1) - y.at(j).at(d1) - y.at(j + 1).at(d1) * a.at(j).at(j1);
                    y.at(j + 1).at(p) += y.at(j).at(d2) - y.at(j + 1).at(d2) * a.at(j).at(j2);
                }

                *begin = y.at(j).at(p);

                p++;
                if (p > order)
                    p = 0;
                begin++;
            }
            this->pos = p;
            this->mem = y;
        }

        void set_coeff(double f0, double fs, double gainIn = 1.0) 
        {
            this->gain = gainIn;
            
            double e0 {2 * M_PI * f0 / fs};
            for (unsigned int k = 0; k < M; k++) {
                double dk {2 * std::sin((2 * k + 1) * M_PI / (4.0 * M))};

                double t {dk * std::sin(e0) / 2};
                double dnm {1 + t};
                
                double beta1 {(1 - t) / dnm / 2};
                double gamma1 {(0.5 + beta1) * std::cos(e0)};
                double alpha1 {(0.5 + beta1 + gamma1) / 4};

                this->gain *= 2 * alpha1;
                //bCoeff.at(k).at(0) = 1.0;
                //bCoeff.at(k).at(1) = -2.0;
                //bCoeff.at(k).at(2) = 1.0;

                this->aCoeff.at(k).at(0) = 1;
                this->aCoeff.at(k).at(1) = -2 * gamma1;
                this->aCoeff.at(k).at(2) = 2 * beta1;
            }
        }
    };

    template<size_t M>
    class casc_2o_IIR_bp : casc_2o_IIR<M>
    {
    public:

        casc_2o_IIR_bp()
        {
            static_assert(M % 2 == 0, "M must be even!");
        }

        template <typename Iter>
        void process(Iter begin, Iter end) 
        {
            constexpr int order {2};
            constexpr int j1 {1};
            constexpr int j2 {2};

            int p {this->pos};

            std::array<std::array<double, 3>, M + 1> y = this->mem;

            std::array<std::array<double, 3>, M> a = this->aCoeff;

            while (begin < end) {
                y.at(0).at(p) = *begin * this->gain;

                int d1 {p - j1};
                if (d1 < 0)
                    d1 += order + 1;

                int d2 {p - j2};
                if (d2 < 0)
                    d2 += order + 1;

                uint j {0};

                for (j = 0; j < M; j++) {
                    y.at(j + 1).at(p) = y.at(j).at(p);

                    y.at(j + 1).at(p) += - y.at(j + 1).at(d1) * a.at(j).at(j1);
                    y.at(j + 1).at(p) += - y.at(j).at(d2) - y.at(j + 1).at(d2) * a.at(j).at(j2);
                }

                *begin = y.at(j).at(p);

                p++;
                if (p > order)
                    p = 0;
                begin++;
            }
            this->pos = p;
            this->mem = y;
        }

        void set_coeff(double f0, double fs, double Q, double gainIn = 1.0) 
        {
            double q2 {2 * Q};
            this->gain = gainIn;

            double e0 {2 * M_PI * f0 / fs};
            double dnm {std::sin(e0)};
            
            double de {2 * std::tan(e0 / q2) / dnm};
            
            for (unsigned int k = 0; k < M / 2; k++) {
                double D {2 * std::sin((2 * k + 1) * M_PI / (2.0 * M))};
                
                double A {(1 + de * de / 4.0) * 2 / D / de};
                double dk {std::sqrt(de * D / (A + std::sqrt(A * A - 1)))};
                
                double B {D * de / dk / 2.0};
                double W {B + std::sqrt(B * B - 1)};
                
                double t {std::tan(e0 / 2.0)};
                
                double e1 {2.0 * std::atan(t / W)};
                double e2 {2.0 * std::atan(W * t)};

                t = dk * std::sin(e1) / 2.0;
                dnm = (1 + t);
                double beta1 {(1 - t) / dnm / 2.0};

                t = dk * std::sin(e2) / 2.0;
                dnm = (1 + t);
                double beta2 {(1 - t) / dnm / 2.0};

                double gamma1 {(0.5 + beta1) * std::cos(e1)};
                double gamma2 {(0.5 + beta2) * std::cos(e2)};

                t = std::sqrt(1 + (W - 1 / W) / dk * (W - 1 / W) / dk);
                double alpha1 {(0.5 - beta1) * t / 2.0};
                double alpha2 {(0.5 - beta2) * t / 2.0};

                this->gain *= 4 * alpha1 * alpha2;
                //bCoeff.at(2 * k).at(0) = 1.0;
                //bCoeff.at(2 * k + 1).at(0) = 1.0;
                //bCoeff.at(2 * k).at(1) = 0;
                //bCoeff.at(2 * k + 1).at(1) = 0;
                //bCoeff.at(2 * k).at(2) = -1.0;
                //bCoeff.at(2 * k + 1).at(2) = -1.0;

                this->aCoeff.at(2 * k).at(0) = 1;
                this->aCoeff.at(2 * k + 1).at(0) = 1;
                this->aCoeff.at(2 * k).at(1) = -2 * gamma1;
                this->aCoeff.at(2 * k + 1).at(1) = -2 * gamma2;
                this->aCoeff.at(2 * k).at(2) = 2 * beta1;
                this->aCoeff.at(2 * k + 1).at(2) = 2 * beta2;
            }
        }
    };
}