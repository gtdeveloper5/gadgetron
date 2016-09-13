/** \file   fatWaterMultiPeakOperator
    \brief  Implement fat water multi-peak signal model
    \author Hui Xue
*/

#pragma once

#include "curveFittingOperator.h"
#include <complex>
#include <cmath>

#ifdef M_PI
    #undef M_PI
#endif // M_PI
#define M_PI 3.14159265358979323846

namespace Gadgetron {

    // y = ( rhoW+rhoF*sum(alpha_p*exp(j*2*pi*fp*te) )*exp(j*2*pi*fb*te), for every fat peak p=1...P
    // b[0] : rhoW
    // b[1] : rhoF
    // b[2] : fB

    template <class ARRAY> class fatWaterMultiPeakOperator : public curveFittingOperator<ARRAY>
    {
    public:

        typedef curveFittingOperator<ARRAY> BaseClass;
        typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
        typedef typename BaseClass::REAL REAL;

        fatWaterMultiPeakOperator();
        virtual ~fatWaterMultiPeakOperator();

        virtual void gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad);

        // x : echo time
        // b : parameters
        virtual void magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y);

        // fat frequencies
        std::vector<REAL> fat_freqs_;

        // relative amplitudes for fat
        std::vector<REAL> fat_rel_amps_;

    protected:
    };

    template <class ARRAY>
    fatWaterMultiPeakOperator<ARRAY>::fatWaterMultiPeakOperator() : BaseClass()
    {
        fat_freqs_.resize(6);
        fat_freqs_[0] = -3.80; // TODO: these may be ppm; need to convert to Hz
        fat_freqs_[1] = -3.40;
        fat_freqs_[2] = -2.60;
        fat_freqs_[3] = -1.94;
        fat_freqs_[4] = -0.39;
        fat_freqs_[5] = 0.60;

        fat_rel_amps_.resize(6);
        fat_rel_amps_[0] = 0.087;
        fat_rel_amps_[1] = 0.693;
        fat_rel_amps_[2] = 0.128;
        fat_rel_amps_[3] = 0.004;
        fat_rel_amps_[4] = 0.039;
        fat_rel_amps_[5] = 0.048;
    }

    template <class ARRAY>
    fatWaterMultiPeakOperator<ARRAY>::~fatWaterMultiPeakOperator()
    {
    }

    template <class ARRAY>
    void fatWaterMultiPeakOperator<ARRAY>::gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad)
    {
        try
        {
            size_t num = b.size();
            if(grad.size()!=num) grad.resize(num, 0);

            size_t num_fat_peaks = fat_freqs_.size();

            ELEMENT_TYPE v(0, 2*M_PI*b[2].real()*xi.real());
            ELEMENT_TYPE exp_v = std::exp(v);

            // to rhoW
            grad[0] = exp_v;

            // to rhoF
            size_t pp;
            ELEMENT_TYPE fat_w(0);
            for (pp=0; pp<num_fat_peaks; pp++)
            {
                ELEMENT_TYPE v(0, 2*M_PI*fat_freqs_[pp]*xi.real());
                fat_w += fat_rel_amps_[pp]*std::exp(v);
            }

            grad[1] = fat_w*exp_v;

            // to fB
            ELEMENT_TYPE v2(0, 2*M_PI*xi.real());
            grad[2] = (b[0]+b[1]*fat_w) * exp_v * v2;
        }
        catch(...)
        {
            GADGET_THROW("Errors happened in fatWaterMultiPeakOperator<ARRAY>::gradient(...) ... ");
        }
    }

    template <class ARRAY>
    void fatWaterMultiPeakOperator<ARRAY>::magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y)
    {
        size_t num = x.size();
        if(y.size()!=x.size()) y.resize(num, 0);

        size_t num_fat_peaks = fat_freqs_.size();

        size_t ii, pp;
        for (ii=0; ii<num; ii++)
        {
            ELEMENT_TYPE v(0, 2*M_PI*b[2].real()*x[ii].real());
            ELEMENT_TYPE exp_v = std::exp(v);

            ELEMENT_TYPE fat_w(0);
            for (pp=0; pp<num_fat_peaks; pp++)
            {
                ELEMENT_TYPE v(0, 2*M_PI*fat_freqs_[pp]*x[ii].real());
                fat_w += fat_rel_amps_[pp]*std::exp(v);
            }

            y[ii] = exp_v * (b[0] + b[1]*fat_w);
        }
    }
}
