/** \file   fatWaterMultiPeakR2StarOperator.h
    \brief  Implement fat water multi-peak signal model with the R2* signal decay
    \author Hui Xue
*/

#pragma once

#include "fatWaterMultiPeakOperator.h"

namespace Gadgetron {

    // y = ( rhoW+rhoF*sum(alpha_p*exp(j*2*pi*fp*te) )*exp(j*2*pi*fb*te)*exp(-R2star*te), for every fat peak p=1...P
    // b[0] : rhoW
    // b[1] : rhoF
    // b[2] : fB
    // b[3] : R2Star

    template <class ARRAY> class fatWaterMultiPeakR2StarOperator : public fatWaterMultiPeakOperator<ARRAY>
    {
    public:

        typedef fatWaterMultiPeakOperator<ARRAY> BaseClass;
        typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
        typedef typename BaseClass::REAL REAL;

        fatWaterMultiPeakR2StarOperator();
        virtual ~fatWaterMultiPeakR2StarOperator();

        virtual void gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad);

        // x : echo time
        // b : parameters
        virtual void magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y);

        using BaseClass::fat_freqs_;
        using BaseClass::fat_rel_amps_;

    protected:
    };

    template <class ARRAY>
    fatWaterMultiPeakR2StarOperator<ARRAY>::fatWaterMultiPeakR2StarOperator() : BaseClass()
    {
    }

    template <class ARRAY>
    fatWaterMultiPeakR2StarOperator<ARRAY>::~fatWaterMultiPeakR2StarOperator()
    {
    }

    template <class ARRAY>
    void fatWaterMultiPeakR2StarOperator<ARRAY>::gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad)
    {
        try
        {
            size_t num = b.size();
            if(grad.size()!=num) grad.resize(num, 0);

            size_t num_fat_peaks = fat_freqs_.size();

            ELEMENT_TYPE v(0, 2*M_PI*b[2].real()*xi.real());
            ELEMENT_TYPE exp_v = std::exp(v);

            ELEMENT_TYPE exp_v2 = std::exp(-b[3].real()*xi.real());

            // to rhoW
            grad[0] = exp_v * exp_v2;

            // to rhoF
            size_t pp;
            ELEMENT_TYPE fat_w(0);
            for (pp=0; pp<num_fat_peaks; pp++)
            {
                ELEMENT_TYPE v(0, 2*M_PI*fat_freqs_[pp]*xi.real());
                fat_w += fat_rel_amps_[pp]*std::exp(v);
            }

            grad[1] = fat_w * exp_v * exp_v2;

            // to fB
            ELEMENT_TYPE v3(0, 2*M_PI*xi.real());
            grad[2] = (b[0]+b[1]*fat_w) * exp_v * exp_v2 * v3;

            // to R2*
            grad[3] = (b[0]+b[1]*fat_w) * exp_v * exp_v2 * -xi;
        }
        catch(...)
        {
            GADGET_THROW("Errors happened in fatWaterMultiPeakR2StarOperator<ARRAY>::gradient(...) ... ");
        }
    }

    template <class ARRAY>
    void fatWaterMultiPeakR2StarOperator<ARRAY>::magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y)
    {
        size_t num = x.size();
        if(y.size()!=x.size()) y.resize(num, 0);

        BaseClass::magnitude(x, b, y);

        size_t ii;
        for (ii=0; ii<num; ii++)
        {
            y[ii] *= std::exp( -b[3].real()*x[ii].real() );
        }
    }
}
