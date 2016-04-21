
#include "GenericReconCartesianNonLinearSpirit2DTGadget.h"
#include "hoSPIRIT2DTOperator.h"

namespace Gadgetron {

    GenericReconCartesianNonLinearSpirit2DTGadget::GenericReconCartesianNonLinearSpirit2DTGadget() : BaseClass()
    {
    }

    GenericReconCartesianNonLinearSpirit2DTGadget::~GenericReconCartesianNonLinearSpirit2DTGadget()
    {
    }

    int GenericReconCartesianNonLinearSpirit2DTGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        // -------------------------------------------------

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        // -------------------------------------------------
        // check the parameters
        if(this->spirit_nl_iter_max.value()==0)
        {
            this->spirit_nl_iter_max.value(15);
            GDEBUG_STREAM("spirit_iter_max: " << this->spirit_nl_iter_max.value());
        }

        if (this->spirit_nl_iter_thres.value()<FLT_EPSILON)
        {
            this->spirit_nl_iter_thres.value(0.004);
            GDEBUG_STREAM("spirit_nl_iter_thres: " << this->spirit_nl_iter_thres.value());
        }

        if (this->spirit_image_reg_lamda.value() < FLT_EPSILON)
        {
            if(this->spirit_reg_proximity_across_cha.value())
            {
                this->spirit_image_reg_lamda.value(0.0002);
            }
            else
            {
                this->spirit_image_reg_lamda.value(0.00005);
            }

            GDEBUG_STREAM("spirit_image_reg_lamda: " << this->spirit_image_reg_lamda.value());
        }

        if (this->spirit_reg_N_weighting_ratio.value() < FLT_EPSILON)
        {
            if(acceFactorE1_[0]<=5)
            {
                this->spirit_reg_N_weighting_ratio.value(10.0);
            }
            else
            {
                this->spirit_reg_N_weighting_ratio.value(20.0);
            }

            GDEBUG_STREAM("spirit_reg_N_weighting_ratio: " << this->spirit_reg_N_weighting_ratio.value());
        }

        return GADGET_OK;
    }

    void GenericReconCartesianNonLinearSpirit2DTGadget::perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            typedef std::complex<float> T;

            size_t RO = recon_bit.data_.data_.get_size(0);
            size_t E1 = recon_bit.data_.data_.get_size(1);
            size_t E2 = recon_bit.data_.data_.get_size(2);
            size_t dstCHA = recon_bit.data_.data_.get_size(3);
            size_t N = recon_bit.data_.data_.get_size(4);
            size_t S = recon_bit.data_.data_.get_size(5);
            size_t SLC = recon_bit.data_.data_.get_size(6);

            hoNDArray< std::complex<float> >& src = recon_obj.ref_calib_;

            size_t ref_RO = src.get_size(0);
            size_t ref_E1 = src.get_size(1);
            size_t ref_E2 = src.get_size(2);
            size_t srcCHA = src.get_size(3);
            size_t ref_N = src.get_size(4);
            size_t ref_S = src.get_size(5);
            size_t ref_SLC = src.get_size(6);

            size_t convkRO = recon_obj.kernel_.get_size(0);
            size_t convkE1 = recon_obj.kernel_.get_size(1);
            size_t convkE2 = recon_obj.kernel_.get_size(2);

            recon_obj.recon_res_.data_.create(RO, E1, E2, 1, N, S, SLC);
            Gadgetron::clear(recon_obj.recon_res_.data_);
            recon_obj.full_kspace_ = recon_bit.data_.data_;
            Gadgetron::clear(recon_obj.full_kspace_);

            std::stringstream os;
            os << "encoding_" << e;
            std::string suffix = os.str();

            // if (!debug_folder_full_path_.empty()) { gt_exporter_.exportArrayComplex(recon_bit.data_.data_, debug_folder_full_path_ + "data_src_" + suffix); }

            // ------------------------------------------------------------------
            // compute effective acceleration factor
            // ------------------------------------------------------------------
            size_t e1, e2, n, s;
            size_t num_readout_lines = 0;
            for (s = 0; s < S; s++)
            {
                for (n = 0; n < N; n++)
                {
                    for (e2 = 0; e2 < E2; e2++)
                    {
                        for (e1 = 0; e1 < E1; e1++)
                        {
                            if (std::abs(recon_bit.data_.data_(RO / 2, e1, e2, 0, n)) > 0)
                            {
                                num_readout_lines++;
                            }
                        }
                    }
                }
            }

            if (num_readout_lines > 0)
            {
                double lenRO = RO;

                size_t start_RO = recon_bit.data_.sampling_.sampling_limits_[0].min_;
                size_t end_RO = recon_bit.data_.sampling_.sampling_limits_[0].max_;

                if ( (start_RO>=0 && start_RO<RO) && (end_RO>=0 && end_RO<RO) && (end_RO - start_RO + 1 < RO) )
                {
                    lenRO = (end_RO - start_RO + 1);
                }
                if (this->verbose.value()) GDEBUG_STREAM("length for RO : " << lenRO << " - " << lenRO / RO);

                double effectiveAcceFactor = (double)(S*N*E1*E2) / (num_readout_lines);
                if (this->verbose.value()) GDEBUG_STREAM("effectiveAcceFactor : " << effectiveAcceFactor);

                double ROScalingFactor = (double)RO / (double)lenRO;

                float fftCompensationRatio = (float)(std::sqrt(ROScalingFactor*effectiveAcceFactor));

                if (this->verbose.value()) GDEBUG_STREAM("fftCompensationRatio : " << fftCompensationRatio);

                Gadgetron::scal(fftCompensationRatio, recon_bit.data_.data_);
            }
            else
            {
                GWARN_STREAM("cannot find any sampled lines ... ");
            }

            Gadgetron::GadgetronTimer timer(false);

            // ------------------------------------------------------------------
            // compute the reconstruction
            // ------------------------------------------------------------------
            if(this->acceFactorE1_[e]<=1 && this->acceFactorE2_[e]<=1)
            {
                recon_obj.full_kspace_ = recon_bit.data_.data_;
            }
            else
            {
                hoNDArray< std::complex<float> >& kspace = recon_bit.data_.data_;
                hoNDArray< std::complex<float> >& res = recon_obj.full_kspace_;

                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_parallel_imaging_lamda : " << this->spirit_parallel_imaging_lamda.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_image_reg_lamda : " << this->spirit_image_reg_lamda.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_data_fidelity_lamda : " << this->spirit_data_fidelity_lamda.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_nl_iter_max : " << this->spirit_nl_iter_max.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_nl_iter_thres : " << this->spirit_nl_iter_thres.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_name : " << this->spirit_reg_name.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_level : " << this->spirit_reg_level.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_keep_approx_coeff : " << this->spirit_reg_keep_approx_coeff.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_keep_redundant_dimension_coeff : " << this->spirit_reg_keep_redundant_dimension_coeff.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_proximity_across_cha : " << this->spirit_reg_proximity_across_cha.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_use_coil_sen_map : " << this->spirit_reg_use_coil_sen_map.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_RO_weighting_ratio : " << this->spirit_reg_RO_weighting_ratio.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_E1_weighting_ratio : " << this->spirit_reg_E1_weighting_ratio.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_N_weighting_ratio : " << this->spirit_reg_N_weighting_ratio.value());

                size_t slc, s;

                for (slc = 0; slc < SLC; slc++)
                {
                    for (s = 0; s < S; s++)
                    {


                        if (this->perform_timing.value()) timer.start("SPIRIT 2D, linear unwrapping ... ");
                        this->perform_spirit_unwrapping(kspace, recon_obj.kernelIm2D_, res);
                        if (this->perform_timing.value()) timer.stop();

                        // if (!debug_folder_full_path_.empty()) { gt_exporter_.exportArrayComplex(res, debug_folder_full_path_ + "res_spirit_2D_" + suffix); }
                    }
                }
            }

            // ---------------------------------------------------------------------
            // compute coil combined images
            // ---------------------------------------------------------------------
            if (this->perform_timing.value()) timer.start("SPIRIT linear, coil combination ... ");

            if (E2>1)
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_obj.full_kspace_, complex_im_recon_buf_);
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_obj.full_kspace_, complex_im_recon_buf_);
            }

            // if (!debug_folder_full_path_.empty()) { gt_exporter_.exportArrayComplex(complex_im_recon_buf_, debug_folder_full_path_ + "complex_im_recon_buf_" + suffix); }

            size_t num = N*S*SLC;
            long long ii;

#pragma omp parallel default(none) private(ii) shared(num, N, S, recon_obj, RO, E1, E2, dstCHA) if(num>1)
            {
                hoNDArray< std::complex<float> > complexImBuf(RO, E1, E2, dstCHA);

#pragma omp for 
                for (ii = 0; ii < num; ii++)
                {
                    size_t slc = ii / (N*S);
                    size_t s = (ii - slc*N*S) / N;
                    size_t n = ii - slc*N*S - s*N;

                    size_t coilMapN = n;
                    if (coilMapN >= recon_obj.coil_map_.get_size(5)) coilMapN = recon_obj.coil_map_.get_size(5) - 1;

                    size_t coilMapS = s;
                    if (coilMapS >= recon_obj.coil_map_.get_size(6)) coilMapS = recon_obj.coil_map_.get_size(6) - 1;

                    hoNDArray< std::complex<float> > complexIm(RO, E1, E2, dstCHA, &(complex_im_recon_buf_(0, 0, 0, 0, n, s, slc)));
                    hoNDArray< std::complex<float> > coilMap(RO, E1, E2, dstCHA, &(recon_obj.coil_map_(0, 0, 0, 0, coilMapN, coilMapS, slc)));
                    hoNDArray< std::complex<float> > combined(RO, E1, E2, 1, &(recon_obj.recon_res_.data_(0, 0, 0, 0, n, s, slc)));

                    Gadgetron::multiplyConj(complexIm, coilMap, complexImBuf);
                    Gadgetron::sum_over_dimension(complexImBuf, combined, 3);
                }
            }

            if (this->perform_timing.value()) timer.stop();

            // if (!debug_folder_full_path_.empty()) { gt_exporter_.exportArrayComplex(recon_obj.recon_res_.data_, debug_folder_full_path_ + "unwrappedIm_" + suffix); }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconCartesianNonLinearSpirit2DTGadget::perform_unwrapping(...) ... ");
        }
    }

    void GenericReconCartesianNonLinearSpirit2DTGadget::perform_nonlinear_spirit_unwrapping(hoNDArray< std::complex<float> >& kspace, hoNDArray< std::complex<float> >& kerIm, hoNDArray< std::complex<float> >& res)
    {
        try
        {
            bool print_iter = this->spirit_print_iter.value();

            size_t RO = kspace.get_size(0);
            size_t E1 = kspace.get_size(1);
            size_t E2 = kspace.get_size(2);
            size_t CHA = kspace.get_size(3);
            size_t N = kspace.get_size(4);
            size_t S = kspace.get_size(5);
            size_t SLC = kspace.get_size(6);

            size_t ref_N = kerIm.get_size(4);
            size_t ref_S = kerIm.get_size(5);

        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconCartesianNonLinearSpirit2DTGadget::perform_nonlinear_spirit_unwrapping(...) ... ");
        }
    }

    GADGET_FACTORY_DECLARE(GenericReconCartesianNonLinearSpirit2DTGadget)
}
