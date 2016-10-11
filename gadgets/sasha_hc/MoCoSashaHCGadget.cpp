
#include "MoCoSashaHCGadget.h"
#include <iomanip>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "cmr_motion_correction.h"

namespace Gadgetron {

    MoCoSashaHCGadget::MoCoSashaHCGadget() : BaseClass()
    {
    }

    MoCoSashaHCGadget::~MoCoSashaHCGadget()
    {
    }

    int MoCoSashaHCGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(),h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        num_encoding_spaces_ = h.encoding.size();
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << num_encoding_spaces_);

        iters_.resize(4);
        iters_[0] = 32;
        iters_[1] = 64;
        iters_[2] = 100;
        iters_[3] = 100;

        return GADGET_OK;
    }

    int MoCoSashaHCGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1)
    {
        if (perform_timing.value()) { gt_timer_.start("MoCoSashaHCGadget::process"); }

        GDEBUG_CONDITION_STREAM(verbose.value(), "MoCoSashaHCGadget::process(...) starts ... ");

        process_called_times_++;

        std::stringstream os;
        os << "called_" << process_called_times_;
        std::string str = os.str();

        IsmrmrdImageArray* recon_res_ = m1->getObjectPtr();

        size_t encoding = (size_t)recon_res_->meta_[0].as_long("encoding", 0);
        GADGET_CHECK_RETURN(encoding<num_encoding_spaces_, GADGET_FAIL);

        std::string imageInfo;

        size_t RO = recon_res_->data_.get_size(0);
        size_t E1 = recon_res_->data_.get_size(1);
        size_t E2 = recon_res_->data_.get_size(2);
        size_t CHA = recon_res_->data_.get_size(3);
        size_t N = recon_res_->data_.get_size(4);
        size_t S = recon_res_->data_.get_size(5);
        size_t SLC = recon_res_->data_.get_size(6);

        GADGET_CHECK_RETURN(S == 2, GADGET_FAIL);

        size_t n, s, slc;

        // ----------------------------------------------------
        // compute diff images
        // ----------------------------------------------------

        hoNDArray<T> diff;
        diff.create(RO, E1, E2, CHA, N, 1, SLC);

        for (slc=0; slc<SLC; slc++)
        {
            T* pFirstN = &(recon_res_->data_(0, 0, 0, 0, 0, 0, slc));
            T* pSecondN = &(recon_res_->data_(0, 0, 0, 0, 0, 1, slc));

            hoNDArray<T> firstArray(RO, E1, E2, CHA, N, pFirstN);
            hoNDArray<T> secondArray(RO, E1, E2, CHA, N, pSecondN);

            T* pDiffN = &(diff(0, 0, 0, 0, 0, 0, slc));
            hoNDArray<T> diffArray(RO, E1, E2, CHA, N, pDiffN);

            Gadgetron::subtract(firstArray, secondArray, diffArray);
        }

        if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.export_array_complex(diff, debug_folder_full_path_ + "MoCoSashaHC_diff_" + str);
        }

        // find the key frame and run the moco
        hoNDArray<real_value_type> diffMag;
        Gadgetron::abs(diff, diffMag);

        if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.export_array(diffMag, debug_folder_full_path_ + "MoCoSashaHC_diffMag_" + str);
        }

        bool warp_input = false;

        hoNDArray<T> moco;
        moco.create(RO, E1, E2, CHA, N, 1, SLC);

        for (slc=0; slc<SLC; slc++)
        {
            std::stringstream os;
            os << "slice_" << slc;
            std::string str_slc = os.str();

            real_value_type* pMag = &diffMag(0, 0, 0, 0, 0, 0, slc);

            hoNDArray<real_value_type> target(RO, E1, N, pMag);
            size_t key_frame(0);
            Gadgetron::find_key_frame_SSD_2DT(target, key_frame);

            // perform moco
            std::string moco_str = "GenericReconCartesianSpiritGadget, perform moco for " + str_slc;
            if (perform_timing.value()) { gt_timer_.start(moco_str.c_str()); }
            Gadgetron::hoImageRegContainer2DRegistration<real_value_type, float, 2, 2> reg;
            Gadgetron::perform_moco_fixed_key_frame_2DT(target, key_frame, (real_value_type)(regularization_hilbert_strength.value()), iters_, bidirectional_moco.value(), warp_input, reg);
            if (perform_timing.value()) { gt_timer_.stop(); }

            // apply the deformation

            moco_str = "GenericReconCartesianSpiritGadget, apply deformation field for " + str_slc;

            T* pFirstN = &(recon_res_->data_(0, 0, 0, 0, 0, 0, slc));
            hoNDArray<T> firstArray(RO, E1, E2, CHA, N, pFirstN);

            hoNDArray<real_value_type> firstArray_real, firstArray_imag;
            Gadgetron::complex_to_real_imag(firstArray, firstArray_real, firstArray_imag);

            firstArray_real.squeeze();
            firstArray_imag.squeeze();

            hoNDArray<float> dx, dy;
            reg.deformation_field_[0].to_NDArray(0, dx);
            reg.deformation_field_[1].to_NDArray(0, dy);

            hoNDArray<real_value_type> firstArray_real_moco, firstArray_imag_moco;
            Gadgetron::apply_deformation_field(firstArray_real, dx, dy, firstArray_real_moco);
            Gadgetron::apply_deformation_field(firstArray_imag, dx, dy, firstArray_imag_moco);

            pFirstN = &(moco(0, 0, 0, 0, 0, 0, slc));
            hoNDArray<T> firstArrayMOCO(RO, E1, N, pFirstN);
            Gadgetron::real_imag_to_complex(firstArray_real_moco, firstArray_imag_moco, firstArrayMOCO);

            if (perform_timing.value()) { gt_timer_.stop(); }

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "called_" << process_called_times_ << "_slc" << slc;
                std::string str = os.str();

                gt_exporter_.export_array_complex(firstArrayMOCO, debug_folder_full_path_ + "MoCoSashaHC_diffMag_MOCO_" + str);

                gt_exporter_.export_array(dx, debug_folder_full_path_ + "MoCoSashaHC_dx_" + str);
                gt_exporter_.export_array(dy, debug_folder_full_path_ + "MoCoSashaHC_dy_" + str);
            }
        }

        GDEBUG_CONDITION_STREAM(verbose.value(), "MoCoSashaHCGadget::process(...) ends ... ");

        if(send_ori.value())
        {
            Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1_sasha = new Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >();
            m1_sasha->getObjectPtr()->data_.create(RO, E1, E2, CHA, N, 1, SLC);
            m1_sasha->getObjectPtr()->headers_.create(N, 1, SLC);
            m1_sasha->getObjectPtr()->meta_.resize(N*1*SLC);

            Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1_sasha_hc = new Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >();
            m1_sasha_hc->getObjectPtr()->data_.create(RO, E1, E2, CHA, N, 1, SLC);
            m1_sasha_hc->getObjectPtr()->headers_.create(N, 1, SLC);
            m1_sasha_hc->getObjectPtr()->meta_.resize(N*1*SLC);

            Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1_sasha_diff = new Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >();
            m1_sasha_diff->getObjectPtr()->data_.create(RO, E1, E2, CHA, N, 1, SLC);
            m1_sasha_diff->getObjectPtr()->headers_.create(N, 1, SLC);
            m1_sasha_diff->getObjectPtr()->meta_.resize(N*1*SLC);

            for (slc=0; slc<SLC; slc++)
            {
                for (n=0; n<N; n++)
                {
                    memcpy( &(m1_sasha->getObjectPtr()->data_(0, 0, 0, 0, n, 0, slc)), 
                        &(recon_res_->data_(0, 0, 0, 0, n, 0, slc)), sizeof(T)*RO*E1*E2*CHA);

                    memcpy( &(m1_sasha_hc->getObjectPtr()->data_(0, 0, 0, 0, n, 0, slc)), 
                        &(recon_res_->data_(0, 0, 0, 0, n, 1, slc)), sizeof(T)*RO*E1*E2*CHA);

                    memcpy( &(m1_sasha_diff->getObjectPtr()->data_(0, 0, 0, 0, n, 0, slc)), 
                        &(diff(0, 0, 0, 0, n, 0, slc)), sizeof(T)*RO*E1*E2*CHA);

                    // ------------------------------

                    m1_sasha->getObjectPtr()->headers_(n, 0, slc) = recon_res_->headers_(n, 0, slc);
                    m1_sasha->getObjectPtr()->headers_(n, 0, slc).image_series_index += 1;

                    m1_sasha->getObjectPtr()->meta_[n+slc*N] = recon_res_->meta_[n+0*N+slc*N*S];
                    m1_sasha->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE, "SASHA");

                    m1_sasha->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, "SASHA");
                    m1_sasha->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "SASHA");
                    // ------------------------------

                    m1_sasha_hc->getObjectPtr()->headers_(n, 0, slc) = recon_res_->headers_(n, 1, slc);
                    m1_sasha_hc->getObjectPtr()->headers_(n, 0, slc).image_series_index += 2;

                    m1_sasha_hc->getObjectPtr()->meta_[n+slc*N] = recon_res_->meta_[n+1*N+slc*N*S];
                    m1_sasha_hc->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE, "SASHA_HC");

                    m1_sasha_hc->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, "SASHA_HC");
                    m1_sasha_hc->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "SASHA_HC");
                    // ------------------------------

                    m1_sasha_diff->getObjectPtr()->headers_(n, 0, slc) = recon_res_->headers_(n, 0, slc);
                    m1_sasha_diff->getObjectPtr()->headers_(n, 0, slc).image_series_index += 3;

                    m1_sasha_diff->getObjectPtr()->meta_[n+slc*N] = recon_res_->meta_[n+0*N+slc*N*S];
                    m1_sasha_diff->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE, "SASHA_HC_DIFF");

                    m1_sasha_diff->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, "SASHA_HC_DIFF");
                    m1_sasha_diff->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "SASHA_HC_DIFF");
                }
            }

            if (this->next()->putq(m1_sasha) == -1)
            {
                GERROR("MoCoSashaHCGadget::process, passing m1_sasha on to next gadget");
                return GADGET_FAIL;
            }

            if (this->next()->putq(m1_sasha_hc) == -1)
            {
                GERROR("MoCoSashaHCGadget::process, passing m1_sasha_hc on to next gadget");
                return GADGET_FAIL;
            }

            if(send_diff_image.value())
            {
                if (this->next()->putq(m1_sasha_diff) == -1)
                {
                    GERROR("MoCoSashaHCGadget::process, passing m1_sasha_diff on to next gadget");
                    return GADGET_FAIL;
                }
            }
        }

        if(send_moco.value())
        {
            Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1_sasha_moco = new Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >();
            m1_sasha_moco->getObjectPtr()->data_.create(RO, E1, E2, CHA, N, 1, SLC);
            m1_sasha_moco->getObjectPtr()->headers_.create(N, 1, SLC);
            m1_sasha_moco->getObjectPtr()->meta_.resize(N*1*SLC);

            for (slc=0; slc<SLC; slc++)
            {
                for (n=0; n<N; n++)
                {
                    memcpy( &(m1_sasha_moco->getObjectPtr()->data_(0, 0, 0, 0, n, 0, slc)), 
                        &(moco(0, 0, 0, 0, n, 0, slc)), sizeof(T)*RO*E1*E2*CHA);

                    // ------------------------------

                    m1_sasha_moco->getObjectPtr()->headers_(n, 0, slc) = recon_res_->headers_(n, 0, slc);
                    m1_sasha_moco->getObjectPtr()->headers_(n, 0, slc).image_series_index += 4;

                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N] = recon_res_->meta_[n+0*N+slc*N*S];
                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE, "SASHA");
                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_DATA_ROLE, "MOCO");

                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, "SASHA");
                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGECOMMENT, "MOCO");

                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_IMAGEPROCESSINGHISTORY, "MOCO");

                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "SASHA");
                    m1_sasha_moco->getObjectPtr()->meta_[n+slc*N].append(GADGETRON_SEQUENCEDESCRIPTION, "MOCO");
                }
            }

            if (this->next()->putq(m1_sasha_moco) == -1)
            {
                GERROR("MoCoSashaHCGadget::process, passing m1_sasha_moco on to next gadget");
                return GADGET_FAIL;
            }
        }

        if (perform_timing.value()) { gt_timer_.stop(); }

        return GADGET_OK;
    }

    int MoCoSashaHCGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "MoCoSashaHCGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        if (flags != 0)
        {
        }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(MoCoSashaHCGadget)

}
