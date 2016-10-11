
#pragma once

#include "sashahclib_export.h"
#include "GenericReconBase.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"

namespace Gadgetron {

    class EXPORTSASHAHC MoCoSashaHCGadget : public GenericReconImageBase
    {
    public:
        GADGET_DECLARE(MoCoSashaHCGadget);

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;

        typedef GenericReconImageBase BaseClass;

        MoCoSashaHCGadget();
        ~MoCoSashaHCGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the moco
        /// ------------------------------------------------------------------------------------
        GADGET_PROPERTY(send_ori, bool, "Whether to send original images", true);
        GADGET_PROPERTY(send_moco, bool, "Whether to send moco images", true);
        GADGET_PROPERTY(send_diff_image, bool, "Whether to send difference images", true);

        GADGET_PROPERTY(regularization_hilbert_strength, double, "Hilbert stregth for moco", 12.0);
        GADGET_PROPERTY(bidirectional_moco, bool, "Whether to apply bidirectional moco", true);

    protected:

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // iteration of moco
        std::vector<unsigned int> iters_;

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1);

        // close call
        int close(unsigned long flags);
    };
}
