//FlagSashaHCGadget.cpp

#include "FlagSashaHCGadget.h"
#include <iomanip>
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"

using namespace Gadgetron;

int FlagSashaHCGadget::process_config( ACE_Message_Block* mb )
{
    return GADGET_OK;
}

int FlagSashaHCGadget::process(
    GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1,
    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 )
{

    bool is_noise = ISMRMRD::FlagBit( ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT ).isSet( m1->getObjectPtr()->flags );

    if (!is_noise)
    {

        // 2016-09-28 Kelvin: Mark high-contrast images with the contrasts dimension
        ISMRMRD::AcquisitionHeader & acqhdr = *m1->getObjectPtr();
        static int16_t iLastLine = -1;
        static int16_t iLastSet = -1;

        // Reset iLastLine logic when a new image starts
        // Line reset to 0 is for ACS lines... FIXME!
        if (((int16_t)acqhdr.idx.set != iLastSet) || ((int16_t)acqhdr.idx.kspace_encode_step_1 == 0))
        {
            iLastSet = (int16_t)acqhdr.idx.set;
            iLastLine = -1;
        }

        // Assume that lines are always acquired in ascending order unless they are HC
        if ((int16_t)acqhdr.idx.kspace_encode_step_1 < iLastLine)
        {
            acqhdr.idx.contrast = 1;
        }
        else {
            iLastLine = (int16_t)acqhdr.idx.kspace_encode_step_1;
        }

        if (acqhdr.idx.contrast)
        {
            GDEBUG_STREAM( "*e1: " << std::setw( 3 ) << (int16_t)acqhdr.idx.kspace_encode_step_1
                << " set: " << (int16_t)acqhdr.idx.set
                << " con: " << (int16_t)acqhdr.idx.contrast
                << " iLastLine: " << iLastLine
                << " nav " << acqhdr.isFlagSet( ISMRMRD::ISMRMRD_ACQ_IS_RTFEEDBACK_DATA ) );
            //		                      << " user[1]: " << (int16_t)acqhdr.idx.user[1]
            //		                      << " userint[7]: " << (int16_t)acqhdr.user_int[7]);
        }
        else {
            GDEBUG_STREAM( " e1: " << std::setw( 3 ) << (int16_t)acqhdr.idx.kspace_encode_step_1
                << " set: " << (int16_t)acqhdr.idx.set
                << " con: " << (int16_t)acqhdr.idx.contrast
                << " iLastLine: " << iLastLine
                << " nav: " << acqhdr.isFlagSet( ISMRMRD::ISMRMRD_ACQ_IS_RTFEEDBACK_DATA ) );
            //		                      << " user[1]: " << (int16_t)acqhdr.idx.user[1]
            //		                      << " userint[0]: " << (int16_t)acqhdr.user_int[0]
            //		                      << " userint[1]: " << (int16_t)acqhdr.user_int[1]
            //		                      << " userint[2]: " << (int16_t)acqhdr.user_int[2]
            //		                      << " userint[3]: " << (int16_t)acqhdr.user_int[3]
            //		                      << " userint[4]: " << (int16_t)acqhdr.user_int[4]
            //		                      << " userint[5]: " << (int16_t)acqhdr.user_int[5]
            //		                      << " userint[6]: " << (int16_t)acqhdr.user_int[6]
            //		                      << " userint[7]: " << (int16_t)acqhdr.user_int[7]);
        }

        // 2016-10-04 Kelvin: Discard navigator lines
        if (acqhdr.isFlagSet( ISMRMRD::ISMRMRD_ACQ_IS_RTFEEDBACK_DATA ))
        {
            return GADGET_OK;
        }
    }

    //Now pass on image
    if (this->next()->putq( m1 ) < 0) {
        return GADGET_FAIL;
    }

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE( FlagSashaHCGadget )