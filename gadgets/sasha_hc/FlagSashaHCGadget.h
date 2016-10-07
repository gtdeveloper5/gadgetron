//ThresholdGadget.h

#ifndef FLAGSASHAHCGADGET_H
#define FLAGSASHAHCGADGET_H

#include "sashahclib_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include <complex>
#include <ismrmrd/ismrmrd.h>

namespace Gadgetron
{

    class EXPORTSASHAHC FlagSashaHCGadget :
        public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:
        GADGET_DECLARE( FlagSashaHCGadget )

    protected:
        virtual int process( GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1,
            GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 );

        virtual int process_config( ACE_Message_Block* mb );

        //  float threshold_level_;

    };

}
#endif //FLAGSASHAHCGADGET_H