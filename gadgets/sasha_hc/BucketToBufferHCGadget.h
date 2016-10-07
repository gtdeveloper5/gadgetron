#ifndef BUCKETTOBUFFERHC_H
#define BUCKETTOBUFFERHC_H

#include "sashahclib_export.h"
#include "BucketToBufferGadget.h"

namespace Gadgetron {

    class EXPORTSASHAHC BucketToBufferHCGadget :
        public BucketToBufferGadget
    {
    public:
        GADGET_DECLARE( BucketToBufferHCGadget );

        BucketToBufferHCGadget();
        virtual ~BucketToBufferHCGadget();

    protected:
        virtual int process( GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1 );
    };
}
#endif //BUCKETTOBUFFERHC_H
