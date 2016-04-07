
#include "FatWaterSeparationGadget.h"

namespace Gadgetron { 

FatWaterSeparationGadget::FatWaterSeparationGadget() : BaseClass()
{
}

FatWaterSeparationGadget::~FatWaterSeparationGadget()
{
}

int FatWaterSeparationGadget::process_config(ACE_Message_Block* mb)
{
    GADGET_CHECK_RETURN(BaseClass::process_config(mb)==GADGET_OK, GADGET_FAIL);
    return GADGET_OK;
}

int FatWaterSeparationGadget::processImageBuffer(ImageBufferType& ori)
{
    GDEBUG_CONDITION_STREAM(verbose.value(), "FatWaterSeparationGadget::processImageBuffer(...) starts ... ");

    std::vector<std::string> processStr;
    std::vector<std::string> dataRole;

    boost::shared_ptr< std::vector<size_t> > dims = ori.get_dimensions();
    GDEBUG_CONDITION_STREAM(verbose.value(), "[Cha Slice E2 Con Phase Rep Set] = [" << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " << (*dims)[3] << " " << (*dims)[4]  << " " << (*dims)[5] << " " << (*dims)[6] << "]");

    // --------------------------------------------------------------------------------
    // ori
    // --------------------------------------------------------------------------------
    if (send_multi_echo_images.value())
    {
        GADGET_CHECK_RETURN(this->send_out_images(ori, processing_result_image_series_num.value(), processStr, dataRole)==GADGET_OK, GADGET_FAIL);
    }

    if (send_water_images.value() || send_fat_images.value() || send_t2_star_map.value() || send_field_map.value())
    {
        ImageBufferType water, fat, t2s_map, field_map;

        /*GADGET_CHECK_RETURN(this->performT2W(ori, oriMagT2W, oriMagT2WNoSCC, oriMagPD), GADGET_FAIL);

        processStr.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);

        if ( send_ori_T2W_ )
        {
            processStr.clear();
            processStr.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);

            dataRole.clear();
            dataRole.push_back(GADGETRON_IMAGE_SURFACECOILCORRECTION);
            dataRole.push_back(GADGETRON_IMAGE_T2W);

            GADGET_CHECK_RETURN(this->sendOutImages(oriMagT2W, image_series_num_+1, processStr, dataRole), GADGET_FAIL);

            if (send_no_scc_T2W_)
            {
                processStr.clear();

                dataRole.clear();
                dataRole.push_back(GADGETRON_IMAGE_T2W);

                GADGET_CHECK_RETURN(this->sendOutImages(oriMagT2WNoSCC, image_series_num_ + 9, processStr, dataRole), GADGET_FAIL);
            }
        }

        if ( send_ori_PD_ )
        {
            processStr.clear();

            dataRole.clear();
            dataRole.push_back(GADGETRON_IMAGE_PD);

            GADGET_CHECK_RETURN(this->sendOutImages(oriMagPD, image_series_num_+2, processStr, dataRole), GADGET_FAIL);
        }

        GADGET_CHECK_RETURN(this->releaseImageBuffer(oriMagT2W), GADGET_FAIL);
        GADGET_CHECK_RETURN(this->releaseImageBuffer(oriMagT2WNoSCC), GADGET_FAIL);
        GADGET_CHECK_RETURN(this->releaseImageBuffer(oriMagPD), GADGET_FAIL);*/
    }

    GDEBUG_CONDITION_STREAM(verbose.value(), "FatWaterSeparationGadget::process(...) ends ... ");

    return GADGET_OK;
}

int FatWaterSeparationGadget::perform_fat_water(ImageBufferType& input, ImageBufferType& water, ImageBufferType& fat, ImageBufferType& t2s_map, ImageBufferType& field_map)
{
    size_t CHA = input.get_size(0);
    size_t SLC = input.get_size(1);
    size_t CON = input.get_size(2);
    size_t PHS = input.get_size(3);
    size_t REP = input.get_size(4);
    size_t SET = input.get_size(5);
    size_t AVE = input.get_size(6);

    //std::vector<size_t> dim;
    //input.get_dimensions(dim);

    //dim[5] = 1; // only one set is outputted

    //magT2W.create(dim);
    //this->fillWithNULL(magT2W);

    //magT2WNoSCC.create(dim);
    //this->fillWithNULL(magT2WNoSCC);

    //magPD.create(dim);
    //this->fillWithNULL(magPD);

    //size_t cha, slc, con, phs, rep, set, ave;

    //ImageContainerType inputContainer, magT2WContainer, magT2WNoSCCContainer, magPDContainer;

    //psir_reconer_.perform_PSIR_ = false;
    //psir_reconer_.compute_PSIR_windowing_  = false;

    //std::vector<size_t> cols(SET, CON);
    //inputContainer.create(cols, false);

    //for ( cha=0; cha<CHA; cha++ )
    //{
    //    for ( ave=0; ave<AVE; ave++ )
    //    {
    //        for ( slc=0; slc<SLC; slc++ )
    //        {
    //            ImageType gmap;
    //            bool gMapFound = this->findGFactorMap(slc, gmap);
    //            float gMap_scale_factor = 1.0f;
    //            if (gMapFound)
    //            {
    //                GDEBUG_CONDITION_STREAM(verbose.value(), "Gfactor map for slc - " << slc << " is found ... ");

    //                if (gmap.attrib_.length(GADGETRON_IMAGE_SCALE_RATIO) > 0)
    //                {
    //                    gMap_scale_factor = (float)gmap.attrib_.as_double(GADGETRON_IMAGE_SCALE_RATIO, 0);
    //                    GDEBUG_CONDITION_STREAM(verbose.value(), "Gfactor map for slc - " << slc << " has scale factor " << gMap_scale_factor);

    //                    Gadgetron::scal((float)(1.0 / gMap_scale_factor), gmap);
    //                }

    //                if (!debugFolder_fullPath_.empty())
    //                {
    //                    std::ostringstream ostr;
    //                    ostr << "gfactor_map_slc" << slc;
    //                    if (!debugFolder_fullPath_.empty()) gt_exporter_.exportArrayComplex(gmap, debugFolder_fullPath_ + ostr.str());
    //                }
    //            }
    //            else
    //            {
    //                GDEBUG_CONDITION_STREAM(true, "Gfactor map for slc - " << slc << " is NOT found ... ");
    //            }

    //            for ( phs=0; phs<PHS; phs++ )
    //            {
    //                for ( rep=0; rep<REP; rep++ )
    //                {
    //                    // compute PSIR for every SET/CON

    //                    for ( con=0; con<CON; con++ )
    //                    {
    //                        for ( set=0; set<SET; set++ )
    //                        {
    //                            inputContainer.set(input(cha, slc, con, phs, rep, set, ave), set, con);
    //                        }
    //                    }

    //                    GADGET_CHECK_RETURN_FALSE(this->performT2W(inputContainer, gmap, magT2WContainer, magT2WNoSCCContainer, magPDContainer));

    //                    for ( con=0; con<CON; con++ )
    //                    {
    //                        magT2W(cha, slc, con, phs, rep, 0, ave) = &magT2WContainer(0, con);
    //                        magT2W(cha, slc, con, phs, rep, 0, ave)->header_ = input(cha, slc, con, phs, rep, 0, ave)->header_;
    //                        magT2W(cha, slc, con, phs, rep, 0, ave)->attrib_ = input(cha, slc, con, phs, rep, 0, ave)->attrib_;
    //                        magT2W(cha, slc, con, phs, rep, 0, ave)->attrib_.set(GADGETRON_IMAGE_SCALE_RATIO, psir_reconer_.scale_factor_after_SCC_);

    //                        magT2WNoSCC(cha, slc, con, phs, rep, 0, ave) = &magT2WNoSCCContainer(0, con);
    //                        magT2WNoSCC(cha, slc, con, phs, rep, 0, ave)->header_ = input(cha, slc, con, phs, rep, 0, ave)->header_;
    //                        magT2WNoSCC(cha, slc, con, phs, rep, 0, ave)->attrib_ = input(cha, slc, con, phs, rep, 0, ave)->attrib_;

    //                        magPD(cha, slc, con, phs, rep, 0, ave) = &magPDContainer(0, con);
    //                        magPD(cha, slc, con, phs, rep, 0, ave)->header_ = input(cha, slc, con, phs, rep, 0, ave)->header_;
    //                        magPD(cha, slc, con, phs, rep, 0, ave)->attrib_ = input(cha, slc, con, phs, rep, 0, ave)->attrib_;
    //                    }

    //                    magT2WContainer.delete_data_on_destruct(false);
    //                    magT2WNoSCCContainer.delete_data_on_destruct(false);
    //                    magPDContainer.delete_data_on_destruct(false);
    //                }
    //            }
    //        }
    //    }
    //}

    return GADGET_OK;
}

int FatWaterSeparationGadget::perform_fat_water(ImageContainerType& input, ImageContainerType& water, ImageContainerType& fat, ImageContainerType& t2s_map, ImageContainerType& field_map)
{
    try
    {
        /*magT2W.clear();
        magT2WNoSCC.clear();
        magPD.clear();

        if (gmap.dimensions_equal(input(0, 0)))
        {
            psir_reconer_.gmap_ = gmap;
        }
        else
        {
            psir_reconer_.gmap_.clear();
        }

        GADGET_CHECK_RETURN_FALSE(psir_reconer_.performPSIRRecon(input));

        GADGET_CHECK_RETURN_FALSE(magT2W.copyFrom(psir_reconer_.magIR_));
        GADGET_CHECK_RETURN_FALSE(magT2WNoSCC.copyFrom(psir_reconer_.magIR_without_scc_));
        GADGET_CHECK_RETURN_FALSE(magPD.copyFrom(psir_reconer_.magPD_));*/
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in FatWaterSeparationGadget::perform_fat_water(ImageContainerType& input, ...) ... ");
        return false;
    }

    return true;
}

int FatWaterSeparationGadget::close(unsigned long flags)
{
    GDEBUG_CONDITION_STREAM(true, "FatWaterSeparationGadget - close(flags) : " << flags);

    if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(FatWaterSeparationGadget)

}
