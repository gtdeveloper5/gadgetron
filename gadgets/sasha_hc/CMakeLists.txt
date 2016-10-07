if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_SASHA_HC__)
endif ()

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/gadgets/pmri
  ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${ARMADILLO_INCLUDE_DIRS}
  ${ISMRMRD_INCLUDE_DIR}
  )

add_library(gadgetron_sasha_hc SHARED 
    SASHA-HC.xml
    FlagSashaHCGadget.h      FlagSashaHCGadget.cpp
    BucketToBufferHCGadget.h BucketToBufferHCGadget.cpp
    sashahclib_export.h)

set_target_properties(gadgetron_sasha_hc PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_sasha_hc
  gadgetron_gadgetbase
  gadgetron_toolbox_log
  gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_mricore
  ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

install (FILES  FlagSashaHCGadget.h
                BucketToBufferHCGadget.h
                sashahclib_export.h
                DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES SASHA-HC.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

install (TARGETS gadgetron_sasha_hc DESTINATION lib COMPONENT main)

#add_subdirectory(config)
