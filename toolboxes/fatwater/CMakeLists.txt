if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_FATWATER__)
endif ()

include_directories(
  ${Boost_INCLUDE_DIR} 
  ${ARMADILLO_INCLUDE_DIRS} 
  ${MKL_INCLUDE_DIR} 
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
  )

set( gadgetron_fatwater_header_files 
    fatwater_export.h 
    fatwater.h 
    fatWaterMultiPeakOperator.h 
    fatWaterMultiPeakR2StarOperator.h 
    )

set( gadgetron_fatwater_src_files fatwater.cpp )

add_library(gadgetron_toolbox_fatwater SHARED 
  ${gadgetron_fatwater_header_files} 
  ${gadgetron_fatwater_src_files} 
  )

set_target_properties(gadgetron_toolbox_fatwater PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_fatwater
  gadgetron_toolbox_log
  gadgetron_toolbox_cpucore 
  gadgetron_toolbox_cpucore_math 
  ${ARMADILLO_LIBRARIES} 
  ${Boost_LIBRARIES}
  )

install(TARGETS gadgetron_toolbox_fatwater DESTINATION lib COMPONENT main)

install(FILES
  ${gadgetron_fatwater_header_files} 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
