
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_sixPointsResection six_points_resection.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_sixPointsResection
  openMVG_image
  openMVG_multiview
  ${OPENMVG_LIBRARY_DEPENDENCIES})
  
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_sixPointsResection_hardcoded six_points_resection_hardcoded.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_sixPointsResection_hardcoded
  openMVG_image
  openMVG_multiview
  ${OPENMVG_LIBRARY_DEPENDENCIES})

