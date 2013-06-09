
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

ADD_EXECUTABLE(openMVG_sample_eightPointsEssential eight_points_essential.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_eightPointsEssential
  openMVG_image
  openMVG_multiview
  ${OPENMVG_LIBRARY_DEPENDENCIES})
  
ADD_EXECUTABLE(openMVG_sample_eightPointsEssential_hardcoded eight_points_essential_hardcoded.cpp)
TARGET_LINK_LIBRARIES(openMVG_sample_eightPointsEssential_hardcoded
  openMVG_image
  openMVG_multiview
  ${OPENMVG_LIBRARY_DEPENDENCIES})

