# Specify the source files for the library
set(SOURCES
  ProfileMakerCore.C
  ProfileMaker.C
  canvas_margins.C
  dEdx_functions.C
)

set(HDRFILES
  ProfileMakerCore.h
  ProfileMaker.h
  canvas_margins.h
  dEdx_functions.h
)

# Create the shared library
add_library(ProfileMakers SHARED ${SOURCES})

# Specify output directory for the shared library
set_target_properties(ProfileMakers PROPERTIES
  PUBLIC_HEADER "${HDRFILES}"
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib
)

target_link_libraries(ProfileMakers
PUBLIC
  ${ROOT_LIBRARIES}
)

install(TARGETS ProfileMakers
    EXPORT LArProfMaker-targets
    PUBLIC_HEADER DESTINATION include/ProfileMakers COMPONENT Development)
