cmake_minimum_required(VERSION 3.19)
project(HiPOP VERSION 1.0.0)


find_package(OpenMP REQUIRED)


# ------------- SRC -------------

set(HIPOP_HEADERS
  include/hipop/create.h
  include/hipop/graph.h
  include/hipop/shortest_path.h
)

add_subdirectory(src)

# -------------------------------



# ------------- TESTS -------------

include(CTest)
add_subdirectory(tests)

# ---------------------------------




# ------------- INSTALL -------------

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include
    FILES_MATCHING PATTERN "*.h"
    PATTERN "*_p.h" EXCLUDE
)


install(TARGETS hipoplib
        EXPORT HipopTarget
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin COMPONENT RuntimeExecutables
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib COMPONENT RuntimeLibraries
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib COMPONENT Development
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/hipop
)


install(EXPORT HipopTarget DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/hipop)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(${CMAKE_BINARY_DIR}/Hipop/HipopConfigVersion.cmake
  VERSION ${HiPOP_VERSION}
  COMPATIBILITY ExactVersion)

configure_package_config_file(${HiPOP_SOURCE_DIR}/cmake/HipopConfig.cmake.in
${CMAKE_BINARY_DIR}/Hipop/HipopConfig.cmake
INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/hipop)

install(
  FILES
  ${CMAKE_BINARY_DIR}/Hipop/HipopConfig.cmake
  ${CMAKE_BINARY_DIR}/Hipop/HipopConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/hipop
  COMPONENT Development)

# --------------------------------