# Default build configuration.
# Can be overriden with `-DCMAKE_BUILD_TYPE=<Debug|RelWithDebInfo|Release|...>` when calling CMake,
# or by using an environment variable `export CMAKE_BUILD_TYPE=...`.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE set to [${CMAKE_BUILD_TYPE}]")

# Default C++ standard for all the targets in the project.
# Can be overriden with `-DCMAKE_CXX_STANDARD=...` when calling CMake (but should be at least the default value).
if(NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
endif()
message(STATUS "CMAKE_CXX_STANDARD set to [${CMAKE_CXX_STANDARD}]")
