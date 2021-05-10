cmake_minimum_required(VERSION 3.10)

project(siga VERSION 0.10)

set(PROJECT_URL "http://bioinfo.ict.ac.cn/siga")
set(PROJECT_BUGREPORT "siga@ict.ac.cn")

set(CMAKE_CXX_STANDARD 11)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Boost COMPONENTS filesystem iostreams program_options unit_test_framework REQUIRED)
find_package(OpenMP QUIET)
find_package(RapidJSON REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(Log4cxx REQUIRED IMPORTED_TARGET liblog4cxx>=0.10.0)
pkg_check_modules(MLPack QUIET IMPORTED_TARGET mlpack>=3.2.1)
pkg_check_modules(RapidJSON REQUIRED IMPORTED_TARGET RapidJSON>=1.0.0)
pkg_check_modules(TCMalloc QUIET IMPORTED_TARGET libtcmalloc>=2.1)

find_program(HELP2MAN NAMES help2man)
find_program(BUNDLE NAMES bundle)

add_definitions(-DHAVE_CONFIG_H)
configure_file(siga.h.in config.h @ONLY)
enable_testing()

add_subdirectory(3rdparty)
add_subdirectory(docs)
add_subdirectory(src)
add_subdirectory(test)
