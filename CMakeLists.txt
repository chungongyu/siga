cmake_minimum_required(VERSION 3.11)

project(siga VERSION 0.10)

set(PROJECT_URL "http://bioinfo.ict.ac.cn/siga")
set(PACKAGE_BUGREPORT "siga@ict.ac.cn")

set(CMAKE_CXX_STANDARD 11)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Boost COMPONENTS filesystem iostreams program_options unit_test_framework REQUIRED)
find_package(OpenMP QUIET)
find_package(RapidJSON REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(LOG4CXX REQUIRED liblog4cxx>=0.10.0)
pkg_check_modules(TCMalloc QUIET libtcmalloc>=2.1)

find_program(HELP2MAN NAMES help2man)
find_program(BUNDLE NAMES bundle)

configure_file(siga.h.in config.h @ONLY)
enable_testing()

add_subdirectory(3rdparty)
add_subdirectory(docs)
add_subdirectory(src)
add_subdirectory(test)
