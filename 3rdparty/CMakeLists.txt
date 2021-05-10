add_library(3rdparty
    bcr.h
    bcr.c
    sais.hxx
    ssw.h
    ssw.c
    ssw_cpp.cpp
    ssw_cpp.h
  )
target_include_directories(3rdparty PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
