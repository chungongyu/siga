add_library(3rdparty
    bcr.h
    bcr.c
  )
target_include_directories(3rdparty PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
