add_library(libsiga
    alphabet.h
    asqg.cpp
    asqg.h
    bigraph.cpp
    bigraph.h
    bigraph_search.cpp
    bigraph_search.h
    bigraph_visitors.cpp
    bigraph_visitors.h
    bwt.cpp
    bwt.h
    constant.h
    correct_processor.cpp
    correct_processor.h
    coord.cpp
    coord.h
    fmindex.cpp
    fmindex.h
    kmerdistr.cpp
    kmerdistr.h
    kseq.cpp
    kseq.h
    mkqs.h
    overlap_builder.cpp
    overlap_builder.h
    primer_screen.cpp
    primer_screen.h
    quality.h
    reads.cpp
    reads.h
    rlstring.h
    runner.h
    sequence_process_framework.h
    suffix_array.cpp
    suffix_array.h
    suffix_array_builder.cpp
    suffix_array_builder.h
    utils.cpp
    utils.h
  )
target_include_directories(libsiga PUBLIC ${CMAKE_BINARY_DIR} ${LOG4CXX_CFLAGS})
target_link_libraries(libsiga 3rdparty ${LOG4CXX_LDFLAGS} Boost::filesystem Boost::iostreams Boost::program_options)
if(TCMalloc_FOUND)
  target_include_directories(libsiga PUBLIC ${TCMalloc_CFLAGS})
  target_link_libraries(libsiga ${TCMalloc_LDFLAGS})
endif(TCMalloc_FOUND)

add_executable(siga
    assembler.cpp
    correct.cpp
    indexer.cpp
    match.cpp
    overlap.cpp
    preprocess.cpp
    preqc.cpp
    rmdup.cpp
    subgraph.cpp
    main.cpp
  )
target_link_libraries(siga libsiga)
target_include_directories(siga PRIVATE ${CMAKE_BINARY_DIR})

if(HELP2MAN)
  add_custom_target(man ALL
      COMMAND help2man -o siga.1 $<TARGET_FILE:siga>
      DEPENDS siga
    )
endif(HELP2MAN)
