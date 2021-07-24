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
    parallel_framework.h
    primer_screen.cpp
    primer_screen.h
    quality.h
    reads.cpp
    reads.h
    rlstring.h
    runner.h
    suffix_array.cpp
    suffix_array.h
    suffix_array_builder.cpp
    suffix_array_builder.h
    utils.cpp
    utils.h
  )
target_compile_options(libsiga PUBLIC ${OpenMP_CXX_FLAGS} ${Log4cxx_CFLAGS} ${MLPack_CFLAGS} ${RapidJSON_CFLAGS} ${TCMalloc_CFLAGS})
target_include_directories(libsiga PUBLIC ${CMAKE_BINARY_DIR})
target_link_libraries(libsiga 3rdparty OpenMP::OpenMP_CXX Boost::filesystem Boost::iostreams Boost::program_options ${Log4cxx_LDFLAGS} ${MLPack_LDFLAGS} ${RapidJSON_LDFLAGS} ${TCMalloc_LDFLAGS})

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
