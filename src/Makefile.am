noinst_LTLIBRARIES=libsiga.la
libsiga_la_CPPFLAGS=\
            -I$(top_srcdir)/3rdparty \
            ${BOOST_CPPFLAGS}
libsiga_la_LIBADD=\
            ${top_builddir}/3rdparty/lib3rdparty.la
libsiga_la_CXXFLAGS=\
            ${LOG4CXX_CFLAGS} \
            ${OPENMP_CFLAGS} \
            ${RAPIDJSON_CFLAGS} \
            ${TCMALLOC_CFLAGS}
libsiga_la_LDFLAGS=\
            ${BOOST_LDFLAGS} \
            ${LOG4CXX_LIBS} \
            ${OPENMP_LIBS} \
            ${RAPIDJSON_LIBS} \
            ${TCMALLOC_LIBS} \
            ${BOOST_FILESYSTEM_LIB} \
            ${BOOST_IOSTREAMS_LIB} \
            ${BOOST_SYSTEM_LIB}
libsiga_la_SOURCES=\
            alphabet.h \
            asqg.cpp \
            asqg.h \
            bigraph.cpp \
            bigraph.h \
            bigraph_search.cpp \
            bigraph_search.h \
            bigraph_visitors.cpp \
            bigraph_visitors.h \
            bwt.cpp \
            bwt.h \
            constant.h \
            correct_processor.cpp \
            correct_processor.h \
            coord.cpp \
            coord.h \
            fmindex.cpp \
            fmindex.h \
            kmerdistr.cpp \
            kmerdistr.h \
            kseq.cpp \
            kseq.h \
            mkqs.h \
            overlap_builder.cpp \
            overlap_builder.h \
            primer_screen.cpp \
            primer_screen.h \
            quality.h \
            reads.cpp \
            reads.h \
            rlstring.h \
            runner.h \
            sequence_process_framework.h \
            suffix_array.cpp \
            suffix_array.h \
            suffix_array_builder.cpp \
            suffix_array_builder.h \
            utils.cpp \
            utils.h

bin_PROGRAMS=siga
siga_SOURCES=\
            assembler.cpp \
            correct.cpp \
            indexer.cpp \
            match.cpp \
            overlap.cpp \
            preprocess.cpp \
            preqc.cpp \
            rmdup.cpp \
            subgraph.cpp \
            main.cpp
siga_LDADD=libsiga.la
siga_CXXFLAGS=\
            ${LOG4CXX_CFLAGS} \
            ${OPENMP_CFLAGS} \
            ${RAPIDJSON_CFLAGS} \
            ${TCMALLOC_CFLAGS}

if HELP2MAN
man_MANS=siga.1
siga.1: siga${EXEEXT}
	help2man -o siga.1 ./siga${EXEEXT}
clean-local:
	rm -rf siga.1
else
man_MANS=
endif !HELP2MAN

EXTRA_DIST=log4cxx.properties \
            CMakeLists.txt
