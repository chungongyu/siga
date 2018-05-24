bin_PROGRAMS=arcs
arcs_CPPFLAGS=-I$(top_srcdir)/3rdparty
arcs_LDADD=${top_builddir}/3rdparty/lib3rdparty.a
arcs_SOURCES=assembler.cpp \
			 asqg.cpp \
			 asqg.h \
			 bigraph.cpp \
			 bigraph.h \
			 bigraph_visitors.cpp \
			 bigraph_visitors.h \
			 bwt.cpp \
			 bwt.h \
			 constant.h \
			 correct.cpp \
			 correct_processor.cpp \
			 correct_processor.h \
			 coord.cpp \
			 coord.h \
			 fmindex.cpp \
			 fmindex.h \
			 indexer.cpp \
			 kmerdistr.cpp \
			 kmerdistr.h \
			 kseq.cpp \
			 kseq.h \
			 match.cpp \
			 mkqs.h \
			 overlap.cpp \
			 overlap_builder.cpp \
			 overlap_builder.h \
			 preprocess.cpp \
			 preqc.cpp \
			 primer_screen.cpp \
			 primer_screen.h \
			 reads.cpp \
			 reads.h \
			 rmdup.cpp \
			 runner.h \
			 sequence_process_framework.h \
			 subgraph.cpp \
			 suffix_array.cpp \
			 suffix_array.h \
			 suffix_array_builder.cpp \
			 suffix_array_builder.h \
			 utils.h \
			 main.cpp
arcs_CXXFLAGS=${LOG4CXX_CFLAGS} \
			  ${OPENMP_CFLAGS} \
			 ${TCMALLOC_CFLAGS}
arcs_LDFLAGS=${LOG4CXX_LIBS} \
			 ${OPENMP_LIBS} \
			 ${TCMALLOC_LIBS} \
			 -lboost_filesystem \
			 -lboost_iostreams \
			 -lboost_system

EXTRA_DIST=log4cxx.properties

