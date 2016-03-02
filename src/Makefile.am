bin_PROGRAMS=moga
moga_SOURCES=assembler.cpp \
			 bwt.cpp \
			 bwt.h \
			 constant.cpp \
			 constant.h \
			 fmindex.cpp \
			 fmindex.h \
			 gapfiller.cpp \
			 indexer.cpp \
			 kseq.cpp \
			 kseq.h \
			 preprocess.cpp \
			 primer_screen.cpp \
			 primer_screen.h \
			 runner.h \
			 scaffold.cpp \
			 suffix_array.cpp \
			 suffix_array.h \
			 suffix_array_builder.cpp \
			 suffix_array_builder.h \
			 utils.cpp \
			 utils.h \
			 main.cpp
moga_CXXFLAGS=${LOG4CXX_CFLAGS} \
			  ${DIVSUFSORT_CFLAGS}
moga_LDFLAGS=${LOG4CXX_LIBS} \
			 ${DIVSUFSORT_LIBS} \
			 -lboost_filesystem \
			 -lboost_system

EXTRA_DIST=log4cxx.properties

