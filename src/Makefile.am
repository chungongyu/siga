bin_PROGRAMS=moga
moga_SOURCES=assembler.cpp \
			 bwt.cpp \
			 bwt.h \
			 constant.cpp \
			 constant.h \
			 indexer.cpp \
			 gapfiller.cpp \
			 kseq.cpp \
			 kseq.h \
			 preprocess.cpp \
			 scaffold.cpp \
			 suffix_array.cpp \
			 suffix_array.h \
			 suffix_array_builder.cpp \
			 suffix_array_builder.h \
			 runner.h \
			 main.cpp
moga_CXXFLAGS=${LOG4CXX_CFLAGS} \
			  ${DIVSUFSORT_CFLAGS}
moga_LDFLAGS=${LOG4CXX_LIBS} \
			 ${DIVSUFSORT_LIBS} \
			 -lboost_filesystem \
			 -lboost_system

EXTRA_DIST=log4cxx.properties

