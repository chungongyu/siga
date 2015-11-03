bin_PROGRAMS=moga
moga_SOURCES=assembler.cpp \
			 assembler.h \
			 constant.cpp \
			 constant.h \
			 indexer.cpp \
			 indexer.h \
			 preprocess.cpp \
			 preprocess.h \
			 runner.h \
			 main.cpp
moga_CXXFLAGS=${LOG4CXX_CFLAGS}
moga_LDFLAGS=${LOG4CXX_LIBS} \
			 -lboost_filesystem \
			 -lboost_system

EXTRA_DIST=log4cxx.properties

