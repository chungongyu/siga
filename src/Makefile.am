bin_PROGRAMS=moga
moga_SOURCES=constant.cpp \
			 constant.h \
			 runner.h \
			 main.cpp
moga_CXXFLAGS=${LOG4CXX_CFLAGS}
moga_LDFLAGS=${LOG4CXX_LIBS} \
			 -lboost_filesystem \
			 -lboost_system

EXTRA_DIST=log4cxx.properties

