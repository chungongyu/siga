check_PROGRAMS=parallel_test preprocess_test index_test overlap_test assemble_test utils_test

parallel_test_CPPFLAGS=\
            -I$(top_srcdir)/src \
	        ${OPENMP_CFLAGS} \
            ${BOOST_CPPFLAGS}
parallel_test_LDADD=\
            ${top_builddir}/src/libsiga.la
parallel_test_CXXFLAGS=\
	        ${OPENMP_CFLAGS}
parallel_test_LDFLAGS=\
	        ${OPENMP_LIBS} \
            ${BOOST_LDFLAGS} \
            ${BOOST_UNIT_TEST_FRAMEWORK_LIB}
parallel_test_SOURCES=\
            parallel_test.cpp

preprocess_test_CPPFLAGS=\
            -I$(top_srcdir)/src \
            ${BOOST_CPPFLAGS}
preprocess_test_LDADD=\
            ${top_builddir}/src/libsiga.la
preprocess_test_LDFLAGS=\
            ${BOOST_LDFLAGS} \
            ${BOOST_UNIT_TEST_FRAMEWORK_LIB}
preprocess_test_SOURCES=\
            preprocess_test.cpp

index_test_CPPFLAGS=\
            -I$(top_srcdir)/3rdparty \
            -I$(top_srcdir)/src \
            ${BOOST_CPPFLAGS}
index_test_LDADD=\
            ${top_builddir}/src/libsiga.la
index_test_LDFLAGS=\
            ${BOOST_LDFLAGS} \
            ${BOOST_UNIT_TEST_FRAMEWORK_LIB}
index_test_SOURCES=\
            index_test.cpp

overlap_test_CPPFLAGS=\
            -I$(top_srcdir)/src \
            ${BOOST_CPPFLAGS}
overlap_test_LDADD=\
            ${top_builddir}/src/libsiga.la
overlap_test_LDFLAGS=\
            ${BOOST_LDFLAGS} \
            ${BOOST_UNIT_TEST_FRAMEWORK_LIB}
overlap_test_SOURCES=\
            overlap_test.cpp

assemble_test_CPPFLAGS=\
            -I$(top_srcdir)/src \
            ${BOOST_CPPFLAGS}
assemble_test_LDADD=\
            ${top_builddir}/src/libsiga.la
assemble_test_LDFLAGS=\
            ${BOOST_LDFLAGS} \
            ${BOOST_UNIT_TEST_FRAMEWORK_LIB}
assemble_test_SOURCES=\
            assemble_test.cpp

utils_test_CPPFLAGS=\
            -I$(top_srcdir)/src \
            ${BOOST_CPPFLAGS}
utils_test_LDADD=\
            ${top_builddir}/src/libsiga.la
utils_test_LDFLAGS=\
            ${BOOST_LDFLAGS} \
            ${BOOST_UNIT_TEST_FRAMEWORK_LIB}
utils_test_SOURCES=\
            utils_test.cpp

TESTS=${check_PROGRAMS}

EXTRA_DIST=CMakeLists.txt
