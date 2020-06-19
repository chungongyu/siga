check_PROGRAMS=index_test preprocess_test

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
            -I$(top_srcdir)/src \
            ${BOOST_CPPFLAGS}
index_test_LDADD=\
            ${top_builddir}/src/libsiga.la
index_test_LDFLAGS=\
            ${BOOST_LDFLAGS} \
            ${BOOST_UNIT_TEST_FRAMEWORK_LIB}
index_test_SOURCES=\
            index_test.cpp
