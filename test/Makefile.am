check_PROGRAMS=siga_test
siga_test_CPPFLAGS=\
            -I$(top_srcdir)/src \
            ${BOOST_CPPFLAGS}
siga_test_LDADD=\
            ${top_builddir}/src/libsiga.la
siga_test_LDFLAGS=\
            ${BOOST_LDFLAGS} \
            ${BOOST_UNIT_TEST_FRAMEWORK_LIB}
siga_test_SOURCES=alphabet_test.cpp
