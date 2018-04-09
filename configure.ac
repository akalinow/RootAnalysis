dnl Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
AC_INIT(RootAnalysis, 1.0)
AM_INIT_AUTOMAKE
: ${CXXFLAGS=-Werror -O3 -std=c++11 -fopenmp}
AC_PROG_LIBTOOL
LT_INIT

AC_PROG_CC
AC_PROG_CXX
AC_LANG(C++)

AC_CONFIG_MACRO_DIR([m4])

#AC_CHECK_HEADER_STDBOOL
ROOT_PATH(6.06, , AC_MSG_ERROR(Your ROOT version is too old))
BOOST_REQUIRE(1.3, AC_MSG_ERROR(Your BOOST version is too old))
BOOST_TOKENIZER
BOOST_HASH
BOOST_REGEX

AC_CONFIG_FILES(Makefile
                test/Makefile
                src/Makefile
                OTFAnalysis/Makefile
		CPwithPytia8/Makefile
		HTTAnalysis/Makefile
    SVFitAnalysis/Makefile)

AC_OUTPUT
