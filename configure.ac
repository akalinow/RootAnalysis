dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(RootAnalysis, 1.0)
: ${CXXFLAGS=-Werror -O3 -std=c++0x -fopenmp}
AC_PROG_CC
AC_PROG_CXX

AC_CANONICAL_SYSTEM
AC_CHECK_FUNCS([mkdir])
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_PROG_RANLIB

ROOT_PATH

dnl Initialize automake
AM_INIT_AUTOMAKE

dnl this allows us specify individual liking flags for each target
AM_PROG_CC_C_O 

dnl Initialize Libtool
LT_INIT

dnl Check if Libtool is present
dnl Libtool is used for building share libraries 
AC_PROG_LIBTOOL

AC_CONFIG_FILES(Makefile
                test/Makefile
                src/Makefile
                OTFAnalysis/Makefile
                PACAnalysis/Makefile
		CPwithPytia8/Makefile
		HTTAnalysis/Makefile)
AC_OUTPUT
