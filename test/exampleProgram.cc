/*
 ============================================================================
 Name        : exampleProgram.c
 Author      : Artur Kalinowski
 Version     :
 Copyright   : Your copyright notice
 Description : Uses shared library to print greeting
               To run the resulting executable the LD_LIBRARY_PATH must be
               set to ${project_loc}/libRootAnalysis/.libs
               Alternatively, libtool creates a wrapper shell script in the
               build directory of this program which can be used to run it.
               Here the script will be called exampleProgram.
 ============================================================================
 */

#include <iostream>

#include "TreeAnalyzer.h"

int main(void) {

  std::cout<<"Test"<<std::endl;

  return 0;
}
