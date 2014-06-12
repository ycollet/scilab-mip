# Copyright (c) 2004-2007 University of Geneva, HEC, Logilab
# 
# OBOE is published under the Common Public License.
#
# Authors :
# Nidhi Sawhney <nsawhney@yahoo.com>
# The OBOE team
#
# 

include_HEADERS = AccpmVector.hpp AccpmGenMatrix.hpp AccpmDynMatrix.hpp AccpmBlasInterface.hpp

#check_PROGRAMS = Test TestSolve
check_PROGRAMS = Test

#Uncomment the following to run test with make check
#TESTS = $(check_PROGRAMS)

lib_LIBRARIES = libaccpmla.a

libaccpmla_a_SOURCES = AccpmVector.cpp AccpmVector.hpp AccpmGenMatrix.cpp AccpmGenMatrix.hpp  \
	AccpmDynMatrix.cpp  AccpmDynMatrix.hpp \
	AccpmLASolve.hpp AccpmLASolve.cpp AccpmBlasInterface.hpp AccpmBlasInterface.cpp

Test_SOURCES = MatrixTest.cpp 
Test_LDADD = libaccpmla.a $(OBOE_LDADD) 
Test_LDFLAGS = -L./ $(OBOE_LDFLAGS) 

#TestSolve_SOURCES = LinSolveTest.cpp
#TestSolve_LDADD = libaccpmla.a $(OBOE_LDADD) 
#TestSolve_LDFLAGS = -L./ $(OBOE_LDFLAGS) 

INCLUDES = -I$(top_srcdir)/src/Utilities \
	$(OBOE_INCLUDES)
