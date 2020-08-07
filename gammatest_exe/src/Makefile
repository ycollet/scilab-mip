# kdtree.h
# Author: Y. Collette
# Date: 19/11/2005

GCC = gcc -O2 -Wall
GPP = g++ -O2 -Wall

OBJS = GammaTest.o \
	GammaTest_Estimator.o \
	kdtree_static.o \
	kdtree_2.o \
	Arguments.o \
	Filter.o \
	MainGT

all: $(OBJS)

GammaTest.o: GammaTest.cpp GammaTest.h
	$(GPP) -c GammaTest.cpp -o GammaTest.o -I.
GammaTest_Estimator.o: GammaTest_Estimator.cpp GammaTest_Estimator.h
	$(GPP) -c GammaTest_Estimator.cpp -o GammaTest_Estimator.o -I.
kdtree_static.o: kdtree_static.cpp kdtree_static.h
	$(GPP) -c kdtree_static.cpp -o kdtree_static.o -I.
kdtree_2.o: kdtree_2.cpp kdtree_2.h kdtree_static.h
	$(GPP) -c kdtree_2.cpp -o kdtree_2.o -I.
Arguments.o: Arguments.cpp Arguments.h
	$(GPP) -c Arguments.cpp -o Arguments.o -I.
Filter.o: Filter.cpp Filter.h
	$(GPP) -c Filter.cpp -o Filter.o -I.
MainGT: MainGT.cpp GammaTest.h GammaTest.cpp GammaTest.o Arguments.o Filter.o kdtree_2.o GammaTest_Estimator.h GammaTest_Estimator.cpp GammaTest_Estimator.o
	$(GPP) MainGT.cpp -o MainGT kdtree_static.o kdtree_2.o GammaTest.o GammaTest_Estimator.o Arguments.o Filter.o -I.
clean:
	rm *.o
