GCC = gcc -O2 -Wall -pg -ggdb -march=pentium4 -fno-rtti
GPP = g++ -O2 -Wall -pg -ggdb -march=pentium4 -fno-rtti

#LIB = ../Math.o
#LIB = ../Math2.o
LIB = ../Math3.o -lgsl
LAPACK = ../SmallLapack/liblapack.a

OBJ = genData\
	test_lolimot\
	test_lolimot_step\
	test_svdfit

all: $(OBJ)

genData: genData.cpp
	$(GPP) -O2 -g genData.cpp -o genData
test_lolimot: test_lolimot.cpp ../libLolimot.a
	$(GPP) test_lolimot.cpp -o test_lolimot -I. -I.. ../libLolimot.a -lm $(LAPACK)
test_lolimot_step: test_lolimot_step.cpp ../libLolimot.a
	$(GPP) test_lolimot_step.cpp -o test_lolimot_step -I. -I.. ../libLolimot.a -lm $(LAPACK)
test_svdfit: test_svdfit.cpp ../Math.o
	$(GPP) test_svdfit.cpp -o test_svdfit -I. -I.. $(LIB) -lm 

clean:
	echo "Rien à faire"
