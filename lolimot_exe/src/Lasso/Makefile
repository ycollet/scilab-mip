GCC = gcc -O2 -g

all: lasso.o test_lasso

lasso.o: lasso.c lasso.h
	$(GCC) -c lasso.c -o lasso.o -I.

test_lasso: test_lasso.c lasso.o
	$(GCC) test_lasso.c -o test_lasso lasso.o -I. -lm

clean:
	\rm *.o
	\rm *.exe
