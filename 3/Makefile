CC = gcc

INFILEPATH=\"/home/jj/Programs/src-C/3/Inputs\"

CFLAGS =  -g -Wall -I/usr/local/include/pari -L/usr/local/lib \
  -DDEG=3 -DINFILES=$(INFILEPATH)
#CFLAGS =  -O4 -Wall -I/ud/jj/lib/pari/include -L/ud/jj/lib/pari/lib \
#  -DULONG_NOT_DEFINED

FILES = c-generics.o c-degree_dep.o c3-main.o wild-file.o inc_primes.o \
        t3.o t21.o \
#        c-registry.o

INCS = c-target.h c-doout.h

all:    c3

c3:	$(FILES)
	$(CC) $(CFLAGS) -o c3 $(FILES) -lpari -lm

test-new-funcs:	test-new-funcs.o
	$(CC) $(CFLAGS) -o test-new-funcs t3.o t21.o \
    inc_primes.o test-new-funcs.o  -lpari -lm

c-target.h:

c-doout.h:

