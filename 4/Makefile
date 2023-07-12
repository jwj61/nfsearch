CC = gcc

INFILEPATH=\"/home/jj/Programs/src-C/4-pari/Inputs\"

CFLAGS =  -g -Wall -I/usr/local/include/pari -L/usr/local/lib \
  -DDEG=4 -DINFILES=$(INFILEPATH)
#CFLAGS =  -O4 -Wall -I/ud/jj/lib/pari/include -L/ud/jj/lib/pari/lib \
#  -DULONG_NOT_DEFINED   
#  -DLONG_IS_64BIT 

FILES = c-generics.o c-degree_dep.o wild-file.o inc_primes.o \
        t4.o t31.o t22.o t211.o \
#        c-registry.o

INCS = c-target.h c-doout.h

all:    c4 c4-stand

c4:	$(FILES) c4-main.o
	$(CC) $(CFLAGS) -o c4 $(FILES) c4-main.o -lpari -lm

c4-stand:	$(FILES) c4-main-stand.o
	$(CC) $(CFLAGS) -o c4-stand $(FILES) c4-main-stand.o -lpari -lm

test-new-funcs:	test-new-funcs.o
	$(CC) $(CFLAGS) -o test-new-funcs t4.o t31.o t22.o t211.o \
    inc_primes.o test-new-funcs.o  -lpari -lm

c-target.h:

c-doout.h:

