CC = gcc

INFILEPATH=\"/home/jj/Programs/src-C/5-pari/Inputs\"

CFLAGS =  -g -Wall -I/usr/local/include/pari -L/usr/local/lib \
  -DDEG=5 -DINFILES=$(INFILEPATH)
#  -DLONG_IS_32BIT -DDEG=5 -DINFILES=$(INFILEPATH)
#CFLAGS =  -O4 -Wall -I/ud/jj/lib/pari/include -L/ud/jj/lib/pari/lib \
#  -DULONG_NOT_DEFINED

FILES = c-generics.o c-degree_dep.o c5-main.o wild-file.o inc_primes.o \
        t5.o t41.o t32.o t311.o t221.o t2111.o 
#        c-registry.o

FILES2 = c-generics.o c-degree_dep.o c5-main-stand.o wild-file.o inc_primes.o \
        t5.o t41.o t32.o t311.o t221.o t2111.o 

FILES3 = slop-generics.o c-degree_dep.o c5-main-slop.o wild-file.o inc_primes.o \
        t5.o t41.o t32.o t311.o t221.o t2111.o 

INCS = c-target.h c-doout.h

all:    c5 c5-stand c5slop

c5:	$(FILES)
	$(CC) $(CFLAGS) -o c5 $(FILES) -lpari -lm

c5-stand:	$(FILES2)
	$(CC) $(CFLAGS) -o c5-stand $(FILES2) -lpari -lm

c5slop:	$(FILES3)
	$(CC) $(CFLAGS) -o c5slop $(FILES3) -lpari -lm

test-new-funcs:	test-new-funcs.o
	$(CC) $(CFLAGS) -o test-new-funcs t5.o t41.o t32.o t311.o t221.o t2111.o \
    inc_primes.o test-new-funcs.o  -lpari -lm

sort-file: sort-file.o
	$(CC) -o sort-file sort-file.o

c-target.h:

c-doout.h:

