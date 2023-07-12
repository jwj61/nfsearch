CC = gcc

INFILEPATH=\"/home/jj/Programs/src-C/6-pari/Inputs\"

CFLAGS =  -O4 -Wall -I/usr/local/include/pari -L/usr/local/lib \
  -DDEG=6 -DINFILES=$(INFILEPATH)
#  -DLONG_IS_32BIT  
#  -DULONG_NOT_DEFINED
#CFLAGS =  -g -Wall -I/usr/local/include/pari -L/usr/local/lib \

FILES = c-generics.o c-degree_dep.o wild-file.o inc_primes.o \
        t6.o t51.o t42.o t411.o t33.o t321.o t3111.o t222.o t2211.o t21111.o \
#        c-registry.o

INCS = c-target.h c-doout.h

all:    c6 c6-stand

c6:	$(FILES) c6-main.o
	$(CC) $(CFLAGS) -o c6 $(FILES) c6-main.o -lpari -lm

c6-stand:	$(FILES) c6-main-stand.o
	$(CC) $(CFLAGS) -o c6-stand $(FILES) c6-main-stand.o -lpari -lm

test-new-funcs:	test-new-funcs.o c-generics.o c-degree_dep.o wild-file.o
	$(CC) -g $(CFLAGS) -o test-new-funcs t6.o t51.o t42.o t411.o t33.o t321.o t3111.o t222.o t2211.o t21111.o wild-file.o c-generics.o c-degree_dep.o inc_primes.o test-new-funcs.o -lpari -lm

c-target.h:

c-doout.h:

