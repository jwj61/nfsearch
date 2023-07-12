CC = gcc

INFILEPATH=\"/home/jj/Programs/src-C/7-pari/Inputs\"

CFLAGS = -O4 -Wall -I/usr/local/include/pari -L/usr/local/lib \
  -DDEG=7  -DINFILES=$(INFILEPATH)
# -ggdb
#  -DULONG_NOT_DEFINED
# -O4

FILES = c-generics.o c-degree_dep.o wild-file.o inc_primes.o \
        t7.o t61.o t52.o t511.o t43.o t421.o t4111.o t331.o t322.o t3211.o \
        t31111.o t211111.o t22111.o t2221.o

INCS = c-target.h c-doout.h

all:    c7 c7-stand eparc7-stand

c7:	$(FILES) c7-main.o
	$(CC) $(CFLAGS) -o c7 $(FILES) c7-main.o -lpari -lm

c7-stand:	$(FILES) c7-main-stand.o
	$(CC) $(CFLAGS) -o c7-stand $(FILES) c7-main-stand.o -lpari -lm

eparc7-stand:	$(FILES) eparc7-main-stand.o
	$(CC) $(CFLAGS) -o eparc7-stand $(FILES) eparc7-main-stand.o -lpari -lm

test-new-funcs: $(FILES) test-new-funcs.o
	$(CC) $(CFLAGS) -o test-new-funcs test-new-funcs.o $(FILES) -lpari -lm

c-target.h:

c-doout.h:

