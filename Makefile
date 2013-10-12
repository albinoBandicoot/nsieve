CC=gcc
CFLAGS= -g -O3 -pedantic -Wall -std=c99 -I./
USE_ASM=0
MATROW_ASM_FILE=
ifneq ($(USE_ASM),0)
CFLAGS+= -DUSE_ASM
endif
ifeq ($(USE_ASM),32)
	MATROW_ASM_FILE=src/matrow_ops_32.s
else ifeq ($(USE_ASM),64)
	MATROW_ASM_FILE=src/matrow_ops_64.s
endif

vpath %.h src/
vpath %.c src/
vpath %.o build/

OBJECTS= build/*.o
HEADERS= src/*.h

all: nsieve bin/numgen bin/prcheck bin/tdiv bin/rho

bin/prcheck: prcheck.c
	$(CC) $(CFLAGS) -o bin/prcheck src/prcheck.c -lgmp

bin/numgen: numgen.c
	$(CC) $(CFLAGS) -o bin/numgen src/numgen.c -lgmp

bin/tdiv: tdiv.c
	$(CC) $(CFLAGS) -o bin/tdiv src/tdiv.c -lgmp

bin/rho: rho.o
	$(CC) $(CFLAGS) -o bin/rho src/rho.c build/rho.o -lgmp

nsieve: poly.o sieve.o common.o filter.o nsieve.o matrix.o rho.o
ifneq ($(USE_ASM),0)
	gcc -c -g $(MATROW_ASM_FILE) -o build/matrow_ops.o
endif
	ar rc build/libnsieve.a ${OBJECTS} 
	$(CC) $(CFLAGS) -o bin/nsieve src/nsieve.c -Lbuild/ -lnsieve -lgmp -lm

poly.o:	poly.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o build/poly.o src/poly.c 
sieve.o: sieve.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o build/sieve.o src/sieve.c 
common.o: common.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o build/common.o src/common.c
filter.o: filter.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o build/filter.o src/filter.c 
nsieve.o: nsieve.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o build/nsieve.o src/nsieve.c 
matrix.o: matrix.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o build/matrix.o src/matrix.c
rho.o: rhofuncs.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o build/rho.o src/rhofuncs.c

clean:
	rm -f -R build/* bin/*
