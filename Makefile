CC=clang
CFLAGS= -g -pedantic -Wall -std=c99 -I./

vpath %.h src/
vpath %.c src/
vpath %.o build/

OBJECTS= build/*.o

all: nsieve numgen prcheck

prcheck: prcheck.c
	$(CC) $(CFLAGS) -o bin/prcheck src/prcheck.c -lgmp

numgen: numgen.c
	$(CC) $(CFLAGS) -o bin/numgen src/numgen.c -lgmp

nsieve: poly.o sieve.o common.o filter.o nsieve.o matrix.o
	ar rc build/libnsieve.a ${OBJECTS} 
	$(CC) $(CFLAGS) -o bin/nsieve src/nsieve.c -Lbuild/ -lnsieve -lgmp -lm

poly.o:	poly.c poly.h 
	$(CC) $(CFLAGS) -c -o build/poly.o src/poly.c 
sieve.o: sieve.c sieve.h 
	$(CC) $(CFLAGS) -c -o build/sieve.o src/sieve.c 
common.o: common.c common.h
	$(CC) $(CFLAGS) -c -o build/common.o src/common.c
filter.o: filter.c filter.h 
	$(CC) $(CFLAGS) -c -o build/filter.o src/filter.c 
nsieve.o: nsieve.c nsieve.h 
	$(CC) $(CFLAGS) -c -o build/nsieve.o src/nsieve.c 
matrix.o: matrix.c matrix.h 
	$(CC) $(CFLAGS) -c -o build/matrix.o src/matrix.c

clean:
	rm -f build/*.o build/libnsieve.a bin/nsieve 
