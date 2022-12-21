CC=gcc
CFLAGS=-I.

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

main: main.o affine.o jacobien.o
	$(CC) -o main main.o affine.o jacobien.o -lgmp

clean:
	rm -f *.o main

run : main
	./$<