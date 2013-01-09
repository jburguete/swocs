all: bateria

bateria: bateria.c gen_bateria.pl makefile
	gcc bateria.c -lm -Os -Wall -o bateria
	rm i-* o-* s*
	./bateria

clean:
	rm i-* o-* s*
