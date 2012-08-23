all: bateria retraso

bateria: bateria-2.c gen_bateria.pl makefile
	gcc bateria-2.c -lm -Os -Wall -o bateria
	rm i-* o-* f-* a-*
	./bateria

retraso: retraso.c makefile
	gcc retraso.c -Os -Wall -o retraso
	rm r-*

clean:
	rm i-* o-* f-* a-* r-* script-1* script-2* script-r*
