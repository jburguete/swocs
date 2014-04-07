all: bateria rankine-hugoniot

bateria: bateria.c gen_bateria.pl makefile
	gcc bateria.c -lm -Os -Wall -o bateria
	rm i-* o-* s*
	./bateria

rankine-hugoniot: rankine-hugoniot.c makefile
	gcc rankine-hugoniot.c -lm -Os -Wall -o rankine-hugoniot
	./rankine-hugoniot

clean:
	rm i-* o-* s*
