LFLAGS=-L/usr/local/Cellar/sfml/2.3/lib -L/usr/local/lib -lsfml-graphics -lsfml-audio -lsfml-window -lsfml-system -lfftw3 -lm
CFLAGS=-Wall -g -O3 -I/usr/local/Cellar/sfml/2.3/include/ -I/usr/local/include/

all: soundwave

soundwave:	soundwave.cpp
	g++ $(CFLAGS) $(LFLAGS) $^ -o $@

clean:
	rm -r soundwave soundwave.dSYM
