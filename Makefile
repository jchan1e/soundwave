# OSX
ifeq "$(shell uname)" "Darwin"
LFLAGS=-L/usr/local/Cellar/sfml/2.4.0/lib -L/usr/local/lib -lsfml-graphics -lsfml-audio -lsfml-window -lsfml-system -lfftw3 -lm
CFLAGS=-Wall -g -O3 -I/usr/local/Cellar/sfml/2.4.0/include/ -I/usr/local/include/

# Otherwise assume Linux
else
LFLAGS=-L/usr/local/lib -lsfml-graphics -lsfml-audio -lsfml-window -lsfml-system -lfftw3 -lm
CFLAGS=-Wall -g -O3 -std=c++11
endif

all: soundwave# Soundwavelet

soundwave:	soundwave.cpp
	g++ $(CFLAGS) $^ -o $@ $(LFLAGS)

Soundwavelet.o: Soundwave.cu
	nvcc -ccbin g++ -std=c++11 -m64 -gencode arch=compute_52,code=sm_52 -g -G -O3 -c $< -o $@

Soundwavelet: Soundwavelet.o
	nvcc -ccbin g++ -std=c++11 -m64 -gencode arch=compute_52,code=sm_52 -g -G -O3 -o $@ $^ $(LFLAGS)

clean:
	rm -r soundwave Soundwavelet *.o *.dSYM
