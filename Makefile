# OSX
ifeq "$(shell uname)" "Darwin"
LFLAGS=-L/usr/local/Cellar/sfml/2.4.0/lib -L/usr/local/lib -lsfml-graphics -lsfml-audio -lsfml-window -lsfml-system -lfftw3 -lm
CFLAGS=-Wall -g -O3 -std=c++11 -I/usr/local/Cellar/sfml/2.4.0/include/ -I/usr/local/include/

# Linux
else
SFLAGS=$(shell sdl2-config --cflags) -DGLEXT_PROTOTYPES
GFLAGS=-lGLU -lGL $(shell sdl2-config --libs)
LFLAGS=-L/usr/local/lib -lsfml-graphics -lsfml-audio -lsfml-window -lsfml-system -lfftw3 -lm
CFLAGS=-Wall -g -O3 -std=c++11
endif

all: soundwave soundwave_gl soundwave-third-octave soundfall

soundwave: soundwave.cpp
	g++ $(CFLAGS) $^ -o $@ $(LFLAGS)

soundwave_gl: soundwaveGL.cpp
	g++ $(CFLAGS) $(SFLAGS) $^ -o $@ $(LFLAGS) $(GFLAGS)

soundfall: soundfall.cpp
	g++ $(CFLAGS) -fopenmp $(SFLAGS) $^ -o $@ $(LFLAGS) -lgomp $(GFLAGS)

soundwave-third-octave: soundwaveThird.cpp
	g++ $(CFLAGS) $(SFLAGS) $^ -o $@ $(LFLAGS) $(GFLAGS)

Soundwavelet.o: Soundwave.cu
	nvcc -ccbin g++ -std=c++11 -m64 -gencode arch=compute_52,code=sm_52 -g -G -O3 -c $< -o $@

Soundwavelet: Soundwavelet.o
	nvcc -ccbin g++ -std=c++11 -m64 -gencode arch=compute_52,code=sm_52 -g -G -O3 -o $@ $^ $(LFLAGS)

clean:
	rm -rf soundwave soundwave_gl soundwave-third-octave Soundwavelet soundfall *.o *.dSYM
