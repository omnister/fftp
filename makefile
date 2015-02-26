all: fftp
oldall: fftorig fft2 ft fft fftp

fftorig: fftorig.c
	cc -ggdb -o fftorig fftorig.c -lfftw3 -lm

fft: fft.c window.c fftlib.c
	cc -ggdb -o fft fft.c fftlib.c window.c -lfftw3 -lm

fftp: fftp.c window.c fftlib.c
	cc -ggdb -o fftp fftp.c fftlib.c window.c -lm
	#cc -ggdb -o fftp fftp.c fftlib.c window.c -lfftw3 -lm

ft: fftlib.c ft.c window.c
	cc -ggdb -Wall -o ft ft.c fftlib.c window.c -lm

ft1: fftlib.c ft1.c window.c
	cc -ggdb -Wall -o ft1 ft1.c fftlib.c window.c -lm

install: fftp
	cp fftp /usr/local/bin
