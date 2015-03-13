all: fftp

fftp: fftp.c window.c fftlib.c fftlib.h window.h
	cc -ggdb -o fftp fftp.c fftlib.c window.c -lm

install: fftp
	cp fftp /usr/local/bin
	cp fftp.1 /usr/local/man/man1
