CFLAGS=-O
LD=cc
LDFLAGS=-lm

all: mdf
	$(LD) $(LDFLAGS) m.o uf.o -o mdf

mdf: m.o uf.o m.h

m.o: m.c
uf.o: uf.c

clean:
	rm -rf *.o mdf
