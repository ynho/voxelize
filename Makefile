EXEC=a
OPTI=-g
PKGS=

CC?=gcc
#CFLAGS = -Wall $(OPTI) `pkg-config $(PKGS) --cflags` --std=c99 -g
CFLAGS=-Wall $(OPTI) -g
#LDFLAGS = `pkg-config $(PKGS) --libs` -lGL -lm
LDFLAGS=-lGL -lm -lSDL2 -lscecore -lsceutils

all: $(EXEC)

$(EXEC): manidc.o
	$(CC) $^ -o $@ $(LDFLAGS)

manidc.o: main.c Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o

mrproper: clean
	rm -f $(EXEC)

.PHONY: all clean mrproper
