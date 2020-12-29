EXEC=a
OPTI=-g
PKGS=
VOX=vox

CC?=gcc
#CFLAGS = -Wall $(OPTI) `pkg-config $(PKGS) --cflags` --std=c99 -g
CFLAGS=-Wall $(OPTI) -g
#LDFLAGS = `pkg-config $(PKGS) --libs` -lGL -lm
LDFLAGS=-lGL -lm -lSDL2 -lscecore -lsceutils `pkg-config --libs glew`
LDFLAGS_VOX=-lm -lscecore -lsceutils

all: $(EXEC) $(VOX)

$(EXEC): manidc.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(VOX): voxelize.o
	$(CC) $^ -o $@ $(LDFLAGS_VOX)

manidc.o: main.c Makefile
	$(CC) $(CFLAGS) -c $< -o $@

voxelize.o: voxelize.c Makefile
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o

mrproper: clean
	rm -f $(EXEC)

.PHONY: all clean mrproper
