CC = gcc
CFLAGS = -g -c -Wall -pedantic -Dlinux --std=c99 -D_XOPEN_SOURCE -D_BSD_SOURCE
LFLAGS = -lcfitsio

RSRC = reduce.c framedata.c reduction.c
ROBJ = $(RSRC:.c=.o)

FSRC = frameadd.c framedata.c reduction.c
FOBJ = $(FSRC:.c=.o)

all: reduction frametool
reduction: $(ROBJ)
	$(CC) $(LFLAGS) $(ROBJ) -o $@

frametool: $(FOBJ)
	$(CC) $(LFLAGS) $(FOBJ) -o $@

clean:
	-rm $(OBJ) reduction frametool reduction.o framedata.o reduce.o frameadd.o

.SUFFIXES: .c
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
