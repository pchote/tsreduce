# Copyright 2010, 2011 Paul Chote
# This file is part of Puoko-nui, which is free software. It is made available
# to you under the terms of version 3 of the GNU General Public License, as
# published by the Free Software Foundation. For more information, see LICENSE.

CC = gcc
LINKER = gcc
CFLAGS = -g -c -Wall -pedantic -Dlinux --std=c99 -D_XOPEN_SOURCE -D_BSD_SOURCE
LFLAGS = -lcfitsio -lxpa  -lcpgplot -lpgplot

# Mac OS X (with gcc, PGPLOT installed via fink)
ifeq ($(shell uname),Darwin)
    LINKER = gfortran
    LFLAGS += -L/usr/X11R6/lib -lX11 -Wl,-framework -Wl,Foundation -lpng
endif

SRC = tsreduce.c framedata.c helpers.c aperture.c fit.c dft_analysis.c reduction.c ts_analysis.c
OBJ = $(SRC:.c=.o)


tsreduce: $(OBJ)
	$(LINKER) -o $@ $(OBJ) $(LFLAGS)

clean:
	-rm $(OBJ) tsreduce

.SUFFIXES: .c
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
