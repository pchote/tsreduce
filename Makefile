# Copyright 2010, 2011, 2012 Paul Chote
# This file is part of Puoko-nui, which is free software. It is made available
# to you under the terms of version 3 of the GNU General Public License, as
# published by the Free Software Foundation. For more information, see LICENSE.

USE_READLINE := TRUE

CFLAGS = -g -c -Wall -pedantic -Dlinux --std=c99 -D_POSIX_C_SOURCE=200112L -D_BSD_SOURCE -DUSE_READLINE
LFLAGS = -lcfitsio -lcpgplot -lpgplot -lsofa_c -lm -lreadline

# Mac OS X (with gcc, PGPLOT installed via fink)
ifeq ($(shell uname),Darwin)
    LINKER = gfortran
	CC = gcc
	LINKER = gcc
	CFLAGS += -D_DARWIN_C_SOURCE -I/usr/local/include/
    LFLAGS += -L/usr/X11R6/lib -lX11 -Wl,-framework -Wl,Foundation -lpng -L/usr/local/lib

    # Requires newer libreadline than osx provides by default - installed with hombrew
    LFLAGS += -L/usr/local/opt/readline/lib
    CFLAGS += -I/usr/local/opt/readline/include
else
	CC = gcc
	LINKER = gcc
endif

SRC = tsreduce.c datafile.c framedata.c helpers.c aperture.c fit.c dft_analysis.c reduction.c plots.c random.c hashmap.c
OBJ = $(SRC:.c=.o)

ifeq ($(MSYSTEM),MINGW32)
    LINKER = gfortran
    CFLAGS += -I/usr/local/include
    LFLAGS += -L/usr/local/lib -lgdi32 -lregex
endif

tsreduce: $(OBJ)
	$(LINKER) -o $@ $(OBJ) $(LFLAGS)

clean:
	-rm $(OBJ) tsreduce

.SUFFIXES: .c
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
