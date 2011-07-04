# Copyright 2010, 2011 Paul Chote
# This file is part of Puoko-nui, which is free software. It is made available
# to you under the terms of version 3 of the GNU General Public License, as
# published by the Free Software Foundation. For more information, see LICENSE.

CC = gcc
CFLAGS = -g -c -Wall -pedantic -Dlinux --std=c99 -D_XOPEN_SOURCE -D_BSD_SOURCE
LFLAGS = -lcfitsio

SRC = tsreduce.c framedata.c helpers.c aperture.c
OBJ = $(SRC:.c=.o)


tsreduce: $(OBJ)
	$(CC) $(LFLAGS) $(OBJ) -o $@

clean:
	-rm $(OBJ) tsreduce

.SUFFIXES: .c
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
