#include ../../make.in

EXEC = sscf
SRC = $(wildcard *.c)
OBJS := $(addsuffix .o, $(basename $(SRC)))

CFLAGS = -O3 -Wall
CC = gcc

LIBS = ${BLAS_LIBS} -lgfortran -latlas -llapack ./erd/liberd.a ./oed/liboed.a

all:${EXEC}

$(EXEC): $(OBJS)
	$(CC) ${CFLAGS} ${LDFLAGS} $(OBJS) -o $(EXEC) ${LIBS}

%.o : %.c Makefile
	$(CC) ${CFLAGS} ${INC} -c $< -o $@

clean: 
	rm -f *.o *.s *.d *~ ${EXEC} 
