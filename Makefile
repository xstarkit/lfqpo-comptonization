CC = gcc
SIM5LIB = .

CFLAGS = -I$(SIM5LIB) -L$(SIM5LIB) -Wall -O3 -fopenmp -w -g -rdynamic
#CFLAGS = -Wall -O3 -fopenmp

LFLAGS = $(CFLAGS) -lm -lgomp


default: lfqpo



%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm -f *.o


lfqpo-src = \
    $(SIM5LIB)/sim5lib.c \
    lfqpo.c \

lfqpo-obj = $(lfqpo-src:.c=.o)

lfqpo: $(lfqpo-obj)
	$(CC) $(lfqpo-obj) -o $@ $(LFLAGS) 




