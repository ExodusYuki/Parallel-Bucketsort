TARGETS=bucketsort

CFLAGS=-Wall -Wextra -Werror -g

all: $(TARGETS) 

bucketsort: bucketsort.c
	mpicc $(CFLAGS) -o $@ $^ -lm 

clean:
	rm -f $(TARGETS) 
