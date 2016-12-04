TARGETS=bucketsort

CFLAGS=-Wall -Wextra -Werror -g

all: $(TARGETS) 

bucketsort: bucketsort.c
	gcc $(CFLAGS) -o $@ $^

clean:
	rm -f $(TARGETS) 
