prog := altgroup
prog_objs := altgroup.o

CC := gcc
CFLAGS := -O2 -Wall -g -std=c11
LDFLAGS := -lm -lgmp

.PHONY: all clean

all: $(prog)

$(prog): $(prog_objs)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

altgroup.o: altgroup.c

clean:
	@rm -rf *.o $(prog)
