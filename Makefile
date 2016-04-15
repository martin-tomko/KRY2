NAME=kry
MODULES=$(NAME) factor gmp_helper sieve elliptic_curves primes
OBJS=$(addsuffix .o, $(MODULES))
CC=g++
CFLAGS=-pedantic -Wall -W -g -std=c++0x -pthread
LDLIBS=-lgmpxx -lgmp

ARCHIVE=xtomko02.tar.gz
FILES=Makefile *.cc *.h # doc.pdf

all: $(NAME)

$(NAME): $(OBJS)
	$(CC) $(OBJS) $(LDLIBS) -o $(NAME)

%.o: %.cc %.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY:	clean

clean:
	rm -f *.o
	rm -f $(NAME)

pack:
	tar czvf $(ARCHIVE) $(FILES)
