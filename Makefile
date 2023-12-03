GCC_FLAGS = Wall -Werror -Wextra -c -ggdb

all: main

main: main.o TP3.o
		gcc -Wall -o main main.o TP3.o -lm

.o:
		gcc $(GCC_FLAGS) $*.c

clean:
		rm -f *.o main	
