BGQ?=0
all: clcg4.h clcg4.c main.c
	gcc $(CFLAGS) -I. -Wall -O3 -c clcg4.c -o clcg4.o
	mpicc $(CFLAGS) -I. -Wall -O3 main.c clcg4.o -o pconway -lpthread -DBGQ=$(BGQ)
