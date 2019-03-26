BGQ?=0
MPICH_CC?=$(CC)
all: clcg4.h clcg4.c main.c
	$(CC) $(CFLAGS) -I. -Wall -O3 -c clcg4.c -o clcg4.o
	mpicc $(CFLAGS) -I. -Wall -O3 main.c clcg4.o -o pconway -lpthread -DBGQ=$(BGQ)

heat: cv.c
	mpicc $(CFLAGS) -I. -Wall -O3 cv.c -o cv
