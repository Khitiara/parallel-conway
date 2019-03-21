/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/* Alex Sieberer, Chris Nero    ********************************************/
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <clcg4.h>

#include <mpi.h>
#include <pthread.h>

#if BGQ == 1
#include <hwi/include/gqc/A2_inlines.h>
#else
#define GetTimeBase(notused) (g_processor_frequency * MPI_Wtime(notused))
#endif

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD 0
#define MAX_ADDITIONAL_THREAD_COUNT 63

#define CHECK(x, y) (chunk[x][y] & 1)
#define BIRTH(x, y) (chunk[x][y] |= 2)
#define KILL(x, y) (chunk[x][y] &= ~2)
#define COMMIT(x, y) (chunk[x][y] >>= 1)
#define PERSIST(x, y) (chunk[x][y] |= (chunk[x][y] << 1))

#ifndef rowlen
#define rowlen 32768
#endif

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles = 0;
unsigned long long g_end_cycles = 0;
int g_ticks;
double g_threshold;

// useful runtime constants
int mpi_myrank;
int mpi_commsize;
int rows_per_chunk;

pthread_barrier_t barrier;
typedef unsigned char row[rowlen];
row* chunk;

/***************************************************************************/
/* Function Decls **********************************************************/
/***************************************************************************/

void tick(int start, int end);
void commit(int start, int end);
void* do_ticks(void* arg);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

/**
 * Arguments:
 *  1 - number of threads per rank (including this one)
 *  2 - number of ticks
 *  3 - threshold
 *  4 - output file name
 */
int main(int argc, char* argv[])
{
    int threads_per_rank;
    pthread_t threads[MAX_ADDITIONAL_THREAD_COUNT];
    int start_end[MAX_ADDITIONAL_THREAD_COUNT + 1][2];
    // Example MPI startup and using CLCG4 RNG
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);

    // check and gather arguments
    if (argc != 5) {
        puts("Usage: pconway <threads per rank> <ticks> <threshold> <output file>");
        exit(1);
    }
    rows_per_chunk = rowlen / mpi_commsize;
    threads_per_rank = atoi(argv[1]);
    g_ticks = atoi(argv[2]);
    g_threshold = strtod(argv[3], NULL);
    pthread_barrier_init(&barrier, NULL, threads_per_rank);

    // Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();

    // Note, used the mpi_myrank to select which RNG stream to use.
    // You must replace mpi_myrank with the right row being used.
    // This just show you how to call the RNG.
    printf("Rank %d of %d has been started and a first Random Value of %lf\n",
        mpi_myrank, mpi_commsize, GenVal(mpi_myrank));

    MPI_Barrier(MPI_COMM_WORLD);

    // create our personal chunk of the universe
    // chunk[-1] is the ghost row at the start,
    // chunk[rows_per_chunk] is the ghost row at the end
    chunk = calloc(rows_per_chunk + 2, sizeof(row));
    chunk++;
    {
        // create and start the correct number of additional threads
        int i;
        int rows_per_thread = rows_per_chunk / threads_per_rank;
        for (i = 1; i < threads_per_rank; i++) {
            start_end[i][0] = i * rows_per_thread;
            start_end[i][1] = (i + 1) * rows_per_thread;
            pthread_create(&threads[i - 1], NULL, do_ticks, &start_end[i]);
        }
        // start this thread's work as well
        start_end[0][0] = 0;
        start_end[0][1] = rows_per_thread;
        do_ticks(&start_end[0]);
    }
    // END -Perform a barrier and then leave MPI
    MPI_Barrier(MPI_COMM_WORLD);
    pthread_barrier_destroy(&barrier);
    free(chunk - 1);
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/

/**
 * Run a single tick of the simulation. Does not commit.
 */
void tick(int start, int end)
{
    int living_neighbors, r, c, dr, dc, r1, c1;
    for (r = start; r < end; ++r) {
        for (c = 0; c < rowlen; ++c) {

            living_neighbors = 0;
            // Loop over 3x3 section centered on i,j
            for (dr = -1; dr <= 1; ++dr) {
                for (dc = -1; dc <= 1; ++dc) {
                    // Exclude current cell
                    if (dr != 0 || dc != 0) {
                        // Get actual coordinates
                        c1 = (c + dc) % rowlen;
                        r1 = r + dr;
                        // Check the cell
                        if (CHECK(r1, c1)) {
                            ++living_neighbors;
                        }
                    }
                }
            }
            // Cell dies if less than 2 or more than 3 neighbors
            if (living_neighbors < 2 || living_neighbors > 3) {
                KILL(r, c);
            } else if (living_neighbors == 3) { // Cell is born with exactly 3 neighbors
                BIRTH(r, c);
            } else {
                PERSIST(r, c); // Nothing changes otherwise
            }
        }
    }
}

/**
 * Commit a section of the chunk.
 */
void commit(int start, int end)
{
    int r, c;
    for (r = start; r < end; ++r) {
        for (c = 0; c < rowlen; ++c) {
            COMMIT(r, c);
        }
    }
}

/**
 * Run current thread - takes the start and end row bounds (inclusive start, exclusive end)
 */
void* do_ticks(void* arg)
{
    int(*bounds)[2] = arg;
    int start = (*bounds)[0], end = (*bounds)[1];
    int ticks = g_ticks, rows_per_rank = rows_per_chunk;
    int i;
    for (i = 0; i < ticks; i++) {
        // TODO: check for main thread and recv ghost rows from MPI
        tick(start, end);
        pthread_barrier_wait(&barrier);
        commit(start, end);
        if (start == 0) {
            printf("This is the main thread: global bounds=[%ld, %ld]\n",
                (long)mpi_myrank * rows_per_rank + start,
                (long)mpi_myrank * rows_per_rank + end);
        }
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}
