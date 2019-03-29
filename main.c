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
#include <hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase(notused) (g_processor_frequency * MPI_Wtime(notused))
#endif

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD 0
#define MAX_ADDITIONAL_THREAD_COUNT 63
#define MAX_TICKS 256
#define MAX_FILE_NAME 128

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
pthread_mutex_t alive_lock;
pthread_t threads[MAX_ADDITIONAL_THREAD_COUNT];
int start_end[MAX_ADDITIONAL_THREAD_COUNT + 2];
int alive_count[MAX_TICKS];
int final_alive_count[MAX_TICKS];
typedef unsigned char row[rowlen];
row* chunk;

/***************************************************************************/
/* Function Decls **********************************************************/
/***************************************************************************/

void tick(int start, int end, int global_index_offset, double threshold);
void tick_normally(int local_row);
void tick_randomly(int local_row, Gen g);
void commit(int start, int end, int* alive_count);
void* do_ticks(void* arg);
void write_universe(char* fpath);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

/**
 * Arguments:
 *  1 - number of threads per rank (including this one)
 *  2 - number of ticks
 *  3 - threshold
 *  4 - output universe and heatmap (0 or 1)
 *  5 - output file name
 */
int main(int argc, char* argv[])
{
    int threads_per_rank;
    char universe_out_name[MAX_FILE_NAME];
    char alive_stats_name[MAX_FILE_NAME];
    // Example MPI startup and using CLCG4 RNG
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);

    // check and gather arguments
    if (argc != 6) {
        puts("Usage: pconway <threads per rank> <ticks> <threshold> <output type> <output file>");
        exit(1);
    }
    rows_per_chunk = rowlen / mpi_commsize;
    threads_per_rank = atoi(argv[1]);
    g_ticks = atoi(argv[2]);
    g_threshold = strtod(argv[3], NULL);
    strcpy(universe_out_name, argv[5]);
    strcpy(alive_stats_name, argv[5]);
    strcat(universe_out_name, ".bin");
    strcat(alive_stats_name, ".csv");
    pthread_barrier_init(&barrier, NULL, threads_per_rank);
    pthread_mutex_init(&alive_lock, NULL);

    // Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();

    MPI_Barrier(MPI_COMM_WORLD);

    //start timer
    if (mpi_myrank == 0) {
        g_start_cycles = GetTimeBase();
    }
    // create our personal chunk of the universe
    // chunk[-1] is the ghost row at the start,
    // chunk[rows_per_chunk] is the ghost row at the end
    chunk = calloc(rows_per_chunk + 2, sizeof(row));
    // set every slot to ALIVE
    memset(chunk, ALIVE, (rows_per_chunk + 2) * sizeof(row));
    chunk++;
    {
        // create and start the correct number of additional threads
        int i;
        int rows_per_thread = rows_per_chunk / threads_per_rank;
        // initialize thread arguments
        for (i = 0; i < threads_per_rank + 1; i++) {
            start_end[i] = i * rows_per_thread;
        }
        for (i = 1; i < threads_per_rank; i++) {
            pthread_create(&threads[i - 1], NULL, do_ticks, &start_end[i]);
        }
        // start this thread's work as well
        do_ticks(&start_end[0]);
        // collect alive counts
        MPI_Reduce(alive_count, final_alive_count, MAX_TICKS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        // stop timer
        if (mpi_myrank == 0) {
            FILE* alive_stats_file;
            double time;
            int ticks = g_ticks;
            g_end_cycles = GetTimeBase();
            time = (g_end_cycles - g_start_cycles) / g_processor_frequency;
            printf("Computation statistics:\n"
                   "    Compute time (s): %f\n"
                   "       Compute ticks: %d\n"
                   "   Avg time/tick (s): %f\n"
                   "           MPI Ranks: %d\n"
                   "    Threads per rank: %d\n"
                   "       Total threads: %d\n"
                   "       Rows per rank: %d\n"
                   "     Rows per thread: %d\n",
                time,
                ticks,
                time / ticks,
                mpi_commsize,
                threads_per_rank,
                threads_per_rank * mpi_commsize,
                rows_per_chunk,
                rows_per_thread);
            printf("Alive statistics (CSV) written to '%s'\n", alive_stats_name);
            alive_stats_file = fopen(alive_stats_name, "w");
            fputs("Tick,Alive\n", alive_stats_file);
            for (i = 0; i < ticks; ++i) {
                fprintf(alive_stats_file, "%d,%d\n", i, final_alive_count[i]);
            }
            fclose(alive_stats_file);
        }
    }
    pthread_barrier_destroy(&barrier);
    // write universe to disk if specified
    if (atoi(argv[4])) {
        // synchronize the ranks before starting timer
        MPI_Barrier(MPI_COMM_WORLD);
        if (mpi_myrank == 0) {
            g_start_cycles = GetTimeBase();
        }
        write_universe(universe_out_name);
        // write_universe is already synchronized
        if (mpi_myrank == 0) {
            double time;
            g_end_cycles = GetTimeBase();
            time = (g_end_cycles - g_start_cycles) / g_processor_frequency;
            printf("Universe data written to '%s'\n"
                   "     IO run time (s): %f\n",
                universe_out_name, time);
        }
    }
    // END - leave MPI
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
void tick(int start, int end, int global_index_offset, double threshold)
{
    int r;
    for (r = start; r < end; ++r) {
        if (GenVal(global_index_offset + r) < threshold) {
            tick_randomly(r, global_index_offset + r);
        } else {
            tick_normally(r);
        }
    }
}

/**
 * Run a single normal tick of the given row. Does not commit.
 */
void tick_normally(int local_row)
{
    int living_neighbors, c, dr, dc, r1, c1;
    for (c = 0; c < rowlen; ++c) {
        living_neighbors = 0;
        // Loop over 3x3 section centered on i,j
        for (dr = -1; dr <= 1; ++dr) {
            for (dc = -1; dc <= 1; ++dc) {
                // Exclude current cell
                if (dr != 0 || dc != 0) {
                    // Get actual coordinates
                    c1 = (c + dc) % rowlen;
                    r1 = local_row + dr;
                    // Check the cell
                    if (CHECK(r1, c1)) {
                        ++living_neighbors;
                    }
                }
            }
        }
        // Cell dies if less than 2 or more than 3 neighbors
        if (living_neighbors < 2 || living_neighbors > 3) {
            KILL(local_row, c);
        } else if (living_neighbors == 3) { // Cell is born with exactly 3 neighbors
            BIRTH(local_row, c);
        } else {
            PERSIST(local_row, c); // Nothing changes otherwise
        }
    }
}

/**
 * Run a single random tick of the given row. Does not commit.
 */
void tick_randomly(int local_row, Gen g)
{
    int c;
    for (c = 0; c < rowlen; c++) {
        // randomly set to alive or dead
        if (GenVal(g) < 0.5) {
            BIRTH(local_row, c);
        } else {
            KILL(local_row, c);
        }
    }
}

/**
 * Commit a section of the chunk.
 */
void commit(int start, int end, int* alive_cells)
{
    int r, c;
    for (r = start; r < end; ++r) {
        for (c = 0; c < rowlen; ++c) {
            *alive_cells += COMMIT(r, c);
        }
    }
}

/**
 * Run current thread - takes the start and end row bounds (inclusive start, exclusive end)
 */
void* do_ticks(void* arg)
{
    int* bounds = arg;
    int start = bounds[0], end = bounds[1];
    int ticks = g_ticks, rows_per_rank = rows_per_chunk;
    int myrank = mpi_myrank, commsize = mpi_commsize;
    int prev_rank = (myrank - 1 + commsize) % commsize;
    int next_rank = (myrank + 1) % commsize;
    int global_index_offset = myrank * rows_per_rank;
    double threshold = g_threshold;
    MPI_Request send[2];
    MPI_Request recv[2];
    int i;
    for (i = 0; i < ticks; i++) {
        int alive_cells = 0;
        tick(start, end, global_index_offset, threshold);
        // wait for all threads to finish ticking before swapping out ghost rows
        pthread_barrier_wait(&barrier);
        // thread 0 receives ghost rows via MPI
        if (start == 0) {
            MPI_Irecv(&chunk[-1], rowlen, MPI_UNSIGNED_CHAR, prev_rank, i, MPI_COMM_WORLD, &recv[0]);
            MPI_Irecv(&chunk[rows_per_rank], rowlen, MPI_UNSIGNED_CHAR, next_rank, i, MPI_COMM_WORLD, &recv[1]);
        }
        commit(start, end, &alive_cells);
        // put alive cell count into alive counts array
        pthread_mutex_lock(&alive_lock);
        alive_count[i] += alive_cells;
        pthread_mutex_unlock(&alive_lock);
        // wait for all threads to finish commiting before sending boundary rows
        pthread_barrier_wait(&barrier);
        // thread 0 sends boundary rows
        if (start == 0) {
            MPI_Isend(&chunk[0], rowlen, MPI_UNSIGNED_CHAR, prev_rank, i, MPI_COMM_WORLD, &send[0]);
            MPI_Isend(&chunk[rows_per_rank - 1], rowlen, MPI_UNSIGNED_CHAR, next_rank, i, MPI_COMM_WORLD, &send[1]);
            // wait for ghost rows to be received and boundary rows to be sent
            MPI_Waitall(2, recv, MPI_STATUSES_IGNORE);
            MPI_Waitall(2, send, MPI_STATUSES_IGNORE);
        }
        // wait for all threads (read: thread 0/MPI) to finish before starting
        // the next tick
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

void write_universe(char* fpath)
{
    int mpi_rank = mpi_myrank, rows_per_rank = rows_per_chunk;
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, fpath, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_preallocate(fh, rowlen * rowlen);
    for (int i = 0; i < rows_per_rank; ++i) {
        MPI_File_write_at(fh, (mpi_rank * rows_per_rank + i) * rowlen, chunk[i], rowlen, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
    }
    MPI_File_close(&fh);
}
