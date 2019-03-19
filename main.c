/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/* Team Names Here              **(*****************************************/
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

#define CHECK(x, y) (chunk[x][y] & 1)
#define BIRTH(x, y) (chunk[x][y] |= 2)
#define KILL(x, y) (chunk[x][y] &= ~2)
#define COMMIT(x, y) (chunk[x][y] >>= 1)
#define PERSIST(x, y) (chunk[x][y] |= (chunk[x][y] << 1))

#define rowlen 32768

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles = 0;
unsigned long long g_end_cycles = 0;

// You define these

typedef unsigned char row[rowlen];
row* chunk;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

void tick(int start, int end);
void commit(int start, int end);

// You define these

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char* argv[])
{
    //    int i = 0;
    int mpi_myrank;
    int mpi_commsize;
    int rows_per_chunk;
    // Example MPI startup and using CLCG4 RNG
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);
    rows_per_chunk = rowlen / mpi_commsize;

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
    chunk = calloc(rows_per_chunk + 2, sizeof(row)) + 1;

    // END -Perform a barrier and then leave MPI
    MPI_Barrier(MPI_COMM_WORLD);
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
    int living_neighbors, x, y, i, j, k, l;
    for (i = 0; i < rowlen; ++i) {
        for (j = start; j < end; ++j) {

            living_neighbors = 0;
            // Loop over 3x3 section centered on i,j
            for (k = -1; k <= 1; ++k) {
                for (l = -1; l <= 1; ++l) {
                    // Exclude current cell
                    if (k != 0 || l != 0) {
                        // Get actual coordinates
                        x = (i + k) % rowlen;
                        y = j + l;
                        // Check the cell
                        if (CHECK(x, y)) {
                            ++living_neighbors;
                        }
                    }
                }
            }
            // Cell dies if less than 2 or more than 3 neighbors
            if (living_neighbors < 2 || living_neighbors > 3) {
                KILL(i, j);
            } else if (living_neighbors == 3) { // Cell is born with exactly 3 neighbors
                BIRTH(i, j);
            } else {
                PERSIST(i, j); // Nothing changes otherwise
            }
        }
    }
}

/**
 * Commit a section of the chunk.
 */
void commit(int start, int end)
{
    int i, j;
    for (i = 0; i < rowlen; ++i) {
        for (j = start; j < end; ++j) {
            COMMIT(i, j);
        }
    }
}
