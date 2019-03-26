/**
 * A small parallel program that reduces the universe data output by pconway
 * into 32x32 cell blocks. This program runs in approx. 12 seconds on a laptop
 * with 8 cores, and presumably faster on the CSCI server.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#define rowlen 32768
#define blocksize 32
#define blocks_per_row (rowlen / blocksize)

unsigned char blockdata[blocksize];
int (*alive_counts)[blocks_per_row];

int main(int argc, char* argv[])
{
    char* in_file_name;
    char* out_file_name;
    int rank, commsize;
    MPI_File in_file, out_file;

    // init MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    // parse arguments
    if (argc != 2) {
        puts("Usage: cv <file>");
        exit(0);
    }
    in_file_name = argv[1];
    {
        int len = strlen(in_file_name);
        out_file_name = malloc(len + sizeof(".out"));
        memcpy(out_file_name, in_file_name, len);
        strcat(out_file_name, ".out");
    }

    // open file handles
    MPI_File_open(MPI_COMM_WORLD, in_file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &in_file);
    MPI_File_open(MPI_COMM_WORLD, out_file_name, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &out_file);

    // compress output into heatmap size
    {
        int r, c;
        // get the number of cells we are handling, and the
        // offset in the file where we start
        int row_count = blocks_per_row / commsize;
        int row_offset = rank * row_count;
        alive_counts = calloc(row_count, sizeof(*alive_counts));
        // loop over all of out blocks
        for (c = 0; c < blocks_per_row; ++c) {
            for (r = 0; r < row_count; ++r) {
                int r1, c1;
                int alive_count = 0;
                // sum all alive cells
                for (r1 = 0; r1 < blocksize; ++r1) {
                    // read the r1'th row of the block at (r, c)
                    MPI_File_read_at(in_file,
                        (row_offset + r * rowlen + c) * blocksize + r1 * rowlen,
                        blockdata,
                        blocksize,
                        MPI_UNSIGNED_CHAR,
                        MPI_STATUS_IGNORE);
                    for (c1 = 0; c1 < blocksize; ++c1) {
                        alive_count += blockdata[c1];
                    }
                }
                alive_counts[r][c] = alive_count;
            }
        }
        // write out data
        MPI_File_write_at(out_file, rank * blocks_per_row * row_count, alive_counts, blocks_per_row * row_count, MPI_INT, MPI_STATUS_IGNORE);
    }

    MPI_File_close(&out_file);
    MPI_File_close(&in_file);
    free(alive_counts);
    free(out_file_name);
    MPI_Finalize();
    return 0;
}
