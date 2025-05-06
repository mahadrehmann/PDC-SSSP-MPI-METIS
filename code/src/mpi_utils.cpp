// code/src/mpi_utils.cpp

#include "mpi_utils.hpp"
#include <mpi.h>
#include <vector>

/// Broadcast a batch of {insert/delete} changes (Delta) from rank 0 to all
void broadcastChanges(std::vector<Change>& D, int rank) {
    // 1) tell everyone how many changes
    int n = (int)D.size();
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 2) resize on non‑root
    if (rank != 0) D.resize(n);

    // 3) pack into int buffer: [ is_insert, u, v,  is_insert, u, v, ... ]
    std::vector<int> buf(n * 3);
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            buf[3*i + 0] = D[i].is_insert ? 1 : 0;
            buf[3*i + 1] = D[i].u;
            buf[3*i + 2] = D[i].v;
        }
    }

    // 4) broadcast the buffer
    MPI_Bcast(buf.data(), n * 3, MPI_INT, 0, MPI_COMM_WORLD);

    // 5) unpack on non‑root
    if (rank != 0) {
        for (int i = 0; i < n; i++) {
            D[i].is_insert = (buf[3*i+0] == 1);
            D[i].u         = buf[3*i+1];
            D[i].v         = buf[3*i+2];
        }
    }
}

/// Very simple “halo” exchange: here just a barrier
void exchangeBoundaryDistances(LocalSSSPState & /*st*/,
                               const Graph &      /*localG*/,
                               int                /*rank*/,
                               int                /*size*/)
{
    MPI_Barrier(MPI_COMM_WORLD);
}

/// Global convergence: everyone’s localDone must be true (1)
bool checkGlobalConvergence(bool localDone) {
    int in  = localDone ? 1 : 0;
    int out = 0;
    MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return (out == 1);
}
