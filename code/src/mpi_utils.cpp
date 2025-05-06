#include "mpi_utils.hpp"
#include <mpi.h>
#include <limits>

void broadcastChanges(std::vector<Change>& D, int rank) {
    int n = (int)D.size();
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) D.resize(n);

    std::vector<int> buf(n * 3);
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            buf[3*i + 0] = D[i].is_insert ? 1 : 0;
            buf[3*i + 1] = D[i].u;
            buf[3*i + 2] = D[i].v;
        }
    }
    MPI_Bcast(buf.data(), n * 3, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        for (int i = 0; i < n; i++) {
            D[i].is_insert = (buf[3*i + 0] == 1);
            D[i].u         = buf[3*i + 1];
            D[i].v         = buf[3*i + 2];
        }
    }
}

void exchangeBoundaryDistances(LocalSSSPState &st,
                               const std::vector<int> &local2g,
                               int globalV) {
    const int INF = std::numeric_limits<int>::max();
    // 1) build a full-size array of INF
    std::vector<int> localDist(globalV, INF);
    for (size_t i = 0; i < local2g.size(); ++i) {
        localDist[ local2g[i] ] = st.dist[i];
    }

    // 2) all-reduce by MIN
    std::vector<int> globalDist(globalV);
    MPI_Allreduce(localDist.data(),
                  globalDist.data(),
                  globalV,
                  MPI_INT,
                  MPI_MIN,
                  MPI_COMM_WORLD);

    // 3) write back to local
    for (size_t i = 0; i < local2g.size(); ++i) {
        st.dist[i] = globalDist[ local2g[i] ];
    }
}

bool checkGlobalConvergence(bool localDone) {
    int in  = localDone ? 1 : 0;
    int out = 0;
    MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return (out == 1);
}
