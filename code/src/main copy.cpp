// code/src/main.cpp

#include "graph.hpp"
#include "metis_partition.hpp"
#include "sssp_local.hpp"
#include "mpi_utils.hpp"

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <limits>

// Broadcast full.Graph from rank 0 to all ranks
void broadcastGraph(Graph &g, int rank) {
    // 1) send numVertices, numEdges
    MPI_Bcast(&g.numVertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g.numEdges,    1, MPI_INT, 0, MPI_COMM_WORLD);

    // 2) ensure xadj, adjncy, weights are sized on non‑root
    if (rank != 0) {
        g.xadj.resize(g.numVertices + 1);
        g.adjncy.resize(2 * g.numEdges);
        g.weights.resize(2 * g.numEdges);
    }

    // 3) broadcast arrays
    MPI_Bcast(g.xadj.data(),    g.numVertices + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g.adjncy.data(),  2 * g.numEdges,   MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g.weights.data(), 2 * g.numEdges,   MPI_INT, 0, MPI_COMM_WORLD);
}

// Broadcast the METIS partition array
void broadcastPartition(std::vector<int> &part, int rank, int V) {
    if (rank != 0) {
        part.resize(V);
    }
    MPI_Bcast(part.data(), V, MPI_INT, 0, MPI_COMM_WORLD);
}

// Extract each rank’s subgraph (including halo neighbors)
Graph extractLocalSubgraph(const Graph &full, 
                           const std::vector<int> &part, 
                           int myRank) {
    int V = full.numVertices;
    std::vector<bool> keep(V, false);

    // 1) mark all vertices owned by me
    for (int v = 0; v < V; ++v)
        if (part[v] == myRank)
            keep[v] = true;

    // 2) mark their neighbors too (halo)
    for (int u = 0; u < V; ++u) {
        if (!keep[u]) continue;
        for (int ei = full.xadj[u]; ei < full.xadj[u+1]; ++ei) {
            keep[ full.adjncy[ei] ] = true;
        }
    }

    // 3) build mapping global→local
    std::vector<int> g2l(V, -1), local2g;
    local2g.reserve(V);
    for (int v = 0; v < V; ++v) {
        if (keep[v]) {
            g2l[v] = (int)local2g.size();
            local2g.push_back(v);
        }
    }

    int localV = (int)local2g.size();
    int localE = 0;
    // count edges in local subgraph
    for (int i = 0; i < localV; ++i) {
        int u = local2g[i];
        for (int ei = full.xadj[u]; ei < full.xadj[u+1]; ++ei) {
            int v = full.adjncy[ei];
            if (keep[v])
                localE++;
        }
    }

    Graph sub;
    sub.numVertices = localV;
    sub.numEdges    = localE / 2;    // each undirected edge counted twice
    sub.xadj.assign(localV+1, 0);
    sub.adjncy.resize(localE);
    sub.weights.resize(localE);

    // fill xadj
    for (int i = 0; i < localV; ++i) {
        int u = local2g[i];
        for (int ei = full.xadj[u]; ei < full.xadj[u+1]; ++ei) {
            int v = full.adjncy[ei];
            if (keep[v]) sub.xadj[i+1]++;
        }
    }
    // prefix‑sum
    for (int i = 1; i <= localV; ++i)
        sub.xadj[i] += sub.xadj[i-1];

    // fill adjncy & weights via cursor
    std::vector<int> cursor = sub.xadj;
    for (int i = 0; i < localV; ++i) {
        int u = local2g[i];
        for (int ei = full.xadj[u]; ei < full.xadj[u+1]; ++ei) {
            int v = full.adjncy[ei];
            if (!keep[v]) continue;
            int idx = cursor[i]++;
            sub.adjncy[idx] = g2l[v];
            sub.weights[idx]= full.weights[ei];
        }
    }

    return sub;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    double t0 = MPI_Wtime();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 1) Load or receive full graph
    Graph full;
    if (rank == 0) {
        if (!full.loadFromFile("../data/soc-LiveJournal1.txt")) {
            std::cerr<<"ERROR: failed to load\n";
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }
    broadcastGraph(full, rank);

    // 2) Partition on rank 0 and broadcast
    std::vector<int> part;
    if (size <= 1) {
        part.assign(full.numVertices, 0);
    } else if (rank == 0) {
        part = partitionGraph(full, size);
    }
    broadcastPartition(part, rank, full.numVertices);

    // 3) Extract each rank’s subgraph
    Graph localGraph = extractLocalSubgraph(full, part, rank);

    if (rank==0) {
        std::cout<<"Full loaded: V="<<full.numVertices
                 <<" E="<<full.numEdges
                 <<"\nPartitioned into "<<size<<" parts\n";
    }
    std::cout<<"Rank "<<rank
             <<": localV="<<localGraph.numVertices
             <<" localE="<<localGraph.numEdges
             <<"\n";

    // 4) Init SSSP state on localGraph
    LocalSSSPState state;
    state.dist.assign(localGraph.numVertices, std::numeric_limits<int>::max());
    state.parent.assign(localGraph.numVertices, -1);
    const int sourceGlobal = 0;
    // map global source→local ID if it belongs here
    if (part[sourceGlobal] == rank) {
        // find local index of global 0
        for (int i = 0; i < localGraph.numVertices; ++i)
            if (localGraph.adjncy.size() && localGraph.xadj[i] < localGraph.xadj[i+1] 
                && localGraph.adjncy[ localGraph.xadj[i] ] == 0) {
                state.dist[i] = 0;
                state.parent[i]= -1;
                break;
            }
    }

    // 5) two-phase update
    bool localDone;
    do {
        localDone = true;
        runLocalSSSP(localGraph, state, 0);
        exchangeBoundaryDistances(state, localGraph, rank, size);
    } while (!checkGlobalConvergence(localDone));

    if (rank==0) std::cout<<"SSSP complete\n";

    double t1 = MPI_Wtime();
    if (rank==0) 
      std::cout<<"Total time: "<<(t1-t0)<<" s\n";

    MPI_Finalize();
    return 0;
}
