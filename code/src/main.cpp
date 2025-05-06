// code/src/main.cpp
#include "graph.hpp"
#include "dynamic_update.hpp"
#include "metis_partition.hpp"
#include "mpi_utils.hpp"
#include "sssp_local.hpp"

#include <mpi.h>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>

static const int INF = std::numeric_limits<int>::max();

/// Load the batch of insert/delete changes from file
std::vector<Change> loadChanges() {
    const std::string path = "../data/changes.txt";
    std::vector<Change> D;
    std::ifstream in(path);
    if (!in) {
        std::cerr << "ERROR: cannot open " << path << "\n";
        return D;
    }
    while (true) {
        std::string op;
        int u, v;
        in >> op >> u >> v;
        if (!in) break;
        D.push_back({ op == "insert", u, v });
    }
    return D;
}

/// Broadcast the full CSR from rank 0 to all ranks
void broadcastGraph(Graph &g, int rank) {
    MPI_Bcast(&g.numVertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&g.numEdges,    1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        g.xadj   .resize(g.numVertices + 1);
        g.adjncy .resize(2 * g.numEdges);
        g.weights.resize(2 * g.numEdges);
    }
    MPI_Bcast(g.xadj   .data(), g.numVertices + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g.adjncy .data(), 2 * g.numEdges,   MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g.weights.data(), 2 * g.numEdges,   MPI_INT, 0, MPI_COMM_WORLD);
}

/// Broadcast METIS partition from rank 0
void broadcastPartition(std::vector<int> &part, int rank, int V) {
    if (rank != 0) part.resize(V);
    MPI_Bcast(part.data(), V, MPI_INT, 0, MPI_COMM_WORLD);
}

/// Extract each rank’s subgraph (1-ring halo), returning (subgraph, global→local)
std::pair<Graph,std::vector<int>>
extractLocalSubgraph(const Graph &full,
                     const std::vector<int> &part,
                     int myRank) {
    int V = full.numVertices;
    std::vector<bool> keep(V,false);
    // mark owned vertices
    for(int v=0; v<V; v++)
        if(part[v]==myRank) keep[v]=true;
    // add 1-ring neighbors
    for(int u=0; u<V; u++){
        if(!keep[u]) continue;
        for(int ei=full.xadj[u]; ei<full.xadj[u+1]; ei++){
            keep[ full.adjncy[ei] ] = true;
        }
    }
    // build global→local map
    std::vector<int> g2l(V,-1), local2g;
    local2g.reserve(V);
    for(int v=0; v<V; v++){
        if(keep[v]){
            g2l[v] = (int)local2g.size();
            local2g.push_back(v);
        }
    }
    int localV = (int)local2g.size();
    // count edges
    long long ec=0;
    for(int i=0;i<localV;i++){
        int u = local2g[i];
        for(int ei=full.xadj[u]; ei<full.xadj[u+1]; ei++){
            if(keep[ full.adjncy[ei] ]) ec++;
        }
    }
    // build local graph CSR
    Graph sub;
    sub.numVertices = localV;
    sub.numEdges    = (int)(ec/2);
    sub.xadj.assign(localV+1,0);
    sub.adjncy.resize(ec);
    sub.weights.resize(ec);
    // fill xadj
    for(int i=0;i<localV;i++){
        int u = local2g[i];
        for(int ei=full.xadj[u]; ei<full.xadj[u+1]; ei++){
            if(keep[ full.adjncy[ei] ]) sub.xadj[i+1]++;
        }
    }
    for(int i=1;i<=localV;i++)
        sub.xadj[i] += sub.xadj[i-1];
    // fill edges
    std::vector<int> cursor = sub.xadj;
    for(int i=0;i<localV;i++){
        int u = local2g[i];
        for(int ei=full.xadj[u]; ei<full.xadj[u+1]; ei++){
            int v = full.adjncy[ei];
            if(!keep[v]) continue;
            int pos = cursor[i]++;
            sub.adjncy[pos]  = g2l[v];
            sub.weights[pos] = full.weights[ei];
        }
    }
    return {sub, std::move(g2l)};
}

int main(int argc, char**argv) {
    MPI_Init(&argc,&argv);
    double t0 = MPI_Wtime();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // STEP 1: Load & broadcast the full graph
    Graph full;
    if(rank==0){
        if(!full.loadFromFile("../data/Datasetupd.txt")){
            std::cerr<<"ERROR: cannot load graph\n";
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }
    broadcastGraph(full, rank);
    int globalV = full.numVertices;

    // STEP 2: Partition + extract local subgraph
    Graph localG;
    std::vector<int> g2l;
    if(size<=1){
      std::vector<int> part(globalV,0);
      std::tie(localG,g2l) = extractLocalSubgraph(full, part, 0);
    } else {
      std::vector<int> part;
      if(rank==0) part = partitionGraph(full, size);
      broadcastPartition(part, rank, globalV);
      std::tie(localG,g2l) = extractLocalSubgraph(full, part, rank);
    }

    // Build local2g: local idx -> global ID
    std::vector<int> local2g(localG.numVertices);
    for(int g=0; g<globalV; ++g){
      int l = g2l[g];
      if(l>=0) local2g[l] = g;
    }

    if(rank==0){
      std::cout<<"Global: V="<<globalV
               <<" E="<<full.numEdges<<"\n"
               <<"Parts="<<size<<"\n";
    }
    std::cout<<"Rank "<<rank
             <<": localV="<<localG.numVertices
             <<" localE="<<localG.numEdges<<"\n";

    // STEP 3: Init SSSP state
    LocalSSSPState st;
    st.dist  .assign(localG.numVertices, INF);
    st.parent.assign(localG.numVertices, -1);
    if(g2l[0]>=0) st.dist[g2l[0]] = 0;  // source=0

    // STEP 4: Load & broadcast changes
    std::vector<Change> changes;
    if(rank==0) changes = loadChanges();
    broadcastChanges(changes, rank);

    // STEP 5: Incremental update in batches + profiling
    const int batchSize = 1000;  // tune as needed
    int total = (int)changes.size();
    int nBatch = (total + batchSize - 1) / batchSize;

    double tStart = MPI_Wtime();
    for(int bi=0; bi<nBatch; ++bi){
        int a = bi*batchSize;
        int b = std::min(a+batchSize, total);
        std::vector<Change> batch(changes.begin()+a,
                                  changes.begin()+b);
        if(rank==0){
          std::cout<<"Batch "<<(bi+1)<<"/"<<nBatch
                   <<" changes="<<(b-a)<<"\n";
        }
        applyChangesWithLogging(localG, st, batch);
        exchangeBoundaryDistances(st, local2g, globalV);
    }
    double tEnd = MPI_Wtime();

    if(rank==0){
      std::cout<<"SSSP-update completed in "
               <<(tEnd-tStart)<<" s\n";
      std::cout<<"Total run time: "<<(MPI_Wtime()-t0)<<" s\n";
    }

    MPI_Finalize();
    return 0;
}
