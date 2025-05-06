#include "metis_partition.hpp"
#include <metis.h>

std::vector<int> partitionGraph(const Graph& graph, int numParts) {
    idx_t nv = graph.numVertices, ncon = 1, np = numParts, objval;
    std::vector<idx_t> xadj(graph.xadj.begin(), graph.xadj.end()),
                      adjncy(graph.adjncy.begin(), graph.adjncy.end()),
                      part(nv);
    METIS_PartGraphKway(&nv, &ncon,
                        xadj.data(), adjncy.data(),
                        nullptr, nullptr, nullptr,
                        &np, nullptr, nullptr, nullptr,
                        &objval, part.data());
    return std::vector<int>(part.begin(), part.end());
}
