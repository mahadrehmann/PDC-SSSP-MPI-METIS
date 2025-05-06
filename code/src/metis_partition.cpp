#include "metis_partition.hpp"
#include <metis.h>
#include <iostream>

std::vector<int> partitionGraph(const Graph& graph, int numParts) {
    std::vector<int> part(graph.numVertices);

    // Log how many parts we're trying to create
    std::cout << "partitionGraph(): k = " << numParts << "\n";

    if (numParts <= 1) {
        std::fill(part.begin(), part.end(), 0);
        return part;
    }

    // Prepare METIS input types
    idx_t nv = graph.numVertices;
    idx_t ncon = 1;
    idx_t nparts = numParts;
    idx_t objval;

    std::vector<idx_t> xadj(graph.xadj.begin(), graph.xadj.end());
    std::vector<idx_t> adjncy(graph.adjncy.begin(), graph.adjncy.end());
    std::vector<idx_t> metis_part(nv);  // result

    int status = METIS_PartGraphKway(&nv, &ncon,
                                     xadj.data(), adjncy.data(),
                                     nullptr, nullptr, nullptr,
                                     &nparts, nullptr, nullptr, nullptr,
                                     &objval, metis_part.data());

    if (status != METIS_OK) {
        std::cerr << "METIS_PartGraphKway failed with status = " << status << "\n";
        std::fill(part.begin(), part.end(), 0);  // fallback
    } else {
        for (int i = 0; i < graph.numVertices; ++i)
            part[i] = static_cast<int>(metis_part[i]);
    }

    return part;
}
