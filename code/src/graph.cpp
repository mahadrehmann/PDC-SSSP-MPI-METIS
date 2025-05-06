#include "graph.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <tuple>

bool Graph::loadFromFile(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "ERROR: cannot open " << filename << "\n";
        return false;
    }

    // Read all edges *with* weights
    std::vector<std::tuple<int,int,int>> edges;
    edges.reserve(1000000);

    std::string line;
    int maxNode = -1;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        int u, v, w;
        if (!(iss >> u >> v >> w))
            continue;   // malformed line?
        edges.emplace_back(u, v, w);
        maxNode = std::max(maxNode, std::max(u, v));
    }

    // Build CSR for an *undirected*, *weighted* graph
    numVertices = maxNode + 1;
    numEdges    = static_cast<int>(edges.size());

    xadj.assign(numVertices + 1, 0);

    // Count degrees
    for (auto &t : edges) {
        int u, v, w;
        std::tie(u,v,w) = t;
        xadj[u + 1]++;
        xadj[v + 1]++;
    }
    // Prefix‐sum to get offsets
    for (int i = 1; i <= numVertices; ++i) {
        xadj[i] += xadj[i - 1];
    }

    // Allocate adjacency and weight arrays (each undirected edge appears twice)
    adjncy .resize(numEdges * 2);
    weights.resize(numEdges * 2);

    // Fill CSR using a cursor
    std::vector<int> ptr = xadj;
    for (auto &t : edges) {
        int u, v, w;
        std::tie(u,v,w) = t;

        // u -> v
        adjncy[ ptr[u] ] = v;
        weights[ ptr[u] ] = w;
        ptr[u]++;

        // v -> u
        adjncy[ ptr[v] ] = u;
        weights[ ptr[v] ] = w;
        ptr[v]++;
    }

    std::cout << "Loaded weighted graph: V=" << numVertices
              << "   E=" << numEdges
              << "   CSR‐entries=" << adjncy.size()
              << "\n";
    return true;
}
