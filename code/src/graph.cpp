// code/src/graph.cpp

#include "graph.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

bool Graph::loadFromFile(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "ERROR: cannot open " << filename << "\n";
        return false;
    }

    std::vector<std::pair<int,int>> edges;
    edges.reserve(1000000);

    std::string line;
    int maxNode = -1;

    // Read line by line
    while (std::getline(in, line)) {
        // skip blank or comment lines
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);
        int u, v;
        if (!(iss >> u >> v))
            continue;  // malformed?

        edges.emplace_back(u, v);
        maxNode = std::max(maxNode, std::max(u, v));
    }

    // Build CSR for an undirected, unweighted graph (default weight = 1)
    numVertices = maxNode + 1;
    numEdges    = static_cast<int>(edges.size());

    // xadj[i] will be the start index of neighbors of vertex i
    xadj.assign(numVertices + 1, 0);

    // count degrees
    for (auto &e : edges) {
        xadj[e.first + 1]++;
        xadj[e.second + 1]++;
    }
    // prefix-sum to get starting offsets
    for (int i = 1; i <= numVertices; ++i) {
        xadj[i] += xadj[i - 1];
    }

    // each undirected edge appears twice in CSR
    adjncy .resize(numEdges * 2);
    weights.resize(numEdges * 2, 1);  // set all edge weights = 1

    // temporary cursor array to track insert positions
    std::vector<int> ptr = xadj;
    for (auto &e : edges) {
        int u = e.first, v = e.second;
        adjncy[ ptr[u]++ ] = v;
        adjncy[ ptr[v]++ ] = u;
    }

    std::cout << "Loaded graph: V=" << numVertices
              << "   E=" << numEdges
              << "   CSR-size=" << adjncy.size()
              << "\n";
    return true;
}
