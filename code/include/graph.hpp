// code/include/graph.hpp
#pragma once
#include <vector>
#include <string>

struct Graph {
    int numVertices = 0, numEdges = 0;
    std::vector<int> xadj, adjncy, weights;

    // Skip “#” lines, parse “u v” per line, weight=1
    bool loadFromFile(const std::string& filename);
};
