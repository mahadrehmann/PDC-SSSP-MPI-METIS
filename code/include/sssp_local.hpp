#pragma once
#include "graph.hpp"
#include <vector>

struct LocalSSSPState {
    std::vector<int> dist;
    std::vector<int> parent;
    std::vector<bool> affected;
};

void runLocalSSSP(const Graph& localGraph, LocalSSSPState& state, int source, int maxIters = 10);
