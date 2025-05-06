#pragma once

#include <vector>
#include "graph.hpp"
#include "sssp_local.hpp"
#include "dynamic_update.hpp" // for Change struct

void broadcastChanges(std::vector<Change>& changes, int rank);

void exchangeBoundaryDistances(LocalSSSPState& st,
                               const Graph& localG,
                               int rank,
                               int size);
