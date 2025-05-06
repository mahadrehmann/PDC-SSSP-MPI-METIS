// dynamic_update.hpp
#pragma once
#include "graph.hpp"
#include "sssp_local.hpp"

struct Change { bool is_insert; int u, v; };

/// Step 1: identify affected vertices (Alg. 2)
void processCE(const Graph &G,
               LocalSSSPState &st,
               const std::vector<Change> &Delta,
               std::vector<bool> &affected_del,
               std::vector<bool> &affected);

/// Step 2: update deleted subtrees & then iteratively fix (Alg. 3)
void updateAffectedSubgraph(const Graph &G,
                            LocalSSSPState &st,
                            std::vector<bool> &affected_del,
                            std::vector<bool> &affected);

/// Top‑level for one batch of changes (Alg. 1)
void applyChanges(const Graph &G,
                  LocalSSSPState &st,
                  const std::vector<Change> &Delta);

void applyChangesWithLogging(const Graph &G,
    LocalSSSPState &st,
    const std::vector<Change> &Delta);