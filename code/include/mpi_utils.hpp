#ifndef MPI_UTILS_HPP
#define MPI_UTILS_HPP

#include <vector>
#include "graph.hpp"
#include "dynamic_update.hpp"   // for struct Change
#include "sssp_local.hpp"       // for struct LocalSSSPState

/// Broadcast a batch of {insert/delete} changes (Delta) from rank 0 to all
void broadcastChanges(std::vector<Change>& D, int rank);

/// Exchange ghost (halo) node distances across all processes.
///  - st.dist is in local indexing [0..local2g.size()).
///  - local2g maps local index → global vertex ID (0..globalV-1).
///  - globalV is the total number of vertices.
void exchangeBoundaryDistances(LocalSSSPState &st,
                               const std::vector<int> &local2g,
                               int globalV);

/// Returns true if *all* processes have localDone==true.
bool checkGlobalConvergence(bool localDone);

#endif // MPI_UTILS_HPP
