#include "sssp_local.hpp"
#include <omp.h>
#include <limits>

void runLocalSSSP(const Graph& g, LocalSSSPState& st, int src, int maxIter) {
    const int INF = std::numeric_limits<int>::max();
    st.dist.assign(g.numVertices, INF);
    st.parent.assign(g.numVertices, -1);
    st.dist[src] = 0;
    for (int it = 0; it < maxIter; ++it) {
        #pragma omp parallel for schedule(dynamic)
        for (int u = 0; u < g.numVertices; ++u) {
            int du = st.dist[u];
            if (du == INF) continue;
            for (int ei = g.xadj[u]; ei < g.xadj[u+1]; ++ei) {
                int v = g.adjncy[ei], w = g.weights[ei];
                int nd = du + w;
                #pragma omp critical
                if (nd < st.dist[v]) {
                    st.dist[v] = nd;
                    st.parent[v] = u;
                }
            }
        }
    }
}
