// dynamic_update.cpp
#include "dynamic_update.hpp"
#include <queue>
#include<iostream>
#include <algorithm>
#include <limits>      // <<— needed for std::numeric_limits

void processCE(const Graph &G,
               LocalSSSPState &st,
               const std::vector<Change> &Delta,
               std::vector<bool> &A_del,
               std::vector<bool> &A)
{
  int V = G.numVertices;
  A_del.assign(V,false);
  A.assign(V,false);

  // handle deletions first
  #pragma omp parallel for
  for(int i=0;i<(int)Delta.size();i++){
    if(!Delta[i].is_insert){
      int u=Delta[i].u, v=Delta[i].v;
      // if edge was in the current tree
      if(st.parent[v]==u || st.parent[u]==v){
        int y = (st.dist[u]>st.dist[v]?u:v);
        st.dist[y] = std::numeric_limits<int>::max();
        A_del[y]=A[y]=true;
      }
    }
  }
  // handle insertions
  #pragma omp parallel for
  for(int i=0;i<(int)Delta.size();i++){
    if(Delta[i].is_insert){
      int u=Delta[i].u, v=Delta[i].v, w=1; // assume weight=1
      // pick x,y so y is the farther
      int x = (st.dist[u]>st.dist[v]?u:v);
      int y = (x==u?v:u);
      if(st.dist[y] > st.dist[x] + w){
        st.dist[y]   = st.dist[x]+w;
        st.parent[y] = x;
        A[y]=true;
      }
    }
  }
}

void updateAffectedSubgraph(const Graph &G,
                            LocalSSSPState &st,
                            std::vector<bool> &A_del,
                            std::vector<bool> &A)
{
  int V = G.numVertices;
  // Part 1: propagate all-infinity down deleted subtrees
  std::queue<int> Q;
  for(int v=0;v<V;v++){
    if(A_del[v]) Q.push(v);
  }
  while(!Q.empty()){
    int u=Q.front(); Q.pop();
    A_del[u]=false;
    // assume we can iterate children: we scan all v and check parent[v]==u
    for(int v=0;v<V;v++){
      if(st.parent[v]==u){
        st.dist[v] = std::numeric_limits<int>::max();
        A_del[v]=A[v]=true;
        Q.push(v);
      }
    }
  }

  // Part 2: iterative relaxations until no more updates (Alg. 3)
  bool again=true;
  while(again){
    again=false;
    #pragma omp parallel for reduction(|:again)
    for(int u=0;u<V;u++){
      if(A[u]){
        A[u]=false;
        // scan neighbors
        for(int ei=G.xadj[u]; ei<G.xadj[u+1]; ei++){
          int v=G.adjncy[ei], w=G.weights[ei];
          if(st.dist[v] > st.dist[u]+w){
            st.dist[v]=st.dist[u]+w;
            st.parent[v]=u;
            again = A[v]=true;
          }
          else if(st.dist[u] > st.dist[v]+w){
            st.dist[u]=st.dist[v]+w;
            st.parent[u]=v;
            again = A[u]=true;
          }
        }
      }
    }
  }
}

void applyChanges(const Graph &G,
                  LocalSSSPState &st,
                  const std::vector<Change> &Delta)
{
  std::vector<bool> A_del, A;
  processCE(G,st,Delta,A_del,A);
  updateAffectedSubgraph(G,st,A_del,A);
}

void applyChangesWithLogging(const Graph &G,
  LocalSSSPState &st,
  const std::vector<Change> &Delta)
{
std::cout << "Applying changes:\n";
for (const auto &change : Delta) {
if (change.is_insert) {
std::cout << "  [INSERT]  Edge (" << change.u << ", " << change.v << ")\n";
} else {
std::cout << "  [DELETE]  Edge (" << change.u << ", " << change.v << ")\n";
}
}

// Then apply normally
applyChanges(G, st, Delta);
}
