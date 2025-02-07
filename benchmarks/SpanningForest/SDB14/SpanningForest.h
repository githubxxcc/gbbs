// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include "contract_sf.h"
#include "benchmarks/LowDiameterDecomposition/MPX13/LowDiameterDecomposition.h"
#include "ligra/ligra.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/pbbslib/dyn_arr.h"

namespace spanning_forest {

  using edge = std::pair<uintE, uintE>;

  template <class W>
  struct LDD_Edge_F {
    uintE* cluster_ids;
    uintE* parents;

    LDD_Edge_F(uintE* _cluster_ids, uintE* _parents)
        : cluster_ids(_cluster_ids), parents(_parents) {}

    inline bool update(const uintE& s, const uintE& d, const W& wgh) {
      cluster_ids[d] = cluster_ids[s];
      parents[d] = s;
      return true;
    }

    inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
      if (pbbslib::atomic_compare_and_swap(&cluster_ids[d], UINT_E_MAX, cluster_ids[s])) {
        parents[d] = s;
        return true;
      }
      return false;
    }

    inline bool cond(uintE d) { return cluster_ids[d] == UINT_E_MAX; }
  };

  // Returns a pair containing the cluster_ids and parents.
  template <class Graph>
  inline std::pair<sequence<uintE>, sequence<uintE>> LDD_edges(Graph& G,
      double beta, bool permute = true, bool pack = false) {
    using W = typename Graph::weight_type;
    size_t n = G.n;

    sequence<uintE> vertex_perm;
    if (permute) {
      vertex_perm = pbbslib::random_permutation<uintE>(n);
    }
    auto shifts = ldd_utils::generate_shifts(n, beta);
    auto cluster_ids = sequence<uintE>(n);
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { cluster_ids[i] = UINT_E_MAX; });

    auto parents = sequence<uintE>(n, UINT_E_MAX);

    size_t round = 0, num_visited = 0;
    vertexSubset frontier(n);  // Initially empty
    size_t num_added = 0;
    while (num_visited < n) {
      size_t start = shifts[round];
      size_t end = std::min(static_cast<size_t>(shifts[round + 1]), n);
      size_t num_to_add = end - start;
      if (num_to_add > 0) {
        assert((num_added + num_to_add) <= n);
        auto candidates_f = [&](size_t i) {
          if (permute)
            return vertex_perm[num_added + i];
          else
            return static_cast<uintE>(num_added + i);
        };
        auto candidates = pbbslib::make_sequence<uintE>(num_to_add, candidates_f);
        auto pred = [&](uintE v) { return cluster_ids[v] == UINT_E_MAX; };
        auto new_centers = pbbslib::filter(candidates, pred);
        add_to_vsubset(frontier, new_centers.begin(), new_centers.size());
        par_for(0, new_centers.size(), pbbslib::kSequentialForThreshold, [&] (size_t i) {
          uintE v = new_centers[i];
          cluster_ids[v] = v;
          parents[v] = v;
        });
        num_added += num_to_add;
      }

      num_visited += frontier.size();
      if (num_visited >= n) break;

      auto ldd_f = LDD_Edge_F<W>(cluster_ids.begin(), parents.begin());
      vertexSubset next_frontier =
          edgeMap(G, frontier, ldd_f, -1, sparse_blocked);
      frontier.del();
      frontier = next_frontier;

      round++;
    }
    return std::make_pair(cluster_ids, parents);
  }

  // edge_mapping: edge -> edge
  template <class Graph>
  inline pbbslib::dyn_arr<edge> SpanningForest_Impl(Graph& G, double beta,
                                            size_t level, std::function<edge(edge)>& edge_mapping, bool
                                            pack = false, bool permute = false)
  {
    permute |= (level > 0);
    timer ldd_t;
    ldd_t.start();
    auto clusters_and_parents = LDD_edges(G, beta, permute, pack);
    auto clusters = std::move(clusters_and_parents.first);
    auto parents = std::move(clusters_and_parents.second);
    ldd_t.stop();
    debug(ldd_t.reportTotal("ldd time"););

    // Filter out tree edges added this round (ids are in the current level)
    auto delayed_edges = pbbs::delayed_seq<edge>(parents.size(), [&] (size_t i) {
        return std::make_pair(parents[i], i); });
    auto edges = pbbs::filter(delayed_edges, [&] (const edge& e) { return e.first != e.second; });
    // Apply the mapping to map
    par_for(0, edges.size(), [&] (size_t i) {
      auto e_i = edges[i];
      edges[i] = edge_mapping(e_i);
    });

    timer relabel_t;
    relabel_t.start();
    size_t num_clusters = contract_sf::RelabelIds(clusters);
    relabel_t.stop();
    debug(relabel_t.reportTotal("relabel time"););

    timer contract_t;
    contract_t.start();

    // The contraction here also returns a mapping from edge --> edge. This is
    // because edges incident to a single contracted vertex can come from
    // multiple original vertices.
    auto GC_and_new_mapping = contract_sf::contract(G, clusters, num_clusters, edge_mapping);
    contract_t.stop();
    debug(contract_t.reportTotal("contract time"););
    auto GC = GC_and_new_mapping.first;
    auto& new_mapping = GC_and_new_mapping.second; // sparse_table<edge, edge>

    size_t edges_size = edges.size();
    if (GC.m == 0) return pbbslib::dyn_arr<edge>(edges.to_array(), edges_size, edges_size, true);

    auto empty_val = std::make_pair(UINT_E_MAX, UINT_E_MAX);
    std::function<edge(edge)> new_edge_mapping = [&] (edge e) {
      uintE u = std::min(e.first, e.second);
      uintE v = std::max(e.first, e.second);
      auto ret = new_mapping.find(std::make_pair(u,v), empty_val);
      assert(ret != empty_val);
      return ret;
    };

    auto rec_edge_arr = SpanningForest_Impl(GC, beta, level + 1, new_edge_mapping);
    rec_edge_arr.copyIn(edges, edges.size());
    GC.del();
    return rec_edge_arr;
  }

  // Algorithm maintains a map from edges in the current graph to original
  // edges (initially just identity).
  template <class Graph>
  inline pbbslib::dyn_arr<edge> SpanningForest(Graph& G, double beta = 0.2,
                                        bool pack = false, bool permute = false) {
    std::function<edge(edge)> identity_mapping = [&] (edge e) {
      return e;
    };
    return SpanningForest_Impl(G, beta, 0, identity_mapping, pack, permute);
  }

}  // namespace spanning_forest
