#pragma once

#include "pbbslib/random.h"
#include "ligra/bridge.h"
#include "ligra/macros.h"

namespace barabasi_albert {
auto generate_updates(
    size_t n,
    size_t edges_per_node=10) {
  size_t m = n*edges_per_node;
  auto edges = pbbs::sequence<std::pair<uintE, uintE>>(m);
  // Each edge picks a random id before it.
  for (size_t i=0; i<edges_per_node; i++) {
    edges[i] = std::make_pair(i, i+1);
  }
  pbbs::random rnd;
  parallel_for(1, n, [&] (size_t i) {
    auto i_rnd = rnd.fork(i);
    size_t mod_by = 2*i*edges_per_node;
    for (size_t j=0; j<edges_per_node; j++) {
      bool idx_ok = false;
      uintE idx;
      while (!idx_ok) {
        idx = i_rnd.rand() % mod_by;
        idx_ok = true;
        for (size_t k=0; k<j; k++) {
          if (edges[i*edges_per_node + k].second == idx) {
            idx_ok = false;
            i_rnd = i_rnd.next();
            break;
          }
        }
      }
      i_rnd = i_rnd.next();

      edges[i*edges_per_node + j] = std::make_pair(i, idx);
    }
  });
  // now force the m-edges_per_node unfilled slots
  // note that we can make this work-efficient and O(log n) depth (and it would
  // probably be pretty practical), but it is too much effort for now.
  for (size_t i=edges_per_node; i<m; i++) {
    auto& ref = edges[i];
    uintE offset = ref.second;
    uintE edge_idx = offset/2;
    if (offset % 2 == 0) { // 0 indexed
      ref.second = edges[edge_idx].first;
    } else {
      ref.second = edges[edge_idx].second;
    }
  }
  return edges;
}
} // namespace barabasi_albert

/* As in Kanat's paper, show a figure showing for each point in the stream the
 * number of connected components in the stream */
