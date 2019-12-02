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

// Usage:
// numactl -i all ./GreedySearch -src 10012 -s -m -rounds 3 twitter_SJ
// flags:
//   required:
//     -src: the source to compute the GreedySearch from
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#include "GreedySearch.h"

template <class Graph>
double GreedySearch_runner(Graph& G, commandLine P) {
  uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
  std::cout << "### Application: GreedySearch" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -src = " << src << std::endl;
  std::cout << "### ------------------------------------" << endl;

  timer t; t.start();
  //FIXME;
  //TODO: parse the query Q
  
  std::string fvec_file  = P.getOptionValue("-fvec", "nope");
  std::string gt_file = P.getOptionValue("-gt", "gt");
  std::string query_file = P.getOptionValue("-query", "guery");
  size_t L = 101;



  // Parse coordinates
  float * data = nullptr;
  size_t n_vecs, n_dims;

  load_Tvecs<float>(fvec_file.c_str(), data, n_vecs, n_dims);
  Vecs<float> fvecs(data, n_vecs, n_dims);

  // Parse the groundtruth
  int *gdata = nullptr;
  size_t g_n_vecs, g_n_dims;
  load_Tvecs<int>(gt_file.c_str(), gdata, g_n_vecs, g_n_dims);
 

  Vecs<int> gts(gdata, g_n_vecs, g_n_dims);

  // Parse the query
  float * qdata = nullptr;
  size_t q_n_vecs, q_n_dims;
  load_Tvecs<float>(query_file.c_str(), qdata, q_n_vecs, q_n_dims);
  Vecs<float> queries(qdata, q_n_vecs, q_n_dims);
 

  auto query = vertex_vecs(0, queries);
  auto results = GreedySearch(G, src, query, fvecs, L);
  auto gt = vertex_vecs(0, gts);
  std::unordered_set<int> gt_set(gt.begin(), gt.end());  
    
  size_t Union = gt_set.size(),Inter = 0;
  for(auto itr = results.begin(); itr != results.end(); itr++) {
      //std::cout << " [ " << itr->first << " , " << itr->second << std::endl;
        if(gt_set.count(static_cast<int>(itr->second))) {
            Inter++;
        } else {
            Union++;
        }
  }

  std::cout << "Recall: " << (1.0 * Inter) / Union << std::endl;

  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

generate_main(GreedySearch_runner, false);
