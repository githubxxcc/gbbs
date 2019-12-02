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

#include "ligra/ligra.h"

#include <utility>
#include <vector>
#include "float.h"
#include <math.h>
#include <unordered_set>
#include <set>
#include <utility>


template<typename T>
class Vecs {
    public:
        T *data;
        size_t n_dim;
        size_t n_vecs;

        Vecs(T* _data, size_t _n_vecs, size_t _n_dim) :data(_data), n_dim(_n_dim), n_vecs(_n_vecs) {}
};



static float l2dist(uintE s, const std::vector<float> &query, const Vecs<float> &fvecs);

template<typename T>
static std::vector<T> vertex_vecs(uintE v, const Vecs<T> &fvecs);

static float l2dist(uintE s, const std::vector<float> &query, const Vecs<float> &fvecs);

template<typename T>
static void load_Tvecs(const char *filename, T *&data, size_t &num, size_t &dim);

template <class Graph>
inline std::set<std::pair<float, uintE>> GreedySearch(Graph& G, uintE src, std::vector<float> query, const Vecs<float> &fvecs, size_t L) {
    using W = typename Graph::weight_type;

    // Visited set
    std::unordered_set<uintE> Visited;

    // Candidate set
    std::set<std::pair<float, uintE>> Candidates;
    Candidates.insert(std::make_pair(l2dist(src, query, fvecs), src));
    std::unordered_set<uintE> CandidateSet;
    CandidateSet.insert(src);

    // Frontier: starting from the source
    vertexSubset Frontier(G.n, src);
    uintE next[L] = {0};

    while (!Frontier.isEmpty()) {

        uintE p_star = Frontier.vtx(0);
        float min_dist = FLT_MAX;
        for(auto i = 0; i < Frontier.size(); i++) {
            uintE v = Frontier.vtx(i);
            float dist = l2dist(v, query, fvecs);
            if(dist < min_dist) {
                min_dist = dist;
                p_star = v;
            }
        }

        //std::cout << "dis: (" << p_star << " , " << min_dist << " )" << std::endl;

        //parallel_for(0, Frontier.size(), [&](size_t i) {
        //  uintE v = Frontier.vtx(i);
        //  float dist = l2dist(v, query);
        //  auto dist_id = dist_and_id(dist, v);
        //  pbbs::write_min(&min_dist, dist_id, less_dist_id);
        //});
        // read value of min_dist to get id.
        // next step maps over all outgoing neighbors of min dist guy.

        // Updated visited
        Visited.insert(p_star);

        // Updated candidates
        //
        auto map_f = [&](const uintE &vtx_id, const uintE &ngh, const W &w) {
            if(CandidateSet.count(ngh)) {
                return;
            }
            // Calculate distance
            auto dist = l2dist(ngh, query, fvecs);
            if(Candidates.size() == L) {
                auto last = std::prev(Candidates.end());
                CandidateSet.erase(last->second);
                Candidates.erase(last);
            }
            // Update Candidates
            Candidates.insert(std::make_pair(dist, ngh));
            CandidateSet.insert(ngh);
        };
        G.get_vertex(p_star).mapOutNgh(p_star, map_f);

        // Find out next frontier 
        size_t idx = 0;
        for(auto itr = CandidateSet.begin(); itr != CandidateSet.end(); itr++) {
            //assert(itr->second < fvecs.n_vecs);
            if(!Visited.count(*itr)) {
                next[idx++] = *itr;
            }
        }

        vertexSubset nextFrontier(G.n, idx, next);

        Frontier = nextFrontier;
    }

    return Candidates;
}

// calculate l2 distance from a vertex s to the query.
static float l2dist(uintE s, const std::vector<float> &query, const Vecs<float> &fvecs) {
    auto sfvec = vertex_vecs<float>(s, fvecs);
    float dis2 = 0;
    for(auto i = 0U; i < fvecs.n_dim; i++) {
        dis2 += std::pow(sfvec[i] - query[i], 2.0);
    }

    return std::sqrt(dis2);
}

template<typename T>
void load_Tvecs(const char *filename, T *&data, size_t &num, size_t &dim) {
    std::ifstream in(filename, std::ios::binary | std::ios::ate);

    if (!in.is_open()) {

        std::cout << "Error opening file: " << filename << std::endl;

        exit(-1);

    }

    std::size_t fsize = in.tellg();

    char* temp = new char[fsize];

    in.close();

    in.open(filename, std::ios::binary);

    in.read(temp, fsize);

    dim = *(unsigned*)(temp);

    size_t vec_size = dim*sizeof(T);

    size_t dvec_size = vec_size + sizeof(unsigned);
    size_t nvecs = fsize / dvec_size;
    num = nvecs;
    //std::cout << "nvecs: " << nvecs << ", ndims: " << dim << "\n";
    data = new T[nvecs * dim];
    for(size_t i=0;i<nvecs;i++){
        memcpy(data + (i * dim), temp + i*dvec_size + sizeof(unsigned), vec_size);
    }
    delete[] temp;
    in.close();

}


////Credit: Suhas Jayaram Subramanya
//template<typename T>
//static void load_Tvecs(const char *filename, T *&data, size_t &num, size_t &dim) {
//    //check validity of file
//    std::ifstream in(filename, std::ios::binary | std::ios::ate);
//    if (!in.is_open()) {
//        std::cout << "Error opening file: " << filename << std::endl;
//        exit(-1);
//    }
//    std::size_t fsize = in.tellg();
//    in.seekg(0, std::ios::beg);
//    unsigned dim_u32;
//    in.read((char *) &dim_u32, sizeof(unsigned));
//    in.close();
//    dim = dim_u32;
//
//    std::size_t ndims = (std::size_t) dim;
//    std::size_t disk_vec_size = ndims * sizeof(T) + sizeof(unsigned);
//    std::size_t mem_vec_size = ndims * sizeof(T);
//    std::size_t npts = fsize / disk_vec_size;
//    num = npts;
//    std::cout << "Tvecs: " << filename << ", npts: " << npts
//        << ", ndims: " << ndims << "\n";
//    // allocate memory
//    data = new T[npts * ndims];
//
//    in.seekg(0, std::ios::beg);
//    unsigned dummy_ndims;
//    for (std::size_t i = 0; i < npts; i++) {
//        T *cur_vec = data + (i * ndims);
//        // read and ignore dummy ndims
//        in.read((char *) &dummy_ndims, sizeof(unsigned));
//
//        // read vec
//        in.read((char *) cur_vec, mem_vec_size);
//    }
//    return;
//}

template<typename T>
static std::vector<T> vertex_vecs(uintE v, const Vecs<T> &fvecs) {
    std::vector<T> tmp;
    tmp.reserve(fvecs.n_dim);

    std::copy(fvecs.data + v * fvecs.n_dim, fvecs.data + (v+1) * fvecs.n_dim, std::back_inserter(tmp));
    return tmp;
}

