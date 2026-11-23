#ifndef GSC_COUNTING_H
#define GSC_COUNTING_H

#include <map>
#include <algorithm>
#include "graph.h"
#include "types.h"

extern int threads_per_block;

class CountingAlgorithm {

public:

    //Triangle
    static long long db_count_triangle_parallel(Graph* augmented_graph);
    static long long db_count_triangle_ei_parallel(Graph* augmented_graph, VertexID* edge_src_vertices);
    static long long db_count_triangle_warp_parallel(Graph* augmented_graph);
    static long long db_count_triangle_parallel_parameterized(Graph* augmented_graph, int thread_cnt);

    //Triangle - Atomic
    static unsigned long long db_count_triangle_parallel_atomic(Graph* augmented_graph);
    static unsigned long long db_count_triangle_ei_parallel_atomic(Graph* augmented_graph, VertexID* edge_src_vertices);
    static unsigned long long db_count_triangle_warp_parallel_atomic(Graph* augmented_graph);


    //Square
    static long long db_count_square_parallel(Graph* main_graph, Graph* augmented_graph);
    static long long db_count_square_ei_parallel (Graph* graph, VertexID* edge_src_vertices);

};


#endif
