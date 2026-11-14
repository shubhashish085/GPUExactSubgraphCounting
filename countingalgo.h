#ifndef GSC_COUNTING_H
#define GSC_COUNTING_H

#include <map>
#include <algorithm>
#include "graph.h"
#include "types.h"

class CountingAlgorithm {

public:

    static long long db_count_triangle_parallel(Graph* augmented_graph);
<<<<<<< HEAD
    static long long db_count_square_parallel(Graph* augmented_graph);
=======
    static long long db_count_triangle_ei_parallel(Graph* augmented_graph, VertexID* edge_src_vertices);
    static long long db_count_triangle_warp_parallel(Graph* augmented_graph);

>>>>>>> main
};


#endif
