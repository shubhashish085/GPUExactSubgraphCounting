#include "countingalgo.h"
#include <cuda_runtime.h>

#define WARP_SIZE 32


__device__ ui array_intersection_gpu(const VertexID* d_neighbors, VertexID start_u, VertexID end_u, VertexID start_v, VertexID end_v){

    ui count = 0;
    VertexID i = start_u, j = start_v;

    VertexID a, b;

    while (i < end_u && j < end_v) {

        a = d_neighbors[i];
        b = d_neighbors[j];

        if (a == b) {
            count++;
            i++;
            j++;
        }
        else if (a < b) {
            i++;
        }
        else {
            j++;
        }
    }

    return count;
}


__device__ __forceinline__ ui array_intersection_gpu (const VertexID* d_neighbors, const ui array_1_start, const ui array_1_end, const ui array_2_start, const ui array_2_end) {
    ui i = array_1_start, j = array_2_start;
    ui count = 0;

    while (i < array_1_end && j < array_2_end) {
        if (d_neighbors[i] == d_neighbors[j]) {
            count++;
            i++;
            j++;
        } 
        else if (d_neighbors[i] < d_neighbors[j]) {
            i++;
        } 
        else {
            j++;
        }
    }

    return count;
}


__device__ bool check_directed_edge_existence_gpu(const ui* d_offsets, const VertexID* d_neighbors, VertexID u, VertexID v){

    ui count = 0;

    int begin = d_offsets[u];
    int end = d_offsets[u + 1];
    
    while (begin <= end) {
    
        int mid = begin + ((end - begin) >> 1);

        if (d_neighbors[mid] == v) {
            return true;
        }
        else if (d_neighbors[mid] > v){
            end = mid - 1;
        }
        else{
            begin = mid + 1;
        }
    }

    return false;
}


__device__ bool is_smaller_gpu(const ui* d_main_degrees, VertexID u, VertexID v){

    if(u == v){
        return false;
    }

    if(d_main_degrees[u] < d_main_degrees[v]){
        return true;
    }

    if((d_main_degrees[u] == d_main_degrees[v]) && (u < v)){
        return true;
    }

    return false;
}


__global__ void square_counting_kernel (const ui* d_vertices_count, const ui* d_edges_count, const ui* d_main_degrees, const ui* d_degrees, const ui* d_offsets, const VertexID* d_neighbors, ui* d_square_counts)
{
    // GPU triangle counting implementation
    long long exact_count = 0;
    ui count_1 = 0, count_2 = 0, count_3 = 0;
    
    VertexID v1 = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(v1 >= *(d_vertices_count)){
        return;
    }

    VertexID v2, v3, v4;
    ui v1_start = d_offsets[v1], v2_start, v3_start;
    ui v1_end = d_offsets[v1 + 1], v2_end, v3_end;

    // 1243 - Square
    

    for (VertexID j = v1_start; j < v1_end; j++){

        v2 = d_neighbors[j];
        v2_start = d_offsets[v2];
        v2_end = d_offsets[v2 + 1];

        for (VertexID k = v1_start; k < v1_end; k++){
            
            count_1 = 0;
            v3 = d_neighbors[k];

            if(!is_smaller_gpu(d_main_degrees, v2, v3)){
                continue;
            }

            v3_start = d_offsets[v3];
            v3_end = d_offsets[v3 + 1];

            count_1 = array_intersection_gpu (d_neighbors, v2_start, v2_end, v3_start, v3_end);
            exact_count += count_1;
        }
    }

    //std::cout << "Count 1 : " << exact_count << std::endl;
    // 1234 - Square

    for (VertexID j = v1_start; j < v1_end; j++)
    {
        v2 = d_neighbors[j];
        v2_start = d_offsets[v2];
        v2_end = d_offsets[v2 + 1];

        for (VertexID k = v2_start; k < v2_end; k++)
        {
            count_2 = 0;
            v3 = d_neighbors[k];
            v3_start = d_offsets[v3];
            v3_end = d_offsets[v3 + 1];
            
            count_2 = array_intersection_gpu (d_neighbors, v1_start, v1_end, v3_start, v3_end);
            exact_count += count_2;
        }
    }

    //std::cout << "Count 2 : " << exact_count << std::endl;
    // 1324 - Square

    std::map<std::pair<VertexID, VertexID>, ui> wedge_map;
    std::pair<VertexID, VertexID> search_pair;

    for (VertexID v1 = 0; v1 < graph->vertices_count; v1++)
    {
        nbrs_1 = graph->getVertexNeighbors(v1, nbrs_1_cnt);

        for (ui j = 0; j < nbrs_1_cnt; j++)
        {
            for (ui k = j+1; k < nbrs_1_cnt; k++)
            {
                search_pair = std::make_pair(std::min(nbrs_1[j], nbrs_1[k]), std::max(nbrs_1[j], nbrs_1[k]));
                auto search = wedge_map.find(search_pair);
                if (search == wedge_map.end())
                {
                    wedge_map[search_pair] = 1;
                }
                else
                {
                    wedge_map[search_pair] = wedge_map[search_pair] + 1;
                }
            }
        }
    }

    for (auto const &[key, value] : wedge_map){
        count_3 += (value * (value - 1)) / 2;
    }

    //std::cout << "Count 3 : " << count_3 << std::endl;

    exact_count += count_3;

    d_square_counts[v] = exact_count;     
}


__global__ void triangle_counting_kernel (const ui* d_vertices_count, const ui* d_edges_count, const ui* d_main_degrees, const ui* d_degrees, const ui* d_offsets, const VertexID* d_neighbors, ui* d_triangle_counts)
{
    // GPU triangle counting implementation
    long long exact_count = 0;
    
    VertexID v = (blockIdx.x * blockDim.x) + threadIdx.x;


    if(v >= *(d_vertices_count)){
        return;
    }

    VertexID u, w;

    ui start = d_offsets[v];
    ui end = d_offsets[v + 1];
    
    for (ui i = start; i < end; i++){
        u = d_neighbors[i];
        for(ui j = start; j < end; j++){
            w = d_neighbors[j];
            if(is_smaller_gpu(d_main_degrees, u, w) && check_directed_edge_existence_gpu(d_offsets, d_neighbors, u, w)){
                exact_count++;
            }
        }
    }

    printf("Vertex ID : %d - Triangle Count : %d\n", v, exact_count);

    d_triangle_counts[v] = exact_count;     
}


__global__ void triangle_counting_ei_kernel (const ui* d_vertices_count, const ui* d_edges_count, const ui* d_offsets, const VertexID* d_src_vertices, const VertexID* d_neighbors, ui* d_triangle_counts)
{
    // GPU triangle counting implementation (EdgeIterator)    
    VertexID e_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(e_idx >= *(d_edges_count)){
        return;
    }

    VertexID u, v, start_u, start_v, end_u, end_v;

    u = d_src_vertices[e_idx];
    v = d_neighbors[e_idx];

    start_u = d_offsets[u];
    end_u   = d_offsets[u + 1];

    start_v = d_offsets[v];
    end_v = d_offsets[v + 1];

    ui intersection_length = array_intersection_gpu(d_neighbors, start_u, end_u, start_v, end_v);

    d_triangle_counts[e_idx] = intersection_length;     
}


__global__ void triangle_counting_warp_kernel(const ui* d_vertices_count, const ui* d_edges_count, const ui* d_main_degrees, const ui* d_degrees, const ui* d_offsets, const VertexID* d_neighbors, ui* d_triangle_counts)
{
    const ui lane = threadIdx.x & (WARP_SIZE - 1);

    const ui warp_id = threadIdx.x / WARP_SIZE;

    const VertexID u = blockIdx.x * (blockDim.x / WARP_SIZE) + warp_id;

    if (u >= *(d_vertices_count)){
        return;
    }

    VertexID start_u = d_offsets[u];
    VertexID end_u   = d_offsets[u + 1];

    unsigned long long local_count = 0;
    VertexID start_v, end_v;
    ui intersection_length;

    for (VertexID i = start_u + lane; i < end_u; i += WARP_SIZE)
    {
        VertexID v = d_neighbors[i];    
        start_v = d_offsets[v];
        end_v   = d_offsets[v + 1];

        intersection_length = array_intersection_gpu(d_neighbors, start_u, end_u, start_v, end_v);
        local_count += intersection_length;
    }

    // Warp-level reduction.
    for (int offset = WARP_SIZE / 2; offset > 0; offset >>= 1){
        local_count += __shfl_down_sync( 0xffffffff, local_count, offset);
    }

    if (lane == 0) {
        d_triangle_counts[u] = local_count;
    }
}


// Thread Per Vertex
long long CountingAlgorithm::db_count_triangle_parallel(Graph* augmented_graph){


    ui* d_degrees = nullptr, *d_main_degrees = nullptr, *d_offsets = nullptr, *d_neighbors = nullptr, *d_triangle_counts = nullptr;
    ui* d_vertices_count = nullptr, *d_edges_count = nullptr;
    ui* triangle_counts = new ui[augmented_graph->vertices_count];

    long long exact_count = 0;
    

    cudaMalloc((void**)&d_vertices_count, sizeof(ui));
    cudaMalloc((void**)&d_edges_count, sizeof(ui));    
    cudaMalloc((void**)&d_degrees, (augmented_graph->vertices_count * sizeof(ui)));
    cudaMalloc((void**)&d_main_degrees, (augmented_graph->vertices_count * sizeof(ui)));
    cudaMalloc((void**)&d_offsets, ((augmented_graph->vertices_count + 1) * sizeof(ui)));
    cudaMalloc((void**)&d_neighbors, (augmented_graph->edges_count * sizeof(VertexID)));
    cudaMalloc((void**)&d_triangle_counts, (augmented_graph->vertices_count * sizeof(ui)));


    cudaMemcpy(d_vertices_count, &(augmented_graph->vertices_count), sizeof(ui), cudaMemcpyHostToDevice);
    cudaMemcpy(d_edges_count, &(augmented_graph->edges_count), sizeof(ui), cudaMemcpyHostToDevice);
    cudaMemcpy(d_degrees, augmented_graph->degrees, (augmented_graph->vertices_count * sizeof(ui)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_main_degrees, augmented_graph->main_degrees, (augmented_graph->vertices_count * sizeof(ui)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_offsets, augmented_graph->offsets, ((augmented_graph->vertices_count + 1) * sizeof(ui)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighbors, augmented_graph->neighbors, ((augmented_graph->edges_count) * sizeof(VertexID)), cudaMemcpyHostToDevice);


    int threads_per_block = 256;
    int num_blocks = (augmented_graph->vertices_count + threads_per_block - 1) / threads_per_block;
    
    triangle_counting_kernel<<<num_blocks, threads_per_block>>>(
        d_vertices_count,
        d_edges_count,
        d_main_degrees,
        d_degrees,
        d_offsets,
        d_neighbors,
        d_triangle_counts
    );
    
    cudaDeviceSynchronize();
    cudaMemcpy(triangle_counts, d_triangle_counts, (augmented_graph->vertices_count * sizeof(ui)), cudaMemcpyDeviceToHost);

    for(ui i = 0; i < augmented_graph->vertices_count; i++){
        exact_count += triangle_counts[i];
    }

    cudaFree(d_degrees);
    cudaFree(d_main_degrees);
    cudaFree(d_offsets);
    cudaFree(d_neighbors);
    cudaFree(d_triangle_counts); 

    return exact_count;
}

//Warp Per Vertex
long long CountingAlgorithm::db_count_triangle_warp_parallel(Graph* augmented_graph){

    ui* d_degrees = nullptr, *d_main_degrees = nullptr, *d_offsets = nullptr, *d_neighbors = nullptr, *d_triangle_counts = nullptr;
    ui* d_vertices_count = nullptr, *d_edges_count = nullptr;
    ui* triangle_counts = new ui[augmented_graph->vertices_count];

    long long exact_count = 0;
    

    cudaMalloc((void**)&d_vertices_count, sizeof(ui));
    cudaMalloc((void**)&d_edges_count, sizeof(ui));    
    cudaMalloc((void**)&d_degrees, (augmented_graph->vertices_count * sizeof(ui)));
    cudaMalloc((void**)&d_main_degrees, (augmented_graph->vertices_count * sizeof(ui)));
    cudaMalloc((void**)&d_offsets, ((augmented_graph->vertices_count + 1) * sizeof(ui)));
    cudaMalloc((void**)&d_neighbors, (augmented_graph->edges_count * sizeof(VertexID)));
    cudaMalloc((void**)&d_triangle_counts, (augmented_graph->vertices_count * sizeof(ui)));


    cudaMemcpy(d_vertices_count, &(augmented_graph->vertices_count), sizeof(ui), cudaMemcpyHostToDevice);
    cudaMemcpy(d_edges_count, &(augmented_graph->edges_count), sizeof(ui), cudaMemcpyHostToDevice);
    cudaMemcpy(d_degrees, augmented_graph->degrees, (augmented_graph->vertices_count * sizeof(ui)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_main_degrees, augmented_graph->main_degrees, (augmented_graph->vertices_count * sizeof(ui)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_offsets, augmented_graph->offsets, ((augmented_graph->vertices_count + 1) * sizeof(ui)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighbors, augmented_graph->neighbors, ((augmented_graph->edges_count) * sizeof(VertexID)), cudaMemcpyHostToDevice);


    int threads_per_block = 256;
    int warps_per_block = threads_per_block / WARP_SIZE;
    int num_blocks = (augmented_graph->vertices_count + warps_per_block - 1) / warps_per_block;
    
    triangle_counting_warp_kernel<<<num_blocks, threads_per_block>>>(
        d_vertices_count,
        d_edges_count,
        d_main_degrees,
        d_degrees,
        d_offsets,
        d_neighbors,
        d_triangle_counts
    );
    
    cudaDeviceSynchronize();
    cudaMemcpy(triangle_counts, d_triangle_counts, (augmented_graph->vertices_count * sizeof(ui)), cudaMemcpyDeviceToHost);

    for(ui i = 0; i < augmented_graph->vertices_count; i++){
        exact_count += triangle_counts[i];
    }

    cudaFree(d_degrees);
    cudaFree(d_main_degrees);
    cudaFree(d_offsets);
    cudaFree(d_neighbors);
    cudaFree(d_triangle_counts); 

    return exact_count;
}


// Thread Per Edge
long long CountingAlgorithm::db_count_triangle_ei_parallel(Graph* augmented_graph, VertexID* edge_src_vertices){

    ui* d_offsets = nullptr, *d_src_vertices = nullptr, *d_neighbors = nullptr, *d_triangle_counts = nullptr;
    ui* d_vertices_count = nullptr, *d_edges_count = nullptr;
    ui* triangle_counts = new ui[augmented_graph->edges_count];

    long long exact_count = 0;
    

    cudaMalloc((void**)&d_vertices_count, sizeof(ui));
    cudaMalloc((void**)&d_edges_count, sizeof(ui));    
    cudaMalloc((void**)&d_offsets, ((augmented_graph->vertices_count + 1) * sizeof(ui)));
    cudaMalloc((void**)&d_src_vertices, (augmented_graph->edges_count * sizeof(VertexID)));
    cudaMalloc((void**)&d_neighbors, (augmented_graph->edges_count * sizeof(VertexID)));
    cudaMalloc((void**)&d_triangle_counts, (augmented_graph->edges_count * sizeof(ui)));


    cudaMemcpy(d_vertices_count, &(augmented_graph->vertices_count), sizeof(ui), cudaMemcpyHostToDevice);
    cudaMemcpy(d_edges_count, &(augmented_graph->edges_count), sizeof(ui), cudaMemcpyHostToDevice);
    cudaMemcpy(d_offsets, augmented_graph->offsets, ((augmented_graph->vertices_count + 1) * sizeof(ui)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_src_vertices, edge_src_vertices, ((augmented_graph->edges_count) * sizeof(VertexID)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighbors, augmented_graph->neighbors, ((augmented_graph->edges_count) * sizeof(VertexID)), cudaMemcpyHostToDevice);


    int threads_per_block = 256;
    int num_blocks = (augmented_graph->edges_count + threads_per_block - 1) / threads_per_block;
    
    triangle_counting_ei_kernel<<<num_blocks, threads_per_block>>>(
        d_vertices_count,
        d_edges_count,
        d_offsets,
        d_src_vertices,
        d_neighbors,
        d_triangle_counts
    );
    
    cudaDeviceSynchronize();
    cudaMemcpy(triangle_counts, d_triangle_counts, (augmented_graph->edges_count * sizeof(ui)), cudaMemcpyDeviceToHost);

    for(ui i = 0; i < augmented_graph->edges_count; i++){
        exact_count += triangle_counts[i];
    }

    cudaFree(d_offsets);
    cudaFree(d_src_vertices);
    cudaFree(d_neighbors);
    cudaFree(d_triangle_counts); 

    return exact_count;
}
