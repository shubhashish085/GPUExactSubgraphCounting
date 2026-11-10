#include "countingalgo.h"
#include <cuda_runtime.h>


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