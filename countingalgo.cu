#include "countingalgo.h"
#include <cuda_runtime.h>

#define WARP_SIZE 32

int threads_per_block = 128;


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

__device__ __forceinline__ int upper_bound_device_alt(const VertexID* arr, const ui start_idx, const ui end_idx, VertexID target) {
    ui idx = start_idx;
    
    while(idx < end_idx){
        if(arr[idx] > target){
            break;
        }
        idx++;
    }   

    return idx;
}


__device__ __forceinline__ int upper_bound_device(const VertexID* arr, const ui start_idx, const ui end_idx, VertexID target) {
    ui left = start_idx, mid;
    ui right = end_idx - 1;

    while (left < right) {

        mid = left + (right - left) / 2;

        if (arr[mid] <= target) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    return left;
}


__device__ __forceinline__ ui array_intersection_gpu(const VertexID* d_neighbors, VertexID start_u, VertexID end_u, VertexID start_v, VertexID end_v){

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


__device__ ui array_intersection_with_ref_vtx_gpu(const ui* d_main_degrees, const VertexID* d_main_neighbors, VertexID start_u, VertexID end_u, VertexID start_v, VertexID end_v, VertexID v1, VertexID v3, VertexID v4){

    ui count = 0;
    VertexID i = start_u, j = start_v;

    VertexID a, b;

    while (i < end_u && j < end_v) {

        a = d_main_neighbors[i];
        b = d_main_neighbors[j];

        if (a == b) {

            if(is_smaller_gpu(d_main_degrees, v1, a) && is_smaller_gpu(d_main_degrees, a, v3) && is_smaller_gpu(d_main_degrees, a, v4)){
                count++;
            }
            
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


/*__device__ __forceinline__ ui array_intersection_gpu (const VertexID* d_neighbors, const ui array_1_start, const ui array_1_end, const ui array_2_start, const ui array_2_end) {
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
}*/


__device__ bool check_directed_edge_existence_gpu(const ui* d_offsets, const VertexID* d_neighbors, VertexID u, VertexID v){

    ui count = 0;

    int begin = d_offsets[u];
    int end = d_offsets[u + 1] - 1;
    
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



__global__ void square_counting_kernel (const ui* d_vertices_count, const ui* d_edges_count, const ui* d_main_degrees, const ui* d_main_offsets, const ui* d_main_neighbors, const ui* d_offsets, const VertexID* d_neighbors, ui* d_square_counts)
{
    // GPU triangle counting implementation
    long long exact_count = 0;
    ui count_1 = 0, count_2 = 0, count_3 = 0;
    
    VertexID v1 = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(v1 >= *(d_vertices_count)){
        return;
    }

    VertexID v2, v3, v4;
    ui v1_start = d_offsets[v1], v2_start, v3_start, v4_start;
    ui v1_end = d_offsets[v1 + 1], v2_end, v3_end, v4_end;



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

    ui v1_degree = d_main_degrees[v1];

    for (ui j = v1_start; j < v1_end; j++){
        
        v3 = d_neighbors[j];
        v3_start = d_main_offsets[v3];
        v3_end = d_main_offsets[v3 + 1];

        for (ui k = j+1; k < v1_end; k++){
            v4 = d_neighbors[k];
            v4_start = d_main_offsets[v4];
            v4_end = d_main_offsets[v4 + 1];

            count_3 = array_intersection_with_ref_vtx_gpu(d_main_degrees, d_main_neighbors, v3_start, v3_end, v4_start, v4_end, v1, v3, v4);
            exact_count += count_3;
        }
    }

    d_square_counts[v1] = exact_count;     
}


__global__ void square_counting_ei_kernel (const ui* d_vertices_count, const ui* d_edges_count, const ui* d_offsets, const VertexID* d_src_vertices, const VertexID* d_neighbors, ui* d_square_counts)
{

    // GPU triangle counting implementation (EdgeIterator)
    ui count = 0;    
    VertexID e_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(e_idx >= ((*(d_edges_count)) * 2) || (d_src_vertices[e_idx] >= d_neighbors[e_idx])){
        return;
    }

    VertexID u, v, start_u, start_v, end_u, end_v;
    ui upper_index_u, upper_index_v;

    u = d_src_vertices[e_idx];
    v = d_neighbors[e_idx];

    start_u = d_offsets[u];
    end_u   = d_offsets[u + 1];

    start_v = d_offsets[v];
    end_v   = d_offsets[v + 1];

    upper_index_u = upper_bound_device_alt(d_neighbors, start_u, end_u, v);
    upper_index_v = upper_bound_device_alt(d_neighbors, start_v, end_v, u);

    //printf("Upper Index u : %d, Upper Index v : %d , end u : %d , end v : %d \n", upper_index_u, upper_index_v, end_u, end_v);

    VertexID v3, v4;

    for (ui i = upper_index_u; i < end_u; i++){
        v3 = d_neighbors[i];
        for (ui j = upper_index_v; j < end_v; j++){
            v4 = d_neighbors[j];
            if((v3 != v4) && check_directed_edge_existence_gpu(d_offsets, d_neighbors, v3, v4)){
                count += 1;
            }
        }
    }

    //printf("Count : %d\n", count);

    d_square_counts[e_idx] = count;     
}

// Thread Per Vertex
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

    //printf("Vertex ID : %d - Triangle Count : %d\n", v, exact_count);

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


// Thread Per Vertex - Triangle
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


    //int threads_per_block = 256;
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


    //int threads_per_block = 256;
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


// Thread Per Vertex - Square
long long CountingAlgorithm::db_count_square_parallel(Graph* main_graph, Graph* augmented_graph){


    ui* d_main_degrees = nullptr, *d_main_offsets = nullptr, *d_offsets = nullptr, *d_main_neighbors = nullptr, *d_neighbors = nullptr, *d_square_counts = nullptr;
    ui* d_vertices_count = nullptr, *d_edges_count = nullptr;
    ui* square_counts = new ui[augmented_graph->vertices_count];

    long long exact_count = 0;
    

    cudaMalloc((void**)&d_vertices_count, sizeof(ui));
    cudaMalloc((void**)&d_edges_count, sizeof(ui));    
    cudaMalloc((void**)&d_main_degrees, (augmented_graph->vertices_count * sizeof(ui)));
    cudaMalloc((void**)&d_main_offsets, ((augmented_graph->vertices_count + 1)  * sizeof(ui)));
    cudaMalloc((void**)&d_offsets, ((augmented_graph->vertices_count + 1) * sizeof(ui)));
    cudaMalloc((void**)&d_main_neighbors, ((augmented_graph->edges_count) * 2 * sizeof(VertexID)));
    cudaMalloc((void**)&d_neighbors, ((augmented_graph->edges_count) * sizeof(VertexID)));
    cudaMalloc((void**)&d_square_counts, ((augmented_graph->vertices_count) * sizeof(ui)));


    cudaMemcpy(d_vertices_count, &(augmented_graph->vertices_count), sizeof(ui), cudaMemcpyHostToDevice);
    cudaMemcpy(d_edges_count, &(augmented_graph->edges_count), sizeof(ui), cudaMemcpyHostToDevice);
    cudaMemcpy(d_main_degrees, augmented_graph->main_degrees, ((augmented_graph->vertices_count) * sizeof(ui)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_main_offsets, main_graph->offsets, ((augmented_graph->vertices_count + 1) * sizeof(ui)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_offsets, augmented_graph->offsets, ((augmented_graph->vertices_count + 1) * sizeof(ui)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_main_neighbors, main_graph->neighbors, ((augmented_graph->edges_count) * 2 * sizeof(VertexID)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighbors, augmented_graph->neighbors, ((augmented_graph->edges_count) * sizeof(VertexID)), cudaMemcpyHostToDevice);


    int num_blocks = (augmented_graph->vertices_count + threads_per_block - 1) / threads_per_block;

    square_counting_kernel<<<num_blocks, threads_per_block>>>(
        d_vertices_count,
        d_edges_count,
        d_main_degrees,
        d_main_offsets,
        d_main_neighbors,
        d_offsets,
        d_neighbors,
        d_square_counts
    );
    
    cudaDeviceSynchronize();
    cudaMemcpy(square_counts, d_square_counts, (augmented_graph->vertices_count * sizeof(ui)), cudaMemcpyDeviceToHost);

    for(ui i = 0; i < augmented_graph->vertices_count; i++){
        exact_count += square_counts[i];
    }

    cudaFree(d_main_degrees);
    cudaFree(d_main_offsets);
    cudaFree(d_offsets);
    cudaFree(d_main_neighbors);
    cudaFree(d_neighbors);
    cudaFree(d_square_counts); 

    return exact_count;
}


// Thread Per Edge - Square
long long CountingAlgorithm::db_count_square_ei_parallel(Graph* graph, VertexID* edge_src_vertices){

    ui* d_offsets = nullptr, *d_src_vertices = nullptr, *d_neighbors = nullptr, *d_square_counts = nullptr;
    ui* d_vertices_count = nullptr, *d_edges_count = nullptr;
    ui* square_counts = new ui[(graph->edges_count) * 2];

    long long exact_count = 0;
    

    cudaMalloc((void**)&d_vertices_count, sizeof(ui));
    cudaMalloc((void**)&d_edges_count, sizeof(ui));    
    cudaMalloc((void**)&d_offsets, ((graph->vertices_count + 1) * sizeof(ui)));
    cudaMalloc((void**)&d_src_vertices, ((graph->edges_count) * 2 * sizeof(VertexID)));
    cudaMalloc((void**)&d_neighbors, ((graph->edges_count) * 2 * sizeof(VertexID)));
    cudaMalloc((void**)&d_square_counts, ((graph->edges_count) * 2 * sizeof(ui)));


    cudaMemcpy(d_vertices_count, &(graph->vertices_count), sizeof(ui), cudaMemcpyHostToDevice);
    cudaMemcpy(d_edges_count, &(graph->edges_count), sizeof(ui), cudaMemcpyHostToDevice);
    cudaMemcpy(d_offsets, graph->offsets, ((graph->vertices_count + 1) * sizeof(ui)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_src_vertices, edge_src_vertices, ((graph->edges_count) * 2 * sizeof(VertexID)), cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighbors, graph->neighbors, ((graph->edges_count) * 2 * sizeof(VertexID)), cudaMemcpyHostToDevice);

    cudaMemset(d_square_counts, 0, ((graph->edges_count) * 2 * sizeof(ui)));


    //int threads_per_block = 256;
    int num_blocks = (((graph->edges_count) * 2) + threads_per_block - 1) / threads_per_block;
    
    square_counting_ei_kernel<<<num_blocks, threads_per_block>>>(
        d_vertices_count,
        d_edges_count,
        d_offsets,
        d_src_vertices,
        d_neighbors,
        d_square_counts
    );
    
    cudaDeviceSynchronize();
    cudaMemcpy(square_counts, d_square_counts, (graph->edges_count * sizeof(ui)), cudaMemcpyDeviceToHost);

    ui end_loop_idx = (graph->edges_count) * 2;

    for(ui i = 0; i < end_loop_idx; i++){
        exact_count += square_counts[i];
    }

    cudaFree(d_offsets);
    cudaFree(d_src_vertices);
    cudaFree(d_neighbors);
    cudaFree(d_square_counts); 

    return exact_count;
}