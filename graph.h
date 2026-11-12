
#ifndef GRAPH_H
#define GRAPH_H

#include <unordered_map>
#include <unordered_set>
#include <random>
#include <vector>
#include <map>
#include "types.h"


class Graph{

public:

    std::vector<VertexID> vertices;
    
    ui vertices_count;
    ui edges_count;
    ui max_degree;

    ui* degrees;
    ui* main_degrees;


    ui* offsets;
    VertexID * neighbors;  
    

    Graph(){
        
        vertices_count = 0;
        edges_count = 0;
        max_degree = 0;

        offsets = NULL;
        neighbors = NULL;
    }

    ~Graph() {
        delete[] offsets;
        delete[] neighbors;
    }

public:

    void loadGraphFromFile(const std::string& file_path);

    void loadGraphFromFileForBothDirectionEdges(const std::string& file_path);

    bool is_smaller(VertexID u, VertexID v);
    
    void transformToAugmentedGraph(Graph* augmented_graph);
    void transformToAugmentedGraphForBothDirectionalEdges(Graph* augmented_graph);
    

    const ui* getOffsets() const {
        return offsets;
    }

    const ui* getNeighbors() const {
        return neighbors;
    }

    const ui getVerticesCount() const {
        return vertices_count;
    }

    const ui getEdgesCount() const {
        return edges_count;
    }

    const ui getVertexDegree(const VertexID id) const {
        return offsets[id + 1] - offsets[id];
    }


    const ui getGraphMaxDegree() const {
        return max_degree;
    }

    ui * getVertexNeighbors(const VertexID id, ui& count) const {
        count = offsets[id + 1] - offsets[id];
        return neighbors + offsets[id];
    }

    void fillSrcVertexArray(VertexID* src_vtx_array){

        ui edge_idx = 0;

        for(ui i = 0; i < vertices_count; i++){
            for(ui j = offsets[i]; j < offsets[i + 1]; j++){
                src_vtx_array[edge_idx] = i;
                edge_idx++;
            }
        }
    }

    void getNeighborCount(const VertexID id, ui& count) const {
        count = offsets[id + 1] - offsets[id];
    }

    void printNeighbors(const VertexID id) const {
        ui count = offsets[id + 1] - offsets[id];

        std::cout << "VID : " << id << " --- " ;

        for(ui i = 0; i < count; i++){
            std::cout << neighbors[offsets[id] + i] << "   "; 
        }

        std::cout << std::endl;
    }

    void print_graph_details(){

        double avg_degree, total_degree = 0;
        long long wedge_count = 0;
        ui max_degree = 0;

        for(ui i = 0; i < vertices_count; i++){
            wedge_count += (degrees[i] * (degrees[i] - 1))/ 2;
            total_degree += degrees[i];
            if(degrees[i] > max_degree){
                max_degree = degrees[i];
            }
        }

        avg_degree = (double) total_degree / vertices_count;

        std::cout << "Vertices Count : " << vertices_count << std::endl;
        std::cout << "Edges Count : " << edges_count << std::endl;
        std::cout << "Max Degree : " << max_degree << std::endl;
        std::cout << "Average Degree : " << avg_degree << std::endl;
        std::cout << "Wedge Count : " << wedge_count << std::endl;
    }

    bool checkDirectedEdgeExistence(VertexID u, VertexID v) const {
        
        ui count = 0;
        const VertexID* neighbors =  getVertexNeighbors(u, count);

        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == v) {
                return true;
            }
            else if (neighbors[mid] > v){
                end = mid - 1;
            }
            else{
                begin = mid + 1;
            }
        }

        return false;
    }


    void deleteAndClear(){

        delete[] offsets;
        delete[] neighbors;
    }

};


#endif