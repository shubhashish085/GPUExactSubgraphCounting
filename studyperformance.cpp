#include <chrono>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>


#include "matchingcommand.h"
#include "types.h"
#include "countingalgo.h"


// Triangle Count - Degree Based ID - For GPU
/*int main(int argc, char** argv){

    std::string input_data_graph_file = argv[1];

    clock_t read_begin_clock = clock();
    Graph* graph = new Graph();
    graph->loadGraphFromFile(input_data_graph_file);
    //graph->loadGraphFromFileForBothDirectionEdges(input_data_graph_file);
    double input_read_time = (double(clock() - read_begin_clock)) / CLOCKS_PER_SEC;    


    clock_t transformation_begin_clock = clock();
    Graph* augmented_graph = new Graph();
    graph->transformToAugmentedGraph(augmented_graph);
    //graph->transformToAugmentedGraphForBothDirectionalEdges(augmented_graph);
    double transformation_time = (double(clock() - transformation_begin_clock)) / CLOCKS_PER_SEC;

    

    clock_t counting_begin_clock = clock();
    long long exact_count = CountingAlgorithm::db_count_triangle_parallel(augmented_graph);  
    double counting_time = (double(clock() - counting_begin_clock)) / CLOCKS_PER_SEC;
    double total_time = input_read_time + transformation_time + counting_time;
   

    std::cout << "==============================================" << std::endl;
    std::cout << "Triangle Counting (GPU) " << std::endl;
    std::cout << "Node Ordering - Degree Based ID " << std::endl;
    std::cout << "Input Graph File : " << input_data_graph_file << std::endl;
    std::cout << "Exact Triangle Count : " << exact_count << std::endl;    
    std::cout << "Input File Reading Time : " << input_read_time << " seconds" << std::endl;
    std::cout << "Transformation Time : " << transformation_time << " seconds" << std::endl;
    std::cout << "Counting Time : " << counting_time << " seconds" << std::endl;
    std::cout << "Total Time : " <<  total_time <<  " seconds" << std::endl;
    std::cout << "==============================================" << std::endl;
    
}*/



int main(int argc, char** argv){

    std::string input_data_graph_file = argv[1];

    clock_t read_begin_clock = clock();
    Graph* graph = new Graph();
    //graph->loadGraphFromFile(input_data_graph_file);
    graph->loadGraphFromFileForBothDirectionEdges(input_data_graph_file);
    double input_read_time = (double(clock() - read_begin_clock)) / CLOCKS_PER_SEC;    


    clock_t transformation_begin_clock = clock();
    Graph* augmented_graph = new Graph();
    //graph->transformToAugmentedGraph(augmented_graph);
    graph->transformToAugmentedGraphForBothDirectionalEdges(augmented_graph);
    double transformation_time = (double(clock() - transformation_begin_clock)) / CLOCKS_PER_SEC;

    

    clock_t counting_begin_clock = clock();
    long long exact_count = CountingAlgorithm::db_count_triangle_warp_parallel(augmented_graph);  
    double counting_time = (double(clock() - counting_begin_clock)) / CLOCKS_PER_SEC;
    double total_time = input_read_time + transformation_time + counting_time;
   

    std::cout << "==============================================" << std::endl;
    std::cout << "Triangle Counting (GPU- Warp Per Vertex) " << std::endl;
    std::cout << "Node Ordering - Degree Based ID " << std::endl;
    std::cout << "Input Graph File : " << input_data_graph_file << std::endl;
    std::cout << "Exact Triangle Count : " << exact_count << std::endl;    
    std::cout << "Input File Reading Time : " << input_read_time << " seconds" << std::endl;
    std::cout << "Transformation Time : " << transformation_time << " seconds" << std::endl;
    std::cout << "Counting Time : " << counting_time << " seconds" << std::endl;
    std::cout << "Total Time : " <<  total_time <<  " seconds" << std::endl;
    std::cout << "==============================================" << std::endl;
    
}



