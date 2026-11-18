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



/*int main(int argc, char** argv){

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
    
}*/

//Parameterized Parallelism
/*int main(int argc, char** argv){

    std::string input_data_graph_file = argv[1];
    std::string is_both_direction = argv[2];
    std::string thread_cnt_str = argv[3];

    int thread_cnt = std::stoi(thread_cnt_str);
    long long exact_count = 0;

    Graph* graph = new Graph();
    Graph* augmented_graph = new Graph();
    double input_read_time, transformation_time, counting_time, total_time;

    if(is_both_direction == "1"){

        clock_t read_begin_clock = clock();
        graph->loadGraphFromFileForBothDirectionEdges(input_data_graph_file);
        input_read_time = (double(clock() - read_begin_clock)) / CLOCKS_PER_SEC; 

        clock_t transformation_begin_clock = clock();
        graph->transformToAugmentedGraphForBothDirectionalEdges(augmented_graph);
        transformation_time = (double(clock() - transformation_begin_clock)) / CLOCKS_PER_SEC;

    } else{

        clock_t read_begin_clock = clock();
        graph->loadGraphFromFile(input_data_graph_file);
        input_read_time = (double(clock() - read_begin_clock)) / CLOCKS_PER_SEC; 

        clock_t transformation_begin_clock = clock();
        graph->transformToAugmentedGraph(augmented_graph);
        transformation_time = (double(clock() - transformation_begin_clock)) / CLOCKS_PER_SEC;
    }

    threads_per_block = thread_cnt;

    std::cout << "==============================================" << std::endl;
    std::cout << "Input Graph File : " << input_data_graph_file << std::endl;
    std::cout << "Thread Count Per Block : " << threads_per_block << std::endl;
    std::cout << "Input File Reading Time : " << input_read_time << " seconds" << std::endl;
    std::cout << "Transformation Time : " << transformation_time << " seconds" << std::endl;
    std::cout << "==============================================" << std::endl;

    clock_t counting_begin_clock = clock();
    exact_count = CountingAlgorithm::db_count_triangle_parallel(augmented_graph);  
    counting_time = (double(clock() - counting_begin_clock)) / CLOCKS_PER_SEC;
    total_time = input_read_time + transformation_time + counting_time;
   
    std::cout << "==============================================" << std::endl;
    std::cout << "Triangle Counting (Thread Per Vertex) " << std::endl;
    std::cout << "Node Ordering - Degree Based ID " << std::endl;    
    std::cout << "Exact Triangle Count : " << exact_count << std::endl;    
    std::cout << "Counting Time : " << counting_time << " seconds" << std::endl;
    std::cout << "Total Time : " <<  total_time <<  " seconds" << std::endl;
    std::cout << "==============================================" << std::endl;


    counting_begin_clock = clock();
    exact_count = CountingAlgorithm::db_count_triangle_warp_parallel(augmented_graph);  
    counting_time = (double(clock() - counting_begin_clock)) / CLOCKS_PER_SEC;
    total_time = input_read_time + transformation_time + counting_time;   

    std::cout << "==============================================" << std::endl;
    std::cout << "Triangle Counting (GPU- Warp Per Vertex) " << std::endl;
    std::cout << "Node Ordering - Degree Based ID " << std::endl;
    std::cout << "Exact Triangle Count : " << exact_count << std::endl;    
    std::cout << "Counting Time : " << counting_time << " seconds" << std::endl;
    std::cout << "Total Time : " <<  total_time <<  " seconds" << std::endl;
    std::cout << "==============================================" << std::endl;

    VertexID* src_vtx_array = new VertexID[augmented_graph->edges_count];
    augmented_graph->fillSrcVertexArray(src_vtx_array);

    counting_begin_clock = clock();
    exact_count = CountingAlgorithm::db_count_triangle_ei_parallel(augmented_graph, src_vtx_array);  
    counting_time = (double(clock() - counting_begin_clock)) / CLOCKS_PER_SEC;
    total_time = input_read_time + transformation_time + counting_time;
   

    std::cout << "==============================================" << std::endl;
    std::cout << "Triangle Counting (GPU - Thread Per Edge) " << std::endl;
    std::cout << "Node Ordering - Degree Based ID " << std::endl;
    std::cout << "Exact Triangle Count : " << exact_count << std::endl;    
    std::cout << "Counting Time : " << counting_time << " seconds" << std::endl;
    std::cout << "Total Time : " <<  total_time <<  " seconds" << std::endl;
    std::cout << "==============================================" << std::endl;
    
}*/


// Thread Per Edge
/*int main(int argc, char** argv){

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

    VertexID* src_vtx_array = new VertexID[augmented_graph->edges_count];
    augmented_graph->fillSrcVertexArray(src_vtx_array);

    double transformation_time = (double(clock() - transformation_begin_clock)) / CLOCKS_PER_SEC;
    

    clock_t counting_begin_clock = clock();
    long long exact_count = CountingAlgorithm::db_count_triangle_ei_parallel(augmented_graph, src_vtx_array);  
    double counting_time = (double(clock() - counting_begin_clock)) / CLOCKS_PER_SEC;
    double total_time = input_read_time + transformation_time + counting_time;
   

    std::cout << "==============================================" << std::endl;
    std::cout << "Triangle Counting (GPU - Thread Per Edge) " << std::endl;
    std::cout << "Node Ordering - Degree Based ID " << std::endl;
    std::cout << "Input Graph File : " << input_data_graph_file << std::endl;
    std::cout << "Exact Triangle Count : " << exact_count << std::endl;    
    std::cout << "Input File Reading Time : " << input_read_time << " seconds" << std::endl;
    std::cout << "Transformation Time : " << transformation_time << " seconds" << std::endl;
    std::cout << "Counting Time : " << counting_time << " seconds" << std::endl;
    std::cout << "Total Time : " <<  total_time <<  " seconds" << std::endl;
    std::cout << "==============================================" << std::endl;
    
}*/


// Square Count - Thread Per Vertex
/*int main(int argc, char** argv){

    std::string input_data_graph_file = argv[1];
    std::string is_both_direction = argv[2];
    std::string thread_cnt_str = argv[3];

    int thread_cnt = std::stoi(thread_cnt_str);
    long long exact_count = 0;

    Graph* graph = new Graph();
    Graph* augmented_graph = new Graph();
    double input_read_time, transformation_time, counting_time, total_time;

    if(is_both_direction == "1"){

        clock_t read_begin_clock = clock();
        graph->loadGraphFromFileForBothDirectionEdges(input_data_graph_file);
        input_read_time = (double(clock() - read_begin_clock)) / CLOCKS_PER_SEC; 

        clock_t transformation_begin_clock = clock();
        graph->transformToAugmentedGraphForBothDirectionalEdges(augmented_graph);
        transformation_time = (double(clock() - transformation_begin_clock)) / CLOCKS_PER_SEC;

    } else{

        clock_t read_begin_clock = clock();
        graph->loadGraphFromFile(input_data_graph_file);
        input_read_time = (double(clock() - read_begin_clock)) / CLOCKS_PER_SEC; 

        clock_t transformation_begin_clock = clock();
        graph->transformToAugmentedGraph(augmented_graph);
        transformation_time = (double(clock() - transformation_begin_clock)) / CLOCKS_PER_SEC;
    }

    threads_per_block = thread_cnt;    

    clock_t counting_begin_clock = clock();
    exact_count = CountingAlgorithm::db_count_square_parallel(graph, augmented_graph);  
    counting_time = (double(clock() - counting_begin_clock)) / CLOCKS_PER_SEC;
    total_time = input_read_time + transformation_time + counting_time;
   

    std::cout << "==============================================" << std::endl;
    std::cout << "Square Counting (GPU - Thread Per Vertex) " << std::endl;
    std::cout << "Node Ordering - Degree Based ID " << std::endl;
    std::cout << "Input Graph File : " << input_data_graph_file << std::endl;
    std::cout << "Exact Square Count : " << exact_count << std::endl;    
    std::cout << "Input File Reading Time : " << input_read_time << " seconds" << std::endl;
    std::cout << "Transformation Time : " << transformation_time << " seconds" << std::endl;
    std::cout << "Counting Time : " << counting_time << " seconds" << std::endl;
    std::cout << "Total Time : " <<  total_time <<  " seconds" << std::endl;
    std::cout << "==============================================" << std::endl;
    
}*/


// Square Count - Thread Per Edge
int main(int argc, char** argv){

    std::string input_data_graph_file = argv[1];
    std::string is_both_direction = argv[2];
    std::string thread_cnt_str = argv[3];

    int thread_cnt = std::stoi(thread_cnt_str);
    long long exact_count = 0;

    Graph* graph = new Graph();
    double input_read_time, counting_time, total_time;

    if(is_both_direction == "1"){

        clock_t read_begin_clock = clock();
        graph->loadGraphFromFileForBothDirectionEdges(input_data_graph_file);
        input_read_time = (double(clock() - read_begin_clock)) / CLOCKS_PER_SEC; 

    } else{

        clock_t read_begin_clock = clock();
        graph->loadGraphFromFile(input_data_graph_file);
        input_read_time = (double(clock() - read_begin_clock)) / CLOCKS_PER_SEC; 
    }

    VertexID* src_vtx_array = new VertexID[(graph->edges_count) * 2];
    graph->fillSrcVertexArray(src_vtx_array);

    threads_per_block = thread_cnt;    

    clock_t counting_begin_clock = clock();
    exact_count = CountingAlgorithm::db_count_square_ei_parallel(graph, src_vtx_array);  
    counting_time = (double(clock() - counting_begin_clock)) / CLOCKS_PER_SEC;
    total_time = input_read_time + counting_time;
   

    std::cout << "==============================================" << std::endl;
    std::cout << "Square Counting (GPU - Thread Per Vertex) " << std::endl;
    std::cout << "Node Ordering - Degree Based ID " << std::endl;
    std::cout << "Input Graph File : " << input_data_graph_file << std::endl;
    std::cout << "Exact Square Count : " << exact_count << std::endl;    
    std::cout << "Input File Reading Time : " << input_read_time << " seconds" << std::endl;
    std::cout << "Counting Time : " << counting_time << " seconds" << std::endl;
    std::cout << "Total Time : " <<  total_time <<  " seconds" << std::endl;
    std::cout << "==============================================" << std::endl;
    
}





