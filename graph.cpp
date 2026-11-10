

#include "graph.h"
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>

void Graph::loadGraphFromFileForBothDirectionEdges(const std::string &file_path)
{
    std::ifstream infile(file_path);

    if (!infile.is_open())
    {
        std::cout << "Can not open the graph file " << file_path << " ." << std::endl;
        exit(-1);
    }

    char type;
    std::string input_line;
    ui label = 0;

    std::cout << "Reading File............ " << std::endl;

    ui line_count = 0, count = 0, comment_line_count = 4;

    while (std::getline(infile, input_line))
    {

        if (input_line.rfind("#", 0) == 0)
        {
            line_count++;

            if (input_line.rfind("# Nodes", 0) == 0)
            {

                std::stringstream ss(input_line);
                std::string token;
                int count = 0;

                while (!ss.eof())
                {
                    std::getline(ss, token, ' ');

                    if (!(token.rfind("#", 0) == 0 || token.rfind("Nodes:", 0) == 0 || token.rfind("Edges:", 0) == 0))
                    {

                        if (count == 0)
                        {

                            vertices_count = stoi(token);
                            std::cout << "Vertex Count : " << vertices_count << std::endl;
                            degrees = new ui[vertices_count];
                            std::fill(degrees, degrees + vertices_count, 0);
                            count = 1;
                        }
                        else
                        {
                            count = 0;
                        }
                        std::cout << "Vertices Count : " << vertices_count << " Edges Count : " << edges_count << std::endl;
                    }
                }
            }
        }

        if (line_count >= comment_line_count){
            break;
        }
    }

    VertexID begin, end;

    edges_count = 0;

    while (infile >> begin)
    {
        infile >> end;

        if (begin != end && begin < vertices_count && end < vertices_count)
        {
            degrees[begin] += 1;
            edges_count++;
        }
    }

    infile.close();

    std::ifstream input_file(file_path);

    offsets = new ui[vertices_count + 1];
    offsets[0] = 0;

    neighbors = new VertexID[edges_count];

    std::cout << "Initialization Finished" << std::endl;

    std::vector<ui> neighbors_offset(vertices_count, 0); // used for adjust neighbors with offset

    for (ui id = 0; id < vertices_count; id++)
    {
        offsets[id + 1] = offsets[id] + degrees[id];
    }

    line_count = 0;

    while (std::getline(input_file, input_line))
    {
        line_count++;
        if (line_count >= comment_line_count)
        {
            break;
        }
    }

    while (input_file >> begin)
    {
        input_file >> end;

        line_count++;
        if (begin >= vertices_count || end >= vertices_count || begin == end)
        {
            continue;
        }

        ui offset = offsets[begin] + neighbors_offset[begin];
        neighbors[offset] = end;

        neighbors_offset[begin] += 1;
    }

    input_file.close();

    for (ui i = 0; i < vertices_count; ++i)
    {
        std::sort(neighbors + offsets[i], neighbors + offsets[i + 1]);
    }
}



bool Graph::is_smaller(VertexID u, VertexID v){


    if(u == v){
        return false;
    }

    if(main_degrees[u] < main_degrees[v]){
        return true;
    }

    if((main_degrees[u] == main_degrees[v]) && (u < v)){
        return true;
    }

    return false;
}


void Graph::transformToAugmentedGraph(Graph* augmented_graph){

    augmented_graph->vertices_count = vertices_count;
    augmented_graph->edges_count = edges_count / 2;
    augmented_graph->offsets = new ui[vertices_count + 1];
    augmented_graph->neighbors = new VertexID[edges_count];

    augmented_graph->offsets[0] = 0;

    augmented_graph->main_degrees = new ui[augmented_graph->vertices_count];
    augmented_graph->degrees = new ui[augmented_graph->vertices_count];
    std::fill(augmented_graph->degrees, augmented_graph->degrees + augmented_graph->vertices_count, 0);

    for (ui i = 0; i < vertices_count; i++)
    {
        for (ui j = offsets[i]; j < offsets[i + 1]; j++)
        {
            if ((degrees[neighbors[j]] > degrees[i]) || ((degrees[neighbors[j]] == degrees[i]) && neighbors[j] > i))
            {
                augmented_graph->degrees[i] += 1;
            }
        }
    }


    for (ui i = 0; i < vertices_count; i++)
    {
        augmented_graph->offsets[i + 1] = augmented_graph->offsets[i] + augmented_graph->degrees[i];
        augmented_graph->main_degrees[i] = degrees[i];
    }

    std::vector<ui> neighbors_offset(vertices_count, 0);
    ui offset;

    for (ui i = 0; i < vertices_count; i++)
    {
        for (ui j = offsets[i]; j < offsets[i + 1]; j++)
        {

            if ((degrees[neighbors[j]] > degrees[i]) || ((degrees[neighbors[j]] == degrees[i]) && neighbors[j] > i))
            {

                offset = augmented_graph->offsets[i] + neighbors_offset[i];
                augmented_graph->neighbors[offset] = neighbors[j];
                neighbors_offset[i] += 1;
            }
        }
    }   
}
