#include <stdio.h>
#include <stdlib.h>

// graph.h

#ifndef GRAPH_H
#define GRAPH_H

// Graph Class

typedef struct
{
    int id;
    int bonus;
    int penalty;
} Node;

typedef struct
{
    int node1;
    int node2;
    int weight;
} Edge;

typedef struct
{
    int size;
    Node *nodes;
    Edge **edges;
} Graph;

Graph *create_graph(int size, int **distance_matrix, int **node_data)
{
    int edge_size = 3 * sizeof(int);
    Graph *graph = (Graph *)calloc(1, sizeof(Graph));
    graph->size = size;
    graph->nodes = (Node *)calloc(size, sizeof(Node));
    graph->edges = (Edge **)calloc(size, sizeof(Edge *));

    // create nodes
    for (int i = 0; i < size; i++)
    {
        graph->edges[i] = (Edge *)calloc(size, edge_size);

        Node node = {i, node_data[i][1], node_data[i][2]};
        graph->nodes[i] = node;
    }

    // create edges
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            int weight = distance_matrix[i][j];
            if (weight > 0)
            {
                graph->edges[i][j].node1 = i;
                graph->edges[i][j].node2 = j;
                graph->edges[i][j].weight = weight;
            }
            else
            {
                graph->edges[i][j].node1 = i;
                graph->edges[i][j].node2 = j;
                graph->edges[i][j].weight = graph->nodes[i].penalty;
            }
        }
    }

    return graph;
}

#endif // GRAPH_H
