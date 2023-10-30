#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "globals.h"
#include "aux_functions.h"

#ifndef FITNESS_H
#define FITNESS_H

route_cost_calls = 0;
// neighboor_1_calls = 0;
// neighboor_2_calls = 0;
// neighboor_3_calls = 0;
// neighboor_4_calls = 0;
// neighboor_5_calls = 0;
// ref_1_calls = 0;
// ref_2_calls = 0;
// ref_3_calls = 0;

int route_cost(int *route, Graph *graph, int route_len)
{
    route_cost_calls++;

    int cost = 0;

    for (int i = 0; i < route_len - 1; i++)
    {
        int node1 = route[i];
        int node2 = route[i + 1];
        cost += graph->edges[node1][node2].weight;
    }
    int last_node_id = route[route_len - 1];
    cost += graph->edges[last_node_id][0].weight;

    for (int i = 0; i < graph->size; i++)
    {
        if (!contains(route, route_len, graph->nodes[i].id))
        {
            cost += graph->nodes[i].penalty;
        }
    }

    return cost;
}

float calculate_fitness(int *route, Graph *graph, int route_len)
{
    int cost = route_cost(route, graph, route_len);
    return 1.0 / (cost + 1);
}

#endif