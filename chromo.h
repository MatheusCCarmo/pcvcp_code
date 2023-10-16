#include <stdio.h>
#include <stdlib.h>
#include "fitness.h"
#include "graph.h"

#ifndef CHROMO_H
#define CHROMO_H

typedef struct
{
    int *route;
    int route_len;
    float fitness;
} Chromo;

Chromo *create_chromo(int *route, Graph *graph)
{
    int chromo_size = (graph->size * sizeof(int)) + sizeof(int) + sizeof(float);
    Chromo *chromo = (Chromo *)malloc(chromo_size);

    chromo->route = route;

    int route_len = get_route_length(route);

    chromo->route_len = route_len;

    chromo->fitness = calculate_fitness(route, graph, route_len);

    return chromo;
}

void destroy_chromo(Chromo *chromo)
{
    free(chromo->route);
    free(chromo);
}

#endif