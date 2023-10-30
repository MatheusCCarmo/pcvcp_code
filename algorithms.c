#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <stdbool.h>
#include "graph.h"
#include "chromo.h"
#include "aux_functions.h"
#include "algorithms.h"

#define MAX_COST_CALLS 10000000

double RECOMBINATION_RATE;
double MUTATION_RATE;
double POPULATION_SIZE;
double TOURNAMENT_SIZE;

typedef int *(*FunctionPointer)();

int compare_tuples(const void *a, const void *b)
{

    int *tuple1 = *(int **)a;

    int *tuple2 = *(int **)b;

    if (tuple1[1] < tuple2[1])
    {
        return 1;
    }
    else if (tuple1[1] > tuple2[1])
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

int compare_chromos(const void *a, const void *b)
{
    Chromo *chromoA = *(Chromo **)a;
    Chromo *chromoB = *(Chromo **)b;

    if (chromoA->fitness < chromoB->fitness)
    {
        return 1;
    }
    else if (chromoA->fitness > chromoB->fitness)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

int calculate_bonus_collected(int *route, Graph *graph)
{
    int bonus = 0;
    int id;
    int i = 0;
    id = route[i];

    while (id != -1)
    {
        bonus += graph->nodes[id].bonus;
        id = route[++i];
    }

    return bonus;
}

int *generate_random_route(Graph *graph, int quota)
{

    int bonus_collected = 0;
    int *route;
    route = (int *)calloc(graph->size, sizeof(int));
    int route_index = 0;
    route[route_index++] = 0;
    while (bonus_collected < quota && route_index < graph->size)
    {
        int random_i = rand() % graph->size;
        if (random_i != 0 && !contains(route, route_index, random_i))
        {
            bonus_collected += graph->nodes[random_i].bonus;
            route[route_index++] = random_i;
        }
    }

    route[route_index] = -1;
    return route;
}

int *swap_2(int *route, int i, int k)
{
    int route_len = get_route_length(route);

    int *new_route = (int *)malloc((route_len + 1) * sizeof(int));

    int a, b, j;

    for (j = 0; j < i; ++j)
    {
        new_route[j] = route[j];
    }

    for (a = i, b = k; a <= k; ++a, --b)
    {
        new_route[a] = route[b];
    }

    for (j = k + 1; j < route_len; ++j)
    {
        new_route[j] = route[j];
    }

    new_route[route_len] = -1;

    return new_route;
}

int *swap_2_opt(int *route, int quota, Graph *graph)
{
    // ref_1_calls++;
    if (route_cost_calls >= MAX_COST_CALLS)
    {
        return route;
    }

    int route_len = get_route_length(route);
    int *new_route = (int *)malloc(route_len * sizeof(int));

    int *best_route = (int *)malloc(route_len * sizeof(int));

    int new_route_len, new_cost;

    bool improved = true;

    best_route = route;
    int best_route_len = get_route_length(best_route);

    int best_cost = route_cost(best_route, graph, best_route_len);

    int no_improvement_count = 1000;
    int counter = 0;

    while (improved && counter < no_improvement_count)
    {
        counter++;
        improved = false;

        for (int i = 1; i < best_route_len - 1; i++)
        {
            for (int j = i + 1; j < best_route_len; j++)
            {
                new_route = swap_2(best_route, i, j);
                new_route_len = get_route_length(new_route);
                new_cost = route_cost(new_route, graph, new_route_len);

                if (new_cost < best_cost)
                {
                    best_cost = new_cost;
                    best_route = new_route;
                    best_route_len = get_route_length(best_route);
                    improved = true;
                    counter = 0;
                }
                if (route_cost_calls > MAX_COST_CALLS)
                {
                    return best_route;
                }
            }
        }
    }
    return best_route;
}

// GENETIC

Chromo *mutation(Chromo *chromo, Graph *graph)
{

    int *new_route = (int *)malloc(graph->size * sizeof(int));
    int chromo_size = (graph->size * sizeof(int)) + sizeof(int) + sizeof(float);
    Chromo *new_chromo = (Chromo *)malloc(chromo_size);

    int i = 1 + rand() % (chromo->route_len - 2);
    int k = i + 1 + rand() % (chromo->route_len - i - 1);

    new_route = swap_2(chromo->route, i, k);

    new_chromo = create_chromo(new_route, graph);

    return new_chromo;
}

Chromo *tournament_selection(Chromo **population, int population_size, int tournament_size)
{
    Chromo *tournament_candidates[tournament_size];
    int random_i;

    for (int i = 0; i < tournament_size; i++)
    {
        random_i = rand() % population_size;
        tournament_candidates[i] = population[random_i];
    }

    Chromo *best_candidate = tournament_candidates[0];
    for (int i = 1; i < tournament_size; i++)
    {
        if (tournament_candidates[i]->fitness > best_candidate->fitness)
        {
            best_candidate = tournament_candidates[i];
        }
    }

    return best_candidate;
}

Chromo **init_population(Graph *graph, int population_size, int quota)
{
    int chromo_size = (graph->size * sizeof(int)) + sizeof(int) + sizeof(float);
    Chromo **population = (Chromo **)malloc(population_size * chromo_size);

    for (int i = 0; i < population_size; i++)
    {
        Chromo *chromo;

        int *route;
        route = (int *)calloc(graph->size, sizeof(int));

        route = generate_random_route(graph, quota);

        chromo = create_chromo(route, graph);

        population[i] = chromo;
    }

    return population;
}

void crossover(Chromo *parent_1, Chromo *parent_2, Chromo **child_1, Chromo **child_2, Graph *graph)
{
    int len_min = fmin(parent_1->route_len, parent_2->route_len);

    int crossover_point_begin = 1 + rand() % (len_min - 2);

    int crossover_point_end = crossover_point_begin + 1 + rand() % (len_min - crossover_point_begin - 1);

    int *route_1 = (int *)malloc(graph->size * sizeof(int));

    int *route_2 = (int *)malloc(graph->size * sizeof(int));

    route_1[0] = 0;
    route_2[0] = 0;
    int route_1_len = 1;
    int route_2_len = 1;

    for (int i = 1; i < crossover_point_begin; i++)
    {
        route_1[route_1_len++] = parent_1->route[i];
        route_2[route_2_len++] = parent_2->route[i];
    }
    for (int i = crossover_point_begin; i < crossover_point_end; i++)
    {
        if (!contains(route_1, route_1_len, parent_2->route[i]))
        {
            route_1[route_1_len++] = parent_2->route[i];
        }
        if (!contains(route_2, route_2_len, parent_1->route[i]))
        {
            route_2[route_2_len++] = parent_1->route[i];
        }
    }
    for (int i = 1; i < parent_1->route_len; i++)
    {
        if (!contains(route_1, route_1_len, parent_1->route[i]))
        {
            route_1[route_1_len++] = parent_1->route[i];
        }
    }
    for (int i = 1; i < parent_2->route_len; i++)
    {
        if (!contains(route_2, route_2_len, parent_2->route[i]))
        {
            route_2[route_2_len++] = parent_2->route[i];
        }
    }

    route_1[route_1_len] = -1;
    route_2[route_2_len] = -1;

    *child_1 = create_chromo(route_1, graph);
    *child_2 = create_chromo(route_2, graph);
}

int *genetic_algorithm(Graph *graph, int quota)
{
    int chromo_size = (graph->size * sizeof(int)) + sizeof(int) + sizeof(float);
    int population_size = POPULATION_SIZE;
    int tournament_size = TOURNAMENT_SIZE;
    int *route = (int *)malloc(graph->size * sizeof(int));
    int children_len = 0;
    Chromo **new_population;
    Chromo **population;
    Chromo **children;
    Chromo *parent_1;
    Chromo *parent_2;
    Chromo *child_1;
    Chromo *child_2;
    Chromo *child_3;
    Chromo *selected_child;

    population = (Chromo **)calloc(population_size, chromo_size);
    parent_1 = (Chromo *)malloc(chromo_size);
    parent_2 = (Chromo *)malloc(chromo_size);
    child_1 = (Chromo *)malloc(chromo_size);
    child_2 = (Chromo *)malloc(chromo_size);
    child_3 = (Chromo *)malloc(chromo_size);
    selected_child = (Chromo *)malloc(chromo_size);
    population = init_population(graph, population_size, quota);

    int children_size = 3 * population_size;
    int new_population_size = children_size + population_size;

    children = (Chromo **)calloc(children_size, chromo_size);
    new_population = (Chromo **)calloc(new_population_size, chromo_size);
    while (route_cost_calls < MAX_COST_CALLS)
    {
        children_len = 0;
        for (int i = 0; i < population_size; i++)
        {
            parent_1 = tournament_selection(population, population_size, tournament_size);
            parent_2 = tournament_selection(population, population_size, tournament_size);

            crossover(parent_1, parent_2, &child_1, &child_2, graph);

            children[children_len++] = child_1;
            children[children_len++] = child_2;

            int random_i = rand() % children_len;
            selected_child = children[random_i];

            child_3 = mutation(selected_child, graph);

            children[children_len++] = child_3;
        }
        for (int i = 0; i < population_size; i++)
        {
            new_population[i] = population[i];
        }

        for (int i = 0; i < children_size; i++)
        {
            new_population[i + population_size] = children[i];
        }

        qsort(new_population, new_population_size, sizeof(Chromo **), compare_chromos);

        for (int i = 0; i < population_size; i++)
        {
            population[i] = new_population[i];
        }
    }

    route = population[0]->route;

    return route;
}

// MEMETIC

// VNS

// Remove vértice com maior economia:
int *neighboor_1(int *route, Graph *graph)
{
    // printf("\n1");
    // neighboor_1_calls++;
    int *new_route = (int *)malloc(graph->size * sizeof(int));

    int route_len = get_route_length(route);

    memcpy(new_route, route, (route_len + 1) * sizeof(int));

    int tuple_size = 2 * sizeof(int);
    int i, j, s;
    int best_economy_id;
    int economy_value;
    int economies_len = 0;

    Edge edge;
    Edge edge1;
    Edge edge2;

    int *economy_tuple;

    int **economies = (int **)malloc(graph->size * sizeof(int *));

    for (int i = 0; i < graph->size; i++)
    {
        economies[i] = (int *)malloc(tuple_size);
    }

    route_len = get_route_length(new_route);

    for (int r = 1; r < route_len; r++)
    {
        economy_tuple = (int *)malloc(tuple_size);

        if (r == route_len - 1)
        {
            i = new_route[r - 1];
            j = new_route[r];
            s = new_route[0];
        }
        else
        {
            i = new_route[r - 1];
            j = new_route[r];
            s = new_route[r + 1];
        }

        edge1 = graph->edges[i][j];
        edge2 = graph->edges[j][s];
        edge = graph->edges[i][s];
        economy_value = edge1.weight + edge2.weight - edge.weight - graph->nodes[j].penalty;

        economy_tuple[0] = j;
        economy_tuple[1] = economy_value;
        economies[economies_len++] = economy_tuple;
    }

    qsort(economies, economies_len, sizeof(int *), compare_tuples);

    best_economy_id = economies[0][0];

    remove_by_value(new_route, route_len--, best_economy_id);

    return new_route;
}

// Insere vértice com maior economia:
int *neighboor_2(int *route, Graph *graph)
{
    // printf("\n2");
    // neighboor_2_calls++;
    int *new_route = (int *)malloc(graph->size * sizeof(int));

    int route_len = get_route_length(route);

    memcpy(new_route, route, (route_len + 1) * sizeof(int));

    int tuple_size = 3 * sizeof(int);
    int i, j, id;
    int best_economy_id;
    int best_economy_insert_index;
    int economy_value;
    int economies_len = 0;

    Edge edge;
    Edge edge1;
    Edge edge2;

    int *economy_tuple;
    int **economies = (int **)malloc(graph->size * graph->size * sizeof(int *));

    for (int i = 0; i < graph->size; i++)
    {
        economies[i] = (int *)malloc(tuple_size);
    }

    route_len = get_route_length(new_route);

    for (int k = 0; k < graph->size; k++)
    {

        id = graph->nodes[k].id;
        if (!contains(new_route, route_len, id))
        {
            for (int r = 1; r < route_len; r++)
            {
                economy_tuple = (int *)malloc(tuple_size);

                i = new_route[r - 1];
                j = new_route[r];
                edge = graph->edges[i][j];
                edge1 = graph->edges[i][id];
                edge2 = graph->edges[id][j];
                economy_value = edge.weight + graph->nodes[id].penalty -
                                edge1.weight - edge2.weight;

                economy_tuple[0] = id;
                economy_tuple[1] = economy_value;
                economy_tuple[2] = r;

                economies[economies_len++] = economy_tuple;
            }
        }
    }

    qsort(economies, economies_len, sizeof(int *), compare_tuples);

    best_economy_id = economies[0][0];

    best_economy_insert_index = economies[0][2];

    if (best_economy_insert_index == 0)
    {
        insert_on_index(new_route, route_len, route_len, best_economy_id);
        route_len++;
    }
    else
    {
        insert_on_index(new_route, route_len, best_economy_insert_index, best_economy_id);
        route_len++;
    }

    return new_route;
}

// Troca dois vertices

int *neighboor_3(int *route, Graph *graph)
{
    // printf("\n3");
    // neighboor_3_calls++;

    int *new_route = (int *)malloc(graph->size * sizeof(int));
    int route_len = get_route_length(route);

    memcpy(new_route, route, (route_len + 1) * sizeof(int));

    int i = 1 + rand() % (route_len - 2);
    int k = i + 1 + rand() % (route_len - i - 1);

    new_route = swap_2(new_route, i, k);

    return new_route;
}

// Remove vertice aleatorio:
int *neighboor_4(int *route, Graph *graph)
{
    // printf("\n4");
    // neighboor_4_calls++;

    int route_len = get_route_length(route);

    int *new_route = (int *)malloc(graph->size * sizeof(int));

    memcpy(new_route, route, (route_len + 1) * sizeof(int));

    int index = 1 + rand() % (route_len - 1);

    remove_by_index(new_route, route_len--, index);

    return new_route;
}

// Insere vertice aleatorio:
int *neighboor_5(int *route, Graph *graph)
{
    // printf("\n5");
    // neighboor_5_calls++;

    int *new_route = (int *)malloc(graph->size * sizeof(int));

    int route_len = get_route_length(route);

    memcpy(new_route, route, (route_len + 1) * sizeof(int));

    int index = 1 + rand() % (route_len - 1);
    int node = 1 + rand() % (graph->size - 1);

    while (contains(new_route, route_len, node))
    {
        node = 1 + rand() % (graph->size - 1);
    }

    insert_on_index(new_route, route_len, index, node);

    return new_route;
}

void add_step(int *route, int quota, Graph *graph)
{

    int best_economy_value = INT_MIN;
    int economy_value, i, j, k_best_economy, r_best_economy;
    Edge edge, edge1, edge2;

    int route_len = get_route_length(route);

    int bonus_collected = calculate_bonus_collected(route, graph);

    do
    {
        best_economy_value = INT_MIN;
        k_best_economy = 0;
        for (int k = 0; k < graph->size; k++)
        {
            if (!contains(route, route_len, graph->nodes[k].id))
            {
                for (int r = 1; r < route_len; r++)
                {
                    i = route[r - 1];
                    j = route[r];

                    edge = graph->edges[i][j];
                    edge1 = graph->edges[i][k];
                    edge2 = graph->edges[k][j];
                    economy_value = edge.weight + graph->nodes[k].penalty -
                                    edge1.weight - edge2.weight;

                    if (economy_value > best_economy_value)
                    {
                        best_economy_value = economy_value;
                        k_best_economy = k;
                        r_best_economy = r;
                    }
                }
            }
        }
        if (best_economy_value > 0 || bonus_collected < quota)
        {
            insert_on_index(route, route_len++, r_best_economy, graph->nodes[k_best_economy].id);
        }
        bonus_collected = calculate_bonus_collected(route, graph);
    } while (bonus_collected < quota || best_economy_value > 0);
}

void drop_step(int *route, int quota, Graph *graph)
{

    int best_economy_value = INT_MIN;
    int economy_value, i, j, k, k_best_economy;
    Edge edge, edge1, edge2;

    int route_len = get_route_length(route);

    int bonus_collected = calculate_bonus_collected(route, graph);

    do
    {

        best_economy_value = INT_MIN;
        k_best_economy = 0;

        for (int r = 1; r < route_len; r++)
        {

            if (r == (route_len - 1))
            {
                i = route[r - 1];
                j = route[0];
                k = route[r];
            }
            else
            {
                i = route[r - 1];
                j = route[r + 1];
                k = route[r];
            }

            edge = graph->edges[i][j];
            edge1 = graph->edges[i][k];
            edge2 = graph->edges[k][j];
            economy_value = edge1.weight + edge2.weight - edge.weight - graph->nodes[j].penalty;

            if (economy_value > best_economy_value && (bonus_collected - graph->nodes[k].bonus) > quota)
            {
                best_economy_value = economy_value;
                k_best_economy = k;
            }
        }

        if (best_economy_value > 0)
        {
            remove_by_value(route, route_len--, graph->nodes[k_best_economy].id);
        }
        bonus_collected = calculate_bonus_collected(route, graph);
    } while (best_economy_value > 0);
}

int *add_drop(int *route, int quota, Graph *graph)
{
    // ref_2_calls++;
    int *new_route = (int *)malloc(graph->size * sizeof(int));

    int route_len = get_route_length(route);

    memcpy(new_route, route, (route_len + 1) * sizeof(int));

    add_step(new_route, quota, graph);
    drop_step(new_route, quota, graph);

    return new_route;
}

int *drop_add(int *route, int quota, Graph *graph)
{
    // ref_3_calls++;
    int *new_route = (int *)malloc(graph->size * sizeof(int));

    int route_len = get_route_length(route);
    memcpy(new_route, route, (route_len + 1) * sizeof(int));

    add_step(new_route, quota, graph);
    drop_step(new_route, quota, graph);

    return new_route;
}

int *vns(int *route, int quota, Graph *graph)
{
    if (route_cost_calls >= MAX_COST_CALLS)
    {
        return route;
    }

    clock_t start_time, current_time;
    double no_improvement_time = 0.0;
    int max_time = 200;
    int *best_route = (int *)malloc(graph->size * sizeof(int));
    best_route = route;

    int route_len = get_route_length(best_route);
    int best_cost = route_cost(best_route, graph, route_len);
    int bonus_collected = 0, new_route_len = 0, new_cost = 0, k = 0;
    int *new_route = (int *)malloc(graph->size * sizeof(int));
    new_route_len = route_len;

    FunctionPointer neighborhood_functions[] = {neighboor_1, neighboor_2, neighboor_3, neighboor_4, neighboor_5};

    int neighboors_len = sizeof(neighborhood_functions) / sizeof(FunctionPointer);

    while (no_improvement_time < max_time)
    {
        start_time = clock();

        k = 0;

        while (k < neighboors_len)
        {

            new_route = neighborhood_functions[k](best_route, graph);

            new_route_len = get_route_length(new_route);

            bonus_collected = calculate_bonus_collected(new_route, graph);

            new_route = swap_2_opt(new_route, quota, graph);

            new_route_len = get_route_length(new_route);

            new_cost = route_cost(new_route, graph, new_route_len);

            if (new_cost < best_cost && bonus_collected >= quota)
            {
                best_cost = new_cost;
                best_route = new_route;
                route_len = get_route_length(best_route);
                k = 0;
            }
            else
            {
                free(new_route);
                k++;
            }

            current_time = clock();
            no_improvement_time = ((double)(current_time - start_time));
            if (route_cost_calls >= MAX_COST_CALLS)
            {
                return best_route;
            }
        }

        return best_route;
    }

    return best_route;
}

int *vnd(int *route, int quota, Graph *graph)
{
    if (route_cost_calls >= MAX_COST_CALLS)
    {
        return route;
    }
    int *best_route = route;

    int route_len = get_route_length(best_route);
    int best_cost = route_cost(best_route, graph, route_len);
    int new_route_len = 0, new_cost = 0, k = 0;
    int *new_route = (int *)malloc(graph->size * sizeof(int));

    FunctionPointer refinements_functions[] = {add_drop, swap_2_opt, drop_add};

    int refinements_len = sizeof(refinements_functions) / sizeof(FunctionPointer);

    k = 0;

    while (k < refinements_len)
    {

        new_route = refinements_functions[k](best_route, quota, graph);

        new_route_len = get_route_length(new_route);
        new_cost = route_cost(new_route, graph, new_route_len);

        if (new_cost < best_cost)
        {
            best_cost = new_cost;
            best_route = new_route;
            k = 0;
        }
        else
        {
            k++;
        }
        if (route_cost_calls >= MAX_COST_CALLS)
        {
            return best_route;
        }
    }

    return best_route;
}

int *vns_vnd(int *route, int quota, Graph *graph)
{
    if (route_cost_calls >= MAX_COST_CALLS)
    {
        return route;
    }
    clock_t start_time, current_time;
    double no_improvement_time = 0.0;
    int max_time = 200;
    int *best_route = route;

    int route_len = get_route_length(best_route);
    int best_cost = route_cost(best_route, graph, route_len);
    int bonus_collected = 0, new_route_len = 0, new_cost = 0, k = 0;
    int *new_route = (int *)malloc(graph->size * sizeof(int));

    FunctionPointer neighborhood_functions[] = {neighboor_1, neighboor_2, neighboor_3, neighboor_4, neighboor_5};

    int neighboors_len = sizeof(neighborhood_functions) / sizeof(FunctionPointer);

    while (no_improvement_time < max_time)
    {
        start_time = clock();

        k = 0;

        while (k < neighboors_len)
        {

            new_route = neighborhood_functions[k](best_route, graph); // s' <- N(s)

            bonus_collected = calculate_bonus_collected(new_route, graph);

            if (bonus_collected >= quota)
            {
                new_route = vnd(new_route, quota, graph); // s'' <- VND(s')
                new_route_len = get_route_length(new_route);
                new_cost = route_cost(new_route, graph, new_route_len);

                if (new_cost < best_cost)
                {
                    best_cost = new_cost;
                    best_route = new_route;
                    k = 0;
                }
                else
                {
                    free(new_route);
                    k++;
                }
            }
            else
            {
                k++;
            }
            current_time = clock();
            no_improvement_time = ((double)(current_time - start_time));
            if (route_cost_calls >= MAX_COST_CALLS)
            {
                return best_route;
            }
        }

        return best_route;
    }
    return best_route;
}

int *memetic_algorithm(Graph *graph, int quota)
{

    int chromo_size = (graph->size * sizeof(int)) + sizeof(int) + sizeof(float);
    int population_size = POPULATION_SIZE;
    int tournament_size = TOURNAMENT_SIZE;
    int *route = (int *)malloc(graph->size * sizeof(int));
    int *new_route = (int *)malloc(graph->size * sizeof(int));
    double random_number;

    int children_len = 0;
    Chromo **new_population;
    Chromo **population;
    Chromo **children;
    Chromo **crossover_children;
    Chromo *parent_1;
    Chromo *parent_2;
    Chromo *child_1;
    Chromo *child_2;
    Chromo *new_chromo;
    Chromo *selected_child;

    population = (Chromo **)calloc(population_size, chromo_size);
    parent_1 = (Chromo *)malloc(chromo_size);
    parent_2 = (Chromo *)malloc(chromo_size);
    child_1 = (Chromo *)malloc(chromo_size);
    child_2 = (Chromo *)malloc(chromo_size);
    new_chromo = (Chromo *)malloc(chromo_size);
    selected_child = (Chromo *)malloc(chromo_size);
    population = init_population(graph, population_size, quota);

    int children_size = 3 * population_size;
    int new_population_size = children_size + population_size;

    children = (Chromo **)calloc(children_size, chromo_size);
    crossover_children = (Chromo **)calloc(2, chromo_size);
    new_population = (Chromo **)calloc(new_population_size, chromo_size);
    while (route_cost_calls < MAX_COST_CALLS)
    {
        children_len = 0;
        for (int i = 0; i < population_size; i++)
        {
            parent_1 = tournament_selection(population, population_size, tournament_size);
            parent_2 = tournament_selection(population, population_size, tournament_size);

            random_number = ((double)rand() / (double)RAND_MAX);

            if (random_number < RECOMBINATION_RATE)
            {
                crossover(parent_1, parent_2, &child_1, &child_2, graph);
            }
            else
            {
                child_1 = parent_1;
                child_2 = parent_2;
            }

            children[children_len++] = child_1;
            children[children_len++] = child_2;

            crossover_children[0] = child_1;
            crossover_children[1] = child_2;

            int random_i = rand() % 2;
            selected_child = crossover_children[random_i];

            if (random_number < MUTATION_RATE)
            {
                selected_child = mutation(selected_child, graph);
            }

            new_route = swap_2_opt(selected_child->route, quota, graph);

            new_chromo = create_chromo(new_route, graph);

            children[children_len++] = new_chromo;
        }
        for (int i = 0; i < population_size; i++)
        {
            new_population[i] = population[i];
        }

        for (int i = 0; i < children_size; i++)
        {
            new_population[i + population_size] = children[i];
        }

        qsort(new_population, new_population_size, sizeof(Chromo **), compare_chromos);

        for (int i = 0; i < population_size; i++)
        {
            population[i] = new_population[i];
        }
    }

    route = population[0]->route;
    return route;
}

int *memetic_vns_algorithm(Graph *graph, int quota)
{
    int chromo_size = (graph->size * sizeof(int)) + sizeof(int) + sizeof(float);
    int population_size = POPULATION_SIZE;
    int tournament_size = TOURNAMENT_SIZE;
    int *route = (int *)malloc(graph->size * sizeof(int));
    int *new_route = (int *)malloc(graph->size * sizeof(int));
    double random_number;

    int children_len = 0;
    Chromo **new_population;
    Chromo **population;
    Chromo **children;
    Chromo **crossover_children;
    Chromo *parent_1;
    Chromo *parent_2;
    Chromo *child_1;
    Chromo *child_2;
    Chromo *new_chromo;
    Chromo *selected_child;

    population = (Chromo **)calloc(population_size, chromo_size);
    parent_1 = (Chromo *)malloc(chromo_size);
    parent_2 = (Chromo *)malloc(chromo_size);
    child_1 = (Chromo *)malloc(chromo_size);
    child_2 = (Chromo *)malloc(chromo_size);
    new_chromo = (Chromo *)malloc(chromo_size);
    selected_child = (Chromo *)malloc(chromo_size);
    population = init_population(graph, population_size, quota);

    int children_size = 3 * population_size;
    int new_population_size = children_size + population_size;

    children = (Chromo **)calloc(children_size, chromo_size);
    crossover_children = (Chromo **)calloc(2, chromo_size);
    new_population = (Chromo **)calloc(new_population_size, chromo_size);
    while (route_cost_calls < MAX_COST_CALLS)
    {
        children_len = 0;
        for (int i = 0; i < population_size; i++)
        {
            parent_1 = tournament_selection(population, population_size, tournament_size);
            parent_2 = tournament_selection(population, population_size, tournament_size);

            random_number = ((double)rand() / (double)RAND_MAX);

            if (random_number < RECOMBINATION_RATE)
            {
                crossover(parent_1, parent_2, &child_1, &child_2, graph);
            }
            else
            {
                child_1 = parent_1;
                child_2 = parent_2;
            }

            children[children_len++] = child_1;
            children[children_len++] = child_2;

            crossover_children[0] = child_1;
            crossover_children[1] = child_2;

            int random_i = rand() % 2;
            selected_child = crossover_children[random_i];

            if (random_number < MUTATION_RATE)
            {
                selected_child = mutation(selected_child, graph);
            }

            new_route = vns(selected_child->route, quota, graph);
            // printf("\nchildren_len: %d | route_cost_calls: %d | i: %d", children_len, route_cost_calls, i);

            new_chromo = create_chromo(new_route, graph);

            children[children_len++] = new_chromo;
        }
        for (int i = 0; i < population_size; i++)
        {
            new_population[i] = population[i];
        }

        for (int i = 0; i < children_size; i++)
        {
            new_population[i + population_size] = children[i];
        }

        qsort(new_population, new_population_size, sizeof(Chromo **), compare_chromos);
        for (int i = 0; i < population_size; i++)
        {
            population[i] = new_population[i];
        }
    }

    route = population[0]->route;
    return route;
}

int *memetic_vnd_algorithm(Graph *graph, int quota)
{
    int chromo_size = (graph->size * sizeof(int)) + sizeof(int) + sizeof(float);
    int population_size = POPULATION_SIZE;
    int tournament_size = TOURNAMENT_SIZE;
    int *route = (int *)malloc(graph->size * sizeof(int));
    int *new_route = (int *)malloc(graph->size * sizeof(int));
    double random_number;

    int children_len = 0;
    Chromo **new_population;
    Chromo **population;
    Chromo **children;
    Chromo **crossover_children;
    Chromo *parent_1;
    Chromo *parent_2;
    Chromo *child_1;
    Chromo *child_2;
    Chromo *new_chromo;
    Chromo *selected_child;

    population = (Chromo **)calloc(population_size, chromo_size);
    parent_1 = (Chromo *)malloc(chromo_size);
    parent_2 = (Chromo *)malloc(chromo_size);
    child_1 = (Chromo *)malloc(chromo_size);
    child_2 = (Chromo *)malloc(chromo_size);
    new_chromo = (Chromo *)malloc(chromo_size);
    selected_child = (Chromo *)malloc(chromo_size);
    population = init_population(graph, population_size, quota);

    int children_size = 3 * population_size;
    int new_population_size = children_size + population_size;

    children = (Chromo **)calloc(children_size, chromo_size);
    crossover_children = (Chromo **)calloc(2, chromo_size);
    new_population = (Chromo **)calloc(new_population_size, chromo_size);
    while (route_cost_calls < MAX_COST_CALLS)
    {
        children_len = 0;
        for (int i = 0; i < population_size; i++)
        {
            parent_1 = tournament_selection(population, population_size, tournament_size);
            parent_2 = tournament_selection(population, population_size, tournament_size);

            random_number = ((double)rand() / (double)RAND_MAX);

            if (random_number < RECOMBINATION_RATE)
            {
                crossover(parent_1, parent_2, &child_1, &child_2, graph);
            }
            else
            {
                child_1 = parent_1;
                child_2 = parent_2;
            }

            children[children_len++] = child_1;
            children[children_len++] = child_2;

            crossover_children[0] = child_1;
            crossover_children[1] = child_2;

            int random_i = rand() % 2;
            selected_child = crossover_children[random_i];

            if (random_number < MUTATION_RATE)
            {
                selected_child = mutation(selected_child, graph);
            }

            new_route = vnd(selected_child->route, quota, graph);
            // printf("\nchildren_len: %d | route_cost_calls: %d | i: %d", children_len, route_cost_calls, i);

            new_chromo = create_chromo(new_route, graph);

            children[children_len++] = new_chromo;
        }
        for (int i = 0; i < population_size; i++)
        {
            new_population[i] = population[i];
        }

        for (int i = 0; i < children_size; i++)
        {
            new_population[i + population_size] = children[i];
        }

        qsort(new_population, new_population_size, sizeof(Chromo **), compare_chromos);
        for (int i = 0; i < population_size; i++)
        {
            population[i] = new_population[i];
        }
    }

    route = population[0]->route;
    return route;
}

int *memetic_vns_vnd_algorithm(Graph *graph, int quota)
{
    int chromo_size = (graph->size * sizeof(int)) + sizeof(int) + sizeof(float);
    int population_size = POPULATION_SIZE;
    int tournament_size = TOURNAMENT_SIZE;
    int *route = (int *)malloc(graph->size * sizeof(int));
    int *new_route = (int *)malloc(graph->size * sizeof(int));
    double random_number;

    int children_len = 0;
    Chromo **new_population;
    Chromo **population;
    Chromo **children;
    Chromo **crossover_children;
    Chromo *parent_1;
    Chromo *parent_2;
    Chromo *child_1;
    Chromo *child_2;
    Chromo *new_chromo;
    Chromo *selected_child;

    population = (Chromo **)calloc(population_size, chromo_size);
    parent_1 = (Chromo *)malloc(chromo_size);
    parent_2 = (Chromo *)malloc(chromo_size);
    child_1 = (Chromo *)malloc(chromo_size);
    child_2 = (Chromo *)malloc(chromo_size);
    new_chromo = (Chromo *)malloc(chromo_size);
    selected_child = (Chromo *)malloc(chromo_size);
    population = init_population(graph, population_size, quota);

    int children_size = 3 * population_size;
    int new_population_size = children_size + population_size;

    children = (Chromo **)calloc(children_size, chromo_size);
    crossover_children = (Chromo **)calloc(children_size, chromo_size);
    new_population = (Chromo **)calloc(new_population_size, chromo_size);
    while (route_cost_calls < MAX_COST_CALLS)
    {
        children_len = 0;
        for (int i = 0; i < population_size; i++)
        {
            parent_1 = tournament_selection(population, population_size, tournament_size);
            parent_2 = tournament_selection(population, population_size, tournament_size);

            random_number = ((double)rand() / (double)RAND_MAX);

            if (random_number < RECOMBINATION_RATE)
            {
                crossover(parent_1, parent_2, &child_1, &child_2, graph);
            }
            else
            {
                child_1 = parent_1;
                child_2 = parent_2;
            }

            children[children_len++] = child_1;
            children[children_len++] = child_2;

            crossover_children[0] = child_1;
            crossover_children[1] = child_2;

            int random_i = rand() % 2;
            selected_child = crossover_children[random_i];

            if (random_number < MUTATION_RATE)
            {
                selected_child = mutation(selected_child, graph);
            }

            new_route = vns_vnd(selected_child->route, quota, graph);

            new_chromo = create_chromo(new_route, graph);

            children[children_len++] = new_chromo;
        }
        for (int i = 0; i < population_size; i++)
        {
            new_population[i] = population[i];
        }

        for (int i = 0; i < children_size; i++)
        {
            new_population[i + population_size] = children[i];
        }

        qsort(new_population, new_population_size, sizeof(Chromo **), compare_chromos);
        for (int i = 0; i < population_size; i++)
        {
            population[i] = new_population[i];
        }
    }

    route = population[0]->route;
    return route;
}

// GRASP

int *grasp_constructor(Graph *graph, int quota, float alpha_grasp)
{

    int tuple_size = 2 * sizeof(int);
    int *route = (int *)malloc(graph->size * sizeof(int));
    int i, j, id, random_i, insert_node;
    int best_economy, worst_economy;
    int economy_value;
    int economies_len, candidates_len, route_len = 0;

    float tsh;

    Edge edge;
    Edge edge1;
    Edge edge2;

    int *economy_tuple = (int *)malloc(tuple_size);
    int **economies = (int **)malloc(graph->size * sizeof(int *));
    int **candidates = (int **)malloc(graph->size * sizeof(int *));

    for (int i = 0; i < graph->size; i++)
    {
        economies[i] = (int *)malloc(tuple_size);
        candidates[i] = (int *)malloc(tuple_size);
    }

    int random_index = 1 + rand() % (graph->size - 1);

    route[route_len++] = graph->nodes[0].id;
    route[route_len++] = graph->nodes[random_index].id;
    route[route_len] = -1;

    int bonus_collected = calculate_bonus_collected(route, graph);

    while (bonus_collected < quota || best_economy > 0)
    {
        // RESET
        best_economy = INT_MIN;
        worst_economy = 0;
        economies_len = 0;
        candidates_len = 0;

        for (int i = 0; i < graph->size; i++)
        {
            economies[i] = (int *)malloc(tuple_size);
            candidates[i] = (int *)malloc(tuple_size);
        }

        for (int k = 0; k < graph->size; k++)
        {

            id = graph->nodes[k].id;
            if (!contains(route, route_len, id))
            {

                i = route[route_len - 2];
                j = route[route_len - 1];
                edge = graph->edges[i][j];
                edge1 = graph->edges[i][k];
                edge2 = graph->edges[k][j];
                economy_value = edge.weight + graph->nodes[k].penalty -
                                edge1.weight - edge2.weight;

                economy_tuple[0] = k;
                economy_tuple[1] = economy_value;
                economies[economies_len++] = economy_tuple;
            }
        }
        qsort(economies, economies_len, sizeof(int *), compare_tuples);
        best_economy = economies[0][1];
        worst_economy = economies[economies_len - 1][1];
        // tsh = best_economy - alpha_grasp * (best_economy - worst_economy);
        tsh = worst_economy + alpha_grasp * (best_economy - worst_economy);
        for (int i = 0; i < economies_len; i++)
        {
            if (economies[i][1] >= tsh)
            {
                candidates[candidates_len++] = economies[i];
            }
            else
            {
                break;
            }
        }
        random_i = rand() % candidates_len;
        if (best_economy > 0 || bonus_collected < quota)
        {
            insert_node = candidates[random_i][0];

            insert_on_index(route, route_len, route_len - 1, insert_node);
            route_len++;
            route[route_len] = -1;
        }
        bonus_collected = calculate_bonus_collected(route, graph);
    }
    return route;
}

int *grasp_algorithm(Graph *graph, int quota)
{
    float alpha_grasp = 0.2;
    int *best_route = (int *)malloc(graph->size * sizeof(int));
    int *route = (int *)malloc(graph->size * sizeof(int));
    int best_cost = INT_MAX;
    int cost, route_len;

    while (route_cost_calls < MAX_COST_CALLS)
    {
        route = grasp_constructor(graph, quota, alpha_grasp);
        route_len = get_route_length(route);

        route = swap_2_opt(route, quota, graph);

        cost = route_cost(route, graph, route_len);

        if (cost < best_cost)
        {
            best_cost = cost;
            best_route = route;
        }
    }

    return best_route;
}

int *grasp_vns_algorithm(Graph *graph, int quota)
{
    float alpha_grasp = 0.2;
    int *best_route = (int *)malloc(graph->size * sizeof(int));
    int *route = (int *)malloc(graph->size * sizeof(int));
    int best_cost = INT_MAX;
    int cost;

    while (route_cost_calls < MAX_COST_CALLS)
    {
        route = grasp_constructor(graph, quota, alpha_grasp);
        int route_len = get_route_length(route);

        route = vns(route, quota, graph);
        route_len = get_route_length(route);

        cost = route_cost(route, graph, route_len);

        if (cost < best_cost)
        {
            best_cost = cost;
            best_route = route;
        }
    }

    return best_route;
}

int *grasp_vnd_algorithm(Graph *graph, int quota)
{
    float alpha_grasp = 0.2;
    int *best_route = (int *)malloc(graph->size * sizeof(int));
    int *route = (int *)malloc(graph->size * sizeof(int));
    int best_cost = INT_MAX;
    int cost;

    while (route_cost_calls < MAX_COST_CALLS)
    {
        route = grasp_constructor(graph, quota, alpha_grasp);
        int route_len = get_route_length(route);

        route = vnd(route, quota, graph);
        route_len = get_route_length(route);

        cost = route_cost(route, graph, route_len);

        if (cost < best_cost)
        {
            best_cost = cost;
            best_route = route;
        }
    }
    return best_route;
}

int *grasp_vns_vnd_algorithm(Graph *graph, int quota)
{
    float alpha_grasp = 0.2;
    int *best_route = (int *)malloc(graph->size * sizeof(int));
    int *route = (int *)malloc(graph->size * sizeof(int));
    int best_cost = INT_MAX;
    int cost;

    while (route_cost_calls < MAX_COST_CALLS)
    {
        route = grasp_constructor(graph, quota, alpha_grasp);
        int route_len = get_route_length(route);

        // printf("\nNEW ROUTE: ");
        // for (int i = 0; i < route_len; i++)
        // {
        //     printf("%d -> ", route[i]);
        // }

        route = vns_vnd(route, quota, graph);

        route_len = get_route_length(route);

        // printf("\nAFTER vns_vnd ROUTE: ");
        // for (int i = 0; i < route_len; i++)
        // {
        //     printf("%d -> ", route[i]);
        // }

        cost = route_cost(route, graph, route_len);

        if (cost < best_cost)
        {
            best_cost = cost;
            best_route = route;
        }
    }

    return best_route;
}
