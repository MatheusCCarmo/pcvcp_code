#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "globals.h"
#include "graph.h"
#include "fitness.h"
#include "aux_functions.h"
#include "algorithms.c"

char *instance_path;
char *alg;

int file_exists(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (file != NULL)
    {
        fclose(file);
        return 1;
    }
    return 0;
}

void process_arguments(int argc, char *argv[])
{
    for (int i = 1; i < argc; i++)
    {
        if (strncmp(argv[i], "--recombination=", 16) == 0)
        {
            RECOMBINATION_RATE = atof(argv[i] + 16);
            // printf("\n | RECOMBINATION_RATE: %f", RECOMBINATION_RATE);
        }
        else if (strncmp(argv[i], "--mutation=", 11) == 0)
        {
            MUTATION_RATE = atof(argv[i] + 11);
            // printf("\n | MUTATION_RATE: %f", MUTATION_RATE);
        }
        else if (strncmp(argv[i], "--population_size=", 18) == 0)
        {
            POPULATION_SIZE = atoi(argv[i] + 18);
            // printf("\n | POPULATION_SIZE: %d", POPULATION_SIZE);
        }
        else if (strncmp(argv[i], "--tournament_size=", 18) == 0)
        {
            TOURNAMENT_SIZE = atoi(argv[i] + 18);
            // printf("\n | TOURNAMENT_SIZE: %d", TOURNAMENT_SIZE);
        }
        else if (strncmp(argv[i], "--instance=", 11) == 0)
        {
            instance_path = argv[i] + 11;
            // printf("\n | instance_path: %s", instance_path);
        }
        else if (strncmp(argv[i], "--alg=", 6) == 0)
        {
            alg = argv[i] + 6;
            // printf("\n | alg: %s", alg);
        }
    }
}

int main(int argc, char *argv[])
{

    srand(time(NULL));
    process_arguments(argc, argv);

    clock_t start_time, end_time;
    double duration;
    int cost;
    int number_of_nodes, n_passengers;

    char instance_full_path[256];

    snprintf(instance_full_path, sizeof(instance_full_path), "./instances/%s.in", instance_path);
    FILE *file = fopen(instance_full_path, "r");

    if (file == NULL)
    {
        printf("Nao foi possível abrir o arquivo.\n");
        return 1;
    }

    fscanf(file, "%d %d %*d\n", &number_of_nodes, &n_passengers);

    int **distance_matrix = (int **)calloc(number_of_nodes, sizeof(int *));

    for (int i = 0; i < number_of_nodes; i++)
    {
        distance_matrix[i] = (int *)calloc(number_of_nodes, sizeof(int));
    }

    // Ler matriz de custos
    for (int i = 0; i < number_of_nodes; i++)
    {
        for (int j = 0; j < number_of_nodes; j++)
        {
            fscanf(file, "%d", &distance_matrix[i][j]);
            // printf("\ndistance_matrix[%d][%d]: %d", i, j, distance_matrix[i][j]);
        }
    }

    // Pular linhas dos passageiros
    int ch;
    for (int i = 0; i < number_of_nodes + n_passengers + 3; i++)
    {
        do
        {
            ch = fgetc(file);
            if (ch == EOF)
                return 1;
        } while (ch != '\n');
    }

    int quota;
    fscanf(file, "%d\n", &quota);

    int **node_data = (int **)calloc(number_of_nodes, sizeof(int *));
    for (int i = 0; i < number_of_nodes; i++)
    {
        node_data[i] = (int *)calloc(3, sizeof(int)); // Vértice, bônus, delay
        fscanf(file, "%d %d %d", &node_data[i][0], &node_data[i][1], &node_data[i][2]);
    }
    node_data[0][2] = INT_MAX;

    fclose(file);

    Graph *graph = create_graph(number_of_nodes, distance_matrix, node_data);

    int *route = (int *)calloc(10, sizeof(int));

    char output_path[256];

    // INICIO ALGORITIMO
    if (strcmp(alg, "memetic") == 0)
    {
        snprintf(output_path, sizeof(output_path), "./output/%s/memetic_algorithm.txt", instance_path);
        start_time = clock();
        route = memetic_algorithm(graph, quota);
    }
    else if (strcmp(alg, "memetic_vns") == 0)
    {
        snprintf(output_path, sizeof(output_path), "./output/%s/memetic_algorithm_vns.txt", instance_path);
        start_time = clock();
        route = memetic_vns_algorithm(graph, quota);
    }
    else if (strcmp(alg, "memetic_vnd") == 0)
    {
        snprintf(output_path, sizeof(output_path), "./output/%s/memetic_algorithm_vnd.txt", instance_path);
        start_time = clock();
        route = memetic_vnd_algorithm(graph, quota);
    }
    else if (strcmp(alg, "memetic_vns_vnd") == 0)
    {
        snprintf(output_path, sizeof(output_path), "./output/%s/memetic_algorithm_vns_vnd.txt", instance_path);
        start_time = clock();
        route = memetic_vns_vnd_algorithm(graph, quota);
    }
    else if (strcmp(alg, "grasp") == 0)
    {
        snprintf(output_path, sizeof(output_path), "./output/%s/grasp.txt", instance_path);
        start_time = clock();
        route = grasp_algorithm(graph, quota);
    }
    else if (strcmp(alg, "grasp_vns") == 0)
    {
        snprintf(output_path, sizeof(output_path), "./output/%s/grasp_vns.txt", instance_path);
        start_time = clock();
        route = grasp_vns_algorithm(graph, quota);
    }
    else if (strcmp(alg, "grasp_vnd") == 0)
    {
        snprintf(output_path, sizeof(output_path), "./output/%s/grasp_vnd.txt", instance_path);
        start_time = clock();
        route = grasp_vnd_algorithm(graph, quota);
    }
    else if (strcmp(alg, "grasp_vns_vnd") == 0)
    {
        snprintf(output_path, sizeof(output_path), "./output/%s/grasp_vns_vnd.txt", instance_path);
        start_time = clock();
        route = grasp_vns_vnd_algorithm(graph, quota);
    }
    else
    {
        printf("Algoritmo desconhecido: %s\n", alg);
        return 1;
    }

    // FIM ALGORITIMO
    end_time = clock();
    duration = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    int route_len = get_route_length(route);

    int bonus = calculate_bonus_collected(route, graph);

    cost = route_cost(route, graph, route_len);

    FILE *file_output;
    if (file_exists(output_path))
    {

        file_output = fopen(output_path, "a");
    }
    else
    {

        file_output = fopen(output_path, "w");

        fprintf(file_output, "duration(s);cost;bonus;route\n");
    }

    fprintf(file_output, "%.2f;%d;%d;[", duration, cost, bonus);
    for (int j = 0; j < route_len; j++)
    {
        fprintf(file_output, "%d", route[j]);
        if (j != route_len - 1)
        {
            fprintf(file_output, ", ");
        }
    }
    fprintf(file_output, "]\n");

    fclose(file_output);

    // printf("\nroute_cost_calls: %d", route_cost_calls);
    // printf("\nneighboor_1_calls: %d", neighboor_1_calls);
    // printf("\nneighboor_2_calls: %d", neighboor_2_calls);
    // printf("\nneighboor_3_calls: %d", neighboor_3_calls);
    // printf("\nneighboor_4_calls: %d", neighboor_4_calls);
    // printf("\nneighboor_5_calls: %d", neighboor_5_calls);
    // printf("\nref_1_calls: %d", ref_1_calls);
    // printf("\nref_2_calls: %d", ref_2_calls);
    // printf("\nref_3_calls: %d", ref_3_calls);
    return 0;
}
