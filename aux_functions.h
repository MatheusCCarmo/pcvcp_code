#include <stdlib.h>

#ifndef AUX_FUNCTIONS_H
#define AUX_FUNCTIONS_H

int contains(int *route, int route_len, int id)
{
    for (int i = 0; i < route_len; i++)
    {
        if (route[i] == id)
        {
            return 1;
        }
    }
    return 0;
}

void insert_on_index(int *array, int size, int index, int value)
{
    int temp = array[size - 1];
    for (int i = size - 1; i > index; i--)
    {
        array[i] = array[i - 1];
    }

    array[index] = value;
    array[size] = temp;
    array[size + 1] = -1;
}

void remove_by_value(int *arr, int size, int value)
{
    int i, j;
    for (i = 0, j = 0; i < size; ++i)
    {
        if (arr[i] != value)
        {
            arr[j++] = arr[i];
        }
    }

    arr[j] = -1;
}

void remove_by_index(int *arr, int size, int index)
{
    // printf("\nREMOVE_by_INDEX");
    // printf("| size: %d ", size);
    // printf("| index: %d ", index);
    int i;
    for (i = index; i < size - 1; ++i)
    {
        arr[i] = arr[i + 1];
    }
    arr[size - 1] = -1;
}

int get_route_length(int *route)
{
    int length = 1;
    // printf("\nROUTE:");
    while (route[length] != -1)
    {
        // printf(" %d -> ", route[length]);
        length++;
    }
    return length;
}

#endif
