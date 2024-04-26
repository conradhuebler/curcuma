#include "hungarian.h"

void assign(int dim, double* distance, int* indices)
{

    Cell **t, **table, x;
    t = malloc(dim * sizeof(Cell*));
    table = malloc(dim * sizeof(Cell*));
    CellPosition* assignment;
    for (int i = 0; i < dim; i++) {
        t[i] = malloc(dim * sizeof(Cell));
        table[i] = malloc(dim * sizeof(Cell));
        for (int j = 0; j < dim; j++)
            table[i][j] = t[i][j] = distance[i * dim + j];
    }

    assignment = kuhn_match(dim, dim, table);

    int(*assigned)[dim][dim];
    assigned = calloc(1, sizeof(ssize_t[dim][dim]));

    if (assignment)
        for (int i = 0; i < dim; i++)
            (*assigned)[assignment[i].row][assignment[i].col] += 1;
    int index = 0;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            if ((*assigned)[i][j]) {
                indices[index] = j;
                index++;
            }
        }
    }

    free(assigned);
    free(t);
    free(table);
    free(assignment);
}
