/**
 * 𝓞(n³) implementation of the Hungarian algorithm
 *
 * Copyright (C) 2011, 2014, 2020  Mattias Andrée
 *
 * This program is free software. It comes without any warranty, to
 * the extent permitted by applicable law. You can redistribute it
 * and/or modify it under the terms of the Do What The Fuck You Want
 * To Public License, Version 2, as published by Sam Hocevar. See
 * http://sam.zoy.org/wtfpl/COPYING for more details.
 */

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "hungarian.h"

/**
 * Calculates the floored binary logarithm of a positive integer
 *
 * @param   value  The integer whose logarithm to calculate
 * @return         The floored binary logarithm of the integer
 */
#if defined(__GNUC__)
__attribute__((__const__))
#endif
size_t
lb(BitSetLimb value)
{
    size_t rc = 0;
    BitSetLimb v = value;

    if (v & (int_fast64_t)0xFFFFFFFF00000000LL) {
        rc |= 32L;
        v >>= 32;
    }
    if (v & (int_fast64_t)0x00000000FFFF0000LL) {
        rc |= 16L;
        v >>= 16;
    }
    if (v & (int_fast64_t)0x000000000000FF00LL) {
        rc |= 8L;
        v >>= 8;
    }
    if (v & (int_fast64_t)0x00000000000000F0LL) {
        rc |= 4L;
        v >>= 4;
    }
    if (v & (int_fast64_t)0x000000000000000CLL) {
        rc |= 2L;
        v >>= 2;
    }
    if (v & (int_fast64_t)0x0000000000000002LL) {
        rc |= 1L;
    }

    return rc;
}

/**
 * Constructor for BitSet
 *
 * @param   size  The (fixed) number of bits to bit set should contain
 * @return        The a unique BitSet instance with the specified size
 */
BitSet*
bitset_create(size_t size)
{
    size_t c = (size >> 6) + !!(size & 63L);
    BitSet* bitset = calloc(1, offsetof(BitSet, _buf) + c * sizeof(BitSetLimb) + 2 * (c + 1) * sizeof(size_t));

    bitset->limbs = (BitSetLimb*)&bitset->_buf[0];
    bitset->prev = (size_t*)&bitset->_buf[c * sizeof(BitSetLimb)];
    bitset->next = (size_t*)&bitset->_buf[c * sizeof(BitSetLimb) + c * sizeof(size_t)];

    return bitset;
}

/**
 * Gets the index of any set bit in a bit set
 *
 * @param   bitset  The bit set
 * @return        The index of any set bit
 */
#if defined(__GNUC__)
__attribute__((__pure__))
#endif
ssize_t
bitset_any(BitSet* bitset)
{
    size_t i;

    if (!bitset->first)
        return -1;

    i = bitset->first - 1;
    return (ssize_t)(lb(bitset->limbs[i] & -bitset->limbs[i]) + (i << 6));
}

/**
 * Turns off a bit in a bit set
 *
 * @param  bitset  The bit set
 * @param  i     The index of the bit to turn off
 */
void bitset_unset(BitSet* bitset, size_t i)
{
    size_t p, n, j = i >> 6;
    BitSetLimb old = bitset->limbs[j];

    bitset->limbs[j] &= ~(1LL << (i & 63L));

    if (!bitset->limbs[j] ^ !old) {
        j++;
        p = bitset->prev[j];
        n = bitset->next[j];
        bitset->prev[n] = p;
        bitset->next[p] = n;
        if (bitset->first == j)
            bitset->first = n;
    }
}

/**
 * Turns on a bit in a bit set
 *
 * @param  bitset  The bit set
 * @param  i     The index of the bit to turn on
 */
void bitset_set(BitSet* bitset, size_t i)
{
    size_t j = i >> 6;
    BitSetLimb old = bitset->limbs[j];

    bitset->limbs[j] |= 1LL << (i & 63L);

    if (!bitset->limbs[j] ^ !old) {
        j++;
        bitset->prev[bitset->first] = j;
        bitset->prev[j] = 0;
        bitset->next[j] = bitset->first;
        bitset->first = j;
    }
}

/**
 * Reduces the values on each rows so that, for each row, the
 * lowest cells value is zero, and all cells' values is decrease
 * with the same value [the minium value in the row].
 *
 * @param  n  The table's height
 * @param  m  The table's width
 * @param  t  The table in which to perform the reduction
 */
void kuhn_reduce_rows(size_t n, size_t m, Cell** t)
{
    size_t i, j;
    Cell min, *ti;

    for (i = 0; i < n; i++) {
        ti = t[i];
        min = *ti;
        for (j = 1; j < m; j++)
            if (min > ti[j])
                min = ti[j];
        for (j = 0; j < m; j++)
            ti[j] -= min;
    }
}

/**
 * Determines whether the marking is complete, that is
 * if each row has a marking which is on a unique column.
 *
 * @param   n            The table's height
 * @param   m            The table's width
 * @param   marks        The marking matrix
 * @param   col_covered  Column cover array
 * @return               Whether the marking is complete
 */
Boolean
kuhn_is_done(size_t n, size_t m, Mark** marks, Boolean col_covered[m])
{
    size_t i, j, count = 0;

    memset(col_covered, 0, m * sizeof(*col_covered));

    for (j = 0; j < m; j++) {
        for (i = 0; i < n; i++) {
            if (marks[i][j] == MARKED) {
                col_covered[j] = 1;
                break;
            }
        }
    }

    for (j = 0; j < m; j++)
        count += (size_t)col_covered[j];

    return count == n;
}

/**
 * Create a matrix with marking of cells in the table whose
 * value is zero [minimal for the row]. Each marking will
 * be on an unique row and an unique column.
 *
 * @param   n  The table's height
 * @param   m  The table's width
 * @param   t  The table in which to perform the reduction
 * @return     A matrix of markings as described in the summary
 */
Mark**
kuhn_mark(size_t n, size_t m, Cell** t)
{
    size_t i, j;
    Mark** marks;
    Boolean *row_covered, *col_covered;

    marks = malloc(n * sizeof(Mark*));
    for (i = 0; i < n; i++)
        marks[i] = calloc(m, sizeof(Mark)); /* UNMARKED == 0 */

    row_covered = calloc(n, sizeof(Boolean));
    col_covered = calloc(m, sizeof(Boolean));

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            if (!row_covered[i] && !col_covered[j] && !t[i][j]) {
                marks[i][j] = MARKED;
                row_covered[i] = 1;
                col_covered[j] = 1;
            }
        }
    }

    free(row_covered);
    free(col_covered);
    return marks;
}

/**
 * Finds a prime
 *
 * @param   n            The table's height
 * @param   m            The table's width
 * @param   t            The table
 * @param   marks        The marking matrix
 * @param   row_covered  Row cover array
 * @param   col_covered  Column cover array
 * @param   primep       Output parameter for the row and column of the found prime
 * @return               1 if a prime was found, 0 otherwise
 */
Boolean
kuhn_find_prime(size_t n, size_t m, Cell** t, Mark** marks, Boolean row_covered[n], Boolean col_covered[m], CellPosition* primep)
{
    size_t i, j, row, col;
    ssize_t p;
    Boolean mark_in_row;
    BitSet* zeroes = bitset_create(n * m);

    for (i = 0; i < n; i++)
        if (!row_covered[i])
            for (j = 0; j < m; j++)
                if (!col_covered[j] && !t[i][j])
                    bitset_set(zeroes, i * m + j);

    for (;;) {
        p = bitset_any(zeroes);
        if (p < 0) {
            free(zeroes);
            return 0;
        }

        row = (size_t)p / m;
        col = (size_t)p % m;

        marks[row][col] = PRIME;

        mark_in_row = 0;
        for (j = 0; j < m; j++) {
            if (marks[row][j] == MARKED) {
                mark_in_row = 1;
                col = j;
            }
        }

        if (mark_in_row) {
            row_covered[row] = 1;
            col_covered[col] = 0;

            for (i = 0; i < n; i++) {
                if (!t[i][col] && row != i) {
                    if (!row_covered[i] && !col_covered[col])
                        bitset_set(zeroes, i * m + col);
                    else
                        bitset_unset(zeroes, i * m + col);
                }
            }

            for (j = 0; j < m; j++) {
                if (!t[row][j] && col != j) {
                    if (!row_covered[row] && !col_covered[j])
                        bitset_set(zeroes, row * m + j);
                    else
                        bitset_unset(zeroes, row * m + j);
                }
            }

            if (!row_covered[row] && !col_covered[col])
                bitset_set(zeroes, row * m + col);
            else
                bitset_unset(zeroes, row * m + col);
        } else {
            free(zeroes);
            primep->row = row;
            primep->col = col;
            return 1;
        }
    }
}

/**
 * Removes all prime marks and modifies the marking
 *
 * @param  n           The table's height
 * @param  m           The table's width
 * @param  marks       The marking matrix
 * @param  alt         Marking modification paths
 * @param  col_marks   Markings in the columns
 * @param  row_primes  Primes in the rows
 * @param  prime       The last found prime
 */
void kuhn_alt_marks(size_t n, size_t m, Mark** marks, CellPosition alt[n * m],
    ssize_t col_marks[m], ssize_t row_primes[n], const CellPosition* prime)
{
    size_t i, j, index = 0;
    ssize_t row, col;
    Mark *markx, *marksi;

    alt[0].row = prime->row;
    alt[0].col = prime->col;

    for (i = 0; i < n; i++)
        row_primes[i] = -1;

    for (i = 0; i < m; i++)
        col_marks[i] = -1;

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            if (marks[i][j] == MARKED)
                col_marks[j] = (ssize_t)i;
            else if (marks[i][j] == PRIME)
                row_primes[i] = (ssize_t)j;
        }
    }

    while ((row = col_marks[alt[index].col]) >= 0) {
        index++;
        alt[index].row = (size_t)row;
        alt[index].col = alt[index - 1].col;

        col = row_primes[alt[index].row];
        index++;
        alt[index].row = alt[index - 1].row;
        alt[index].col = (size_t)col;
    }

    for (i = 0; i <= index; i++) {
        markx = &marks[alt[i].row][alt[i].col];
        *markx = *markx == MARKED ? UNMARKED : MARKED;
    }

    for (i = 0; i < n; i++) {
        marksi = marks[i];
        for (j = 0; j < m; j++)
            if (marksi[j] == PRIME)
                marksi[j] = UNMARKED;
    }
}

/**
 * Depending on whether the cells' rows and columns are covered,
 * the the minimum value in the table is added, subtracted or
 * neither from the cells.
 *
 * @param  n            The table's height
 * @param  m            The table's width
 * @param  t            The table to manipulate
 * @param  row_covered  Array that tell whether the rows are covered
 * @param  col_covered  Array that tell whether the columns are covered
 */
void kuhn_add_and_subtract(size_t n, size_t m, Cell** t, Boolean row_covered[n], Boolean col_covered[m])
{
    size_t i, j;
    Cell min = 0x7FFFFFFFL;

    for (i = 0; i < n; i++)
        if (!row_covered[i])
            for (j = 0; j < m; j++)
                if (!col_covered[j] && min > t[i][j])
                    min = t[i][j];

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            if (row_covered[i])
                t[i][j] += min;
            if (!col_covered[j])
                t[i][j] -= min;
        }
    }
}

/**
 * Creates a list of the assignment cells
 *
 * @param   n      The table's height
 * @param   m      The table's width
 * @param   marks  Matrix markings
 * @return         The assignment, an array of row–coloumn pairs
 */
CellPosition*
kuhn_assign(size_t n, size_t m, Mark** marks)
{
    CellPosition* assignment = malloc(n * sizeof(CellPosition));
    size_t i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            if (marks[i][j] == MARKED) {
                assignment[i].row = i;
                assignment[i].col = j;
            }
        }
    }

    return assignment;
}

/**
 * Calculates an optimal bipartite minimum weight matching using an
 * O(n³)-time implementation of The Hungarian Algorithm, also known
 * as Kuhn's Algorithm.
 *
 * @param   n      The height of the table
 * @param   m      The width of the table
 * @param   table  The table in which to perform the matching
 * @return         The optimal assignment, an array of row–coloumn pairs
 */

CellPosition*
kuhn_match(size_t n, size_t m, Cell** table)
{
    size_t i;
    ssize_t *row_primes, *col_marks;
    Mark** marks;
    Boolean *row_covered, *col_covered;
    CellPosition *ret, prime, *alt;

    /* Not copying table since it will only be used once. */

    row_covered = calloc(n, sizeof(Boolean));
    col_covered = calloc(m, sizeof(Boolean));

    row_primes = malloc(n * sizeof(ssize_t));
    col_marks = malloc(m * sizeof(ssize_t));

    alt = malloc(n * m * sizeof(CellPosition));

    kuhn_reduce_rows(n, m, table);
    marks = kuhn_mark(n, m, table);

    while (!kuhn_is_done(n, m, marks, col_covered)) {
        while (!kuhn_find_prime(n, m, table, marks, row_covered, col_covered, &prime))
            kuhn_add_and_subtract(n, m, table, row_covered, col_covered);
        kuhn_alt_marks(n, m, marks, alt, col_marks, row_primes, &prime);
        memset(row_covered, 0, n * sizeof(*row_covered));
        memset(col_covered, 0, m * sizeof(*col_covered));
    }

    free(row_covered);
    free(col_covered);
    free(alt);
    free(row_primes);
    free(col_marks);

    ret = kuhn_assign(n, m, marks);

    for (i = 0; i < n; i++)
        free(marks[i]);
    free(marks);

    return ret;
}
