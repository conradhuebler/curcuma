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

/**
 * Cell markings
 **/
enum {
    UNMARKED = 0,
    MARKED,
    PRIME
};

/**
 *  Value type for marking
 */
typedef int_fast8_t Mark;

/**
 *  Value type for cells
 */
typedef double Cell;

typedef int_fast8_t Boolean;

typedef int_fast64_t BitSetLimb;

/**
 * Bit set, a set of fixed number of bits/booleans
 */
typedef struct {
    /**
     * The set of all limbs, a limb consist of 64 bits
     */
    BitSetLimb* limbs;

    /**
     * Singleton array with the index of the first non-zero limb
     */
    size_t first;

    /**
     * Array the the index of the previous non-zero limb for each limb
     */
    size_t* prev;

    /**
     * Array the the index of the next non-zero limb for each limb
     */
    size_t* next;

    char _buf[];
} BitSet;

typedef struct {
    size_t row;
    size_t col;
} CellPosition;

/**
 * Calculates the floored binary logarithm of a positive integer
 *
 * @param   value  The integer whose logarithm to calculate
 * @return         The floored binary logarithm of the integer
 */
#if defined(__GNUC__)
__attribute__((__const__))
#endif
static size_t
lb(BitSetLimb value);

/**
 * Constructor for BitSet
 *
 * @param   size  The (fixed) number of bits to bit set should contain
 * @return        The a unique BitSet instance with the specified size
 */
static BitSet*
bitset_create(size_t size);

/**
 * Gets the index of any set bit in a bit set
 *
 * @param   bitset  The bit set
 * @return        The index of any set bit
 */
#if defined(__GNUC__)
__attribute__((__pure__))
#endif
static ssize_t
bitset_any(BitSet* bitset);

/**
 * Turns off a bit in a bit set
 *
 * @param  bitset  The bit set
 * @param  i     The index of the bit to turn off
 */
static void
bitset_unset(BitSet* bitset, size_t i);

/**
 * Turns on a bit in a bit set
 *
 * @param  bitset  The bit set
 * @param  i     The index of the bit to turn on
 */
static void
bitset_set(BitSet* bitset, size_t i);

/**
 * Reduces the values on each rows so that, for each row, the
 * lowest cells value is zero, and all cells' values is decrease
 * with the same value [the minium value in the row].
 *
 * @param  n  The table's height
 * @param  m  The table's width
 * @param  t  The table in which to perform the reduction
 */
static void
kuhn_reduce_rows(size_t n, size_t m, Cell** t);

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
static Boolean
kuhn_is_done(size_t n, size_t m, Mark** marks, Boolean col_covered[m]);

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
static Mark**
kuhn_mark(size_t n, size_t m, Cell** t);

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
static Boolean
kuhn_find_prime(size_t n, size_t m, Cell** t, Mark** marks, Boolean row_covered[n], Boolean col_covered[m], CellPosition* primep);

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
static void
kuhn_alt_marks(size_t n, size_t m, Mark** marks, CellPosition alt[n * m],
    ssize_t col_marks[m], ssize_t row_primes[n], const CellPosition* prime);

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
static void
kuhn_add_and_subtract(size_t n, size_t m, Cell** t, Boolean row_covered[n], Boolean col_covered[m]);

/**
 * Creates a list of the assignment cells
 *
 * @param   n      The table's height
 * @param   m      The table's width
 * @param   marks  Matrix markings
 * @return         The assignment, an array of row–coloumn pairs
 */
static CellPosition*
kuhn_assign(size_t n, size_t m, Mark** marks);

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
extern CellPosition*
kuhn_match(size_t n, size_t m, Cell** table);
