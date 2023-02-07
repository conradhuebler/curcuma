/*
 * <Munkress / Hungarian approach to solve assignment problems>
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

// Hungarian / Munkres' assignment algorithm, adopted as header only library taken from
// https://brc2.com/the-algorithm-workshop/
#pragma once

#include <Eigen/Dense>

#include <iostream>

typedef Eigen::MatrixXd Matrix;
static double epsilon = 1e-8;

void PrintMatrix(const Matrix& m)
{
    //  std::cout << m << std::endl;
}

void PrintVector(const std::vector<int>& vector)
{
    for (auto i : vector)
        std::cout << i << " ";
    std::cout << std::endl;
}

int Step1(Matrix& m, double& min2, const std::vector<int>& covered_cols, const std::vector<int>& covered_rows)
{
    min2 = m.maxCoeff();
    for (int i = 0; i < m.rows(); ++i) {
        if (covered_rows[i])
            continue;
        double min = m.row(i).minCoeff();
        for (int j = 0; j < m.cols(); ++j) {
            if (covered_cols[j])
                continue;
            m(i, j) -= min;
            if (m(i, j) > epsilon)
                min2 = std::min(min2, m(i, j));
        }
    }
    return 2;
}

int Step2(Matrix& starred, const Matrix& working, std::vector<int>& covered_cols, std::vector<int>& covered_rows)
{
    for (int r = 0; r < working.rows(); r++)
        for (int c = 0; c < working.cols(); c++) {
            if (working(r, c) == 0 && covered_rows[r] == 0 && covered_cols[c] == 0) {
                starred(r, c) = 1;
                covered_rows[r] = 1;
                covered_cols[c] = 1;
            }
        }
    for (int r = 0; r < covered_rows.size(); r++)
        covered_rows[r] = 0;
    for (int c = 0; c < covered_cols.size(); c++)
        covered_cols[c] = 0;
    return 3;
}

int Step3(const Matrix& starred, std::vector<int>& covered_cols)
{
    int colcount;
    for (int r = 0; r < starred.rows(); r++)
        for (int c = 0; c < starred.cols(); c++)
            if (starred(r, c) == 1)
                covered_cols[c] = 1;

    colcount = 0;
    for (int c = 0; c < covered_cols.size(); c++)
        if (covered_cols[c] == 1)
            colcount += 1;
    if (colcount >= starred.cols() || colcount >= starred.rows())
        return 7;
    else
        return 4;
}
void find_a_zero(const Matrix& working, int& row, int& col, const std::vector<int>& covered_rows, const std::vector<int>& covered_cols)
{
    int r = 0;
    int c;
    bool done;
    row = -1;
    col = -1;
    done = false;
    while (!done) {
        c = 0;
        while (true) {
            if (working(r, c) == 0 && covered_rows[r] == 0 && covered_cols[c] == 0) {
                row = r;
                col = c;
                done = true;
            }
            c += 1;
            if (c >= covered_cols.size() || done)
                break;
        }
        r += 1;
        if (r >= covered_rows.size())
            done = true;
    }
}

bool star_in_row(const Matrix& starred, int row)
{
    bool tmp = false;
    for (int c = 0; c < starred.cols(); c++)
        if (starred(row, c) == 1)
            tmp = true;
    return tmp;
}

void find_star_in_row(const Matrix& starred, int row, int& col)
{
    col = -1;
    for (int c = 0; c < starred.cols(); c++)
        if (starred(row, c) == 1)
            col = c;
}

int Step4(Matrix& starred, Matrix& working, int& col, int& row, std::vector<int>& covered_rows, std::vector<int>& covered_cols, int& path_row_0, int& path_col_0)
{
    // int row = -1;
    // int col = -1;
    bool done;

    done = false;
    while (!done) {
        find_a_zero(working, row, col, covered_rows, covered_cols);
        if (row == -1) {
            done = true;
            return 6;
        } else {
            starred(row, col) = 2;
            if (star_in_row(starred, row)) {
                find_star_in_row(starred, row, col);
                covered_rows[row] = 1;
                covered_cols[col] = 0;
            } else {
                done = true;
                path_row_0 = row;
                path_col_0 = col;
                return 5;
            }
        }
    }
    return 6;
}
void find_star_in_col(const Matrix& starred, int c, int& r)
{
    r = -1;
    for (int i = 0; i < starred.rows(); i++)
        if (starred(i, c) == 1)
            r = i;
}

void find_prime_in_row(const Matrix& starred, int r, int& c)
{
    for (int j = 0; j < starred.cols(); j++)
        if (starred(r, j) == 2)
            c = j;
}

void augment_path(Matrix& starred, const Matrix& path, int path_count)
{
    for (int p = 0; p < path_count; p++) {
        const int i = path(p, 0);
        const int j = path(p, 1);
        if (starred(i, j) == 1)
            starred(i, j) = 0;
        else
            starred(i, j) = 1;
    }
}

void clear_covers(std::vector<int>& covered_rows, std::vector<int>& covered_cols)
{
    for (int r = 0; r < covered_rows.size(); r++)
        covered_rows[r] = 0;
    for (int c = 0; c < covered_cols.size(); c++)
        covered_cols[c] = 0;
}

void erase_primes(Matrix& starred)
{
    for (int r = 0; r < starred.rows(); r++)
        for (int c = 0; c < starred.cols(); c++)
            if (starred(r, c) == 2)
                starred(r, c) = 0;
}

int Step5(Matrix& starred, int& path_count, Matrix& path, int& path_row_0, int& path_col_0, std::vector<int>& covered_rows, std::vector<int>& covered_cols)
{
    bool done;
    int r = -1;
    int c = -1;

    path_count = 1;
    path(path_count - 1, 0) = path_row_0;
    path(path_count - 1, 1) = path_col_0;
    done = false;
    while (!done) {
        find_star_in_col(starred, path(path_count - 1, 1), r);
        if (r > -1 && path_count < path.cols()) {
            path_count += 1;
            path(path_count - 1, 0) = r;
            path(path_count - 1, 1) = path(path_count - 2, 1);
        } else
            done = true;
        if (!done) {
            find_prime_in_row(starred, path(path_count - 1, 0), c);
            path_count += 1;
            path(path_count - 1, 0) = path(path_count - 2, 0);
            path(path_count - 1, 1) = c;
        }
    }
    augment_path(starred, path, path_count);
    clear_covers(covered_rows, covered_cols);
    erase_primes(starred);
    return 3;
}

double find_smallest(const Matrix& working, const std::vector<int>& covered_rows, const std::vector<int>& covered_cols)
{
    double minval = working.maxCoeff();
    for (int r = 0; r < covered_rows.size(); r++)
        for (int c = 0; c < covered_cols.size(); c++)
            if (covered_rows[r] == 0 && covered_cols[c] == 0)
                if (minval > working(r, c))
                    minval = working(r, c);
    return minval;
}
// Add the value found in Step 4 to every element of each covered row, and subtract
// it from every element of each uncovered column.  Return to Step 4 without
// altering any stars, primes, or covered lines.
int Step6(Matrix& working, const std::vector<int>& covered_rows, const std::vector<int>& covered_cols)
{
    double minval = find_smallest(working, covered_rows, covered_cols);
    for (int r = 0; r < covered_rows.size(); r++)
        for (int c = 0; c < covered_cols.size(); c++) {
            if (covered_rows[r] == 1)
                working(r, c) += minval;
            if (covered_cols[c] == 0)
                working(r, c) -= minval;
        }
    return 4;
}

Matrix MunkressAssign(const Matrix& m)
{
    int step = 1;
    double min2 = 0;
    Matrix working = m;
    Matrix starred = Matrix::Zero(m.cols(), m.rows());
    Matrix path = Matrix::Zero(m.cols(), m.rows());

    bool loop = true;
    std::vector<int> covered_rows(m.cols(), 0), covered_cols(m.rows(), 0);
    int path_count = 0;
    double min;
    int path_row_0 = 0;
    int path_col_0 = 0;

    while (loop) {
        int col = -1, row = -1;

        switch (step) {
        case 1:
            step = Step1(working, min2, covered_cols, covered_rows);
            break;
        case 2:
            step = Step2(starred, working, covered_cols, covered_rows);
            break;

        case 3:
            step = Step3(starred, covered_cols);
            break;

        case 4:
            step = Step4(starred, working, col, row, covered_rows, covered_cols, path_row_0, path_col_0);
            break;
        case 5:
            step = Step5(starred, path_count, path, path_row_0, path_col_0, covered_rows, covered_cols);
            break;

        case 6:
            step = Step6(working, covered_rows, covered_cols);
            break;
        case 7:
            loop = false;
            break;
        }
    }
    return starred;
}
