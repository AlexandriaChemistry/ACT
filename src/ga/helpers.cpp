#include <stdlib.h>

#include "helpers.h"

#include "aliases.h"


/*!
 * Allocate memory for a matrix
 * @param n     number of rows
 * @param m     number of columns
 * @return      pointer to the matrix
 */
matrix allocateMatrix(const int n, const int m) {
    matrix mat(n);
    for (int i = 0; i < n; i++) mat[i] = vector(m)
    return mat;
}


/*!
 * Copy values from one vector to another
 * @param arr1      vector to get the values from
 * @param arr2      vector to write the values to
 * @param left      left index (inclusive)
 * @param right     right index (exclusive)
 */
void copyVectorValues(const vector arr1, const vector arr2, const int left, const int right) {
    for (int i = left; i < right; i++) {
        arr2[i] = arr1[i];
    }
}


/*!
 * Copy values from one matrix to another
 * @param mat1      matrix to get the values from
 * @param mat2      matrix to write the values to
 * @param i1        left index for the rows (inclusive)
 * @param i2        right index for the rows (exclusive)
 * @param j1        left index for the columns (inclusive)
 * @param j2        right index for the columns (exclusive)
 */
void copyMatrixValues(const matrix mat1, const matrix mat2, const int i1, const int i2, const int j1, const int j2) {
    int i, j;
    for (i = i1; i < i2; i++) {
        for (j = j1; j < j2; j++) {
            mat2[i][j] = mat1[i][j];
        }
    }
}
