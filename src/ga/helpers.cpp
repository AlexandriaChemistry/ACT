#include <stdlib.h>

#include "helpers.h"


/*!
 * Allocate memory for a matrix
 * @param n     number of rows
 * @param m     number of columns
 * @return      pointer to the matrix
 */
double** allocateMatrix(const int n, const int m) {
    double** mat = (double**) malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        mat[i] = (double*) malloc(m * sizeof(double));
    }
    return mat;
}


/*!
 * Free memory occupied by a matrix
 * @param mat   the matrix
 * @param n     number of rows
 */
void freeMatrix(double** const mat, const int n) {
    for (int i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
}


/*!
 * Allocate memory for an array
 * @param n     length of the array
 * @return      pointer to the array
 */
double* allocateArray(const int n) {
    return arr = (double*) malloc(n * sizeof(double))
}


/*!
 * Generate a random number in range [0,1). srand() must be called before this method is used.
 * @return a double in [0,1)
 */
double rand01() {
    return (double) rand() / RAND_MAX;
}


/*!
 * Copy values from one array to another
 * @param arr1      array to get the values from
 * @param arr2      array to write the values to
 * @param length    length of the arrays
 */
void copyArrayValues(double* const arr1, double* const arr2, const int length) {
    for (int i = 0; i < length; i++) {
        arr2[i] = arr1[i];
    }
}
