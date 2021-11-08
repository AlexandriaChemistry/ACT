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
 * Generate a random number in range [0,1). srand() must be called before this method is used.
 * @return a double in [0,1)
 */
double rand01() {
    return (double) rand() / RAND_MAX;
}


/*!
 * Copy values from one vector to another
 * @param arr1      vector to get the values from
 * @param arr2      vector to write the values to
 * @param length    length of the vectors
 */
void copyArrayValues(const vector arr1, const vector arr2, const int length) {
    for (int i = 0; i < length; i++) {
        arr2[i] = arr1[i];
    }
}
