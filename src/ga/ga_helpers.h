#ifndef GA_HELPERS_H
#define GA_HELPERS_H


#include "aliases.h"


namespace ga
{


/*!
 * Allocate memory for a matrix
 * @param n     number of rows
 * @param m     number of columns
 * @return      pointer to the matrix
 */
matrix allocateMatrix(const int n,
                      const int m);

/*!
 * Copy values from one vector to another
 * @param arr1      vector to get the values from
 * @param arr2      pointer to vector to write the values to
 * @param left      left index (inclusive)
 * @param right     right index (exclusive)
 */
void copyVectorValues(const vector &arr1,
                            vector *arr2,
                      const int     left,
                      const int     right);

/*!
 * Copy values from one matrix to another
 * @param mat1      matrix to get the values from
 * @param mat2      pointer to matrix to write the values to
 * @param i1        left index for the rows (inclusive)
 * @param i2        right index for the rows (exclusive)
 * @param j1        left index for the columns (inclusive)
 * @param j2        right index for the columns (exclusive)
 */
void copyMatrixValues(const matrix &mat1,
                            matrix *mat2,
                      const int     i1,
                      const int     j1,
                      const int     i2,
                      const int     j2);

/*!
 * Find index of the maximum value of a vector
 * @param vec the list that we are finding the maximum of
 * @param len the length of that list
 * @return the index of the maximum element
 */
int findMaximumIndex(const vector  &vec,
                     const int      len);

/*!
 * Print a vector to console
 * @param vec   the vector to print
 */
void printVector(const vector &vec);

/*!
 * Print matrix to console
 * @param mat   the matrix to print
 */
void printMatrix(const matrix &mat);

/*!
 * Sum all entries in a vector
 * @param vec   the vector
 * @return      the sum
 */
double sumVector(const vector &vec);

/*!
 * Compute the mean value in a vector
 * @param vec       the vector
 * @param length    length of the vector
 * @return          the mean
 */
double vectorMEAN(const vector &vec,
                  const int     length);

/*!
 * Compute the standard deviation
 * @param vec    the vector of values
 * @param mean   the mean value
 * @param length the size of the vector
 * @return       the standard deviation
 */
double vectorSTD(const vector  &vec,
                 const double   mean,
                 const int      length);


} //namespace ga


#endif //GA_HELPERS_H
