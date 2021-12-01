#ifndef ACT_GA_HELPERS_H
#define ACT_GA_HELPERS_H

#include "aliases.h"


namespace ga
{

    matrix allocateMatrix(const int n,
                          const int m);

    void copyVectorValues(const vector &arr1,
                                vector *arr2,
                          const int     left,
                          const int     right);

    void copyMatrixValues(const matrix &mat1,
                                matrix *mat2,
                          const int     i1,
                          const int     i2,
                          const int     j1,
                          const int     j2);

    int findMaximumIndex(const vector  &vec,
                         const int      len);

    void printVector(const vector &vec);

    void printMatrix(const matrix &mat);

    double sumVector(const vector &vec);

    double vectorMEAN(const vector &vec,
                      const int     length);

    double vectorSTD(const vector  &vec,
                     const double   mean,
                     const int      length);

}

#endif //ACT_GA_HELPERS_H
