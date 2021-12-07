#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#include "ga_helpers.h"

#include "aliases.h"


namespace ga
{


matrix allocateMatrix(const int n,
                        const int m)
{
    matrix mat(n);
    for (int i = 0; i < n; i++) mat[i] = vector(m);
    return mat;
}


void copyVectorValues(const vector &arr1,
                            vector *arr2,
                        const int     left,
                        const int     right)
{
    for (int i = left; i < right; i++)
    {
        (*arr2)[i] = arr1[i];
    }
}


void copyMatrixValues(const matrix &mat1,
                            matrix *mat2,
                        const int     i1,
                        const int     i2,
                        const int     j1,
                        const int     j2)
{
    int i, j;
    for (i = i1; i < i2; i++)
    {
        for (j = j1; j < j2; j++)
        {
            (*mat2)[i][j] = mat1[i][j];
        }
    }
}


int findMaximumIndex(const vector  &vec,
                        const int      len)
{
    int maxIndex = 0;
    for (int current = 0; current < len; current++)
        if (vec[current] > vec[maxIndex])
            maxIndex = current;
    return maxIndex;
}


void printVector(const vector &vec)
{
    printf("[ ");
    for (double value: vec) printf("%f ", value);
    printf("]\n");
}


void printMatrix(const matrix &mat)
{
    for (vector vec: mat) printVector(vec);
}


double sumVector(const vector &vec)
{
    double sum = 0;
    for (double ele : vec) sum += ele;
    return sum;
}


double vectorMEAN(const vector &vec,
                  const int     length)
{
    return sumVector(vec) / length;
}


double vectorSTD(const vector  &vec,
                 const double   mean,
                 const int      length)
{
    double sum = 0;
    for (double ele : vec) sum += (ele - mean)*(ele - mean);
    return sqrt(sum / length);
}


} //namespace ga
