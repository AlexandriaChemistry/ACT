#ifndef ACT_HELPERS_H
#define ACT_HELPERS_H

double** allocateMatrix(const int n, const int m);
void freeMatrix(double** const mat, const int n);
double* allocateArray(const int n);
double rand01();
void copyArrayValues(double* const arr1, double* const arr2, const int length);

#endif //ACT_HELPERS_H
