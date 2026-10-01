#include "sobol.h"

#define MAX_SOBOL_BITS 30
#define MAX_SOBOL_DEGREE 6

SobolSequence::SobolSequence()
{
    const int          degree[MAX_SOBOL_DEGREE]     = { 1, 2, 3, 3, 4, 4 };
    const unsigned int polynomial[MAX_SOBOL_DEGREE] = { 0, 1, 1, 2, 1, 4 };
    // Initiate the int_vec array
    int_vec.resize(MAX_SOBOL_BITS*MAX_SOBOL_DEGREE,0);
    std::vector<unsigned int> int_vec_help = {
        1, 1, 1, 1, 1, 1, 3, 1, 3, 3, 1, 1, 5, 7, 7, 3, 3, 5, 15, 11, 5, 15, 13, 9
    };
    for(size_t i = 0; i < int_vec_help.size(); i++)
    {
        int_vec[i] = int_vec_help[i];
    }
    
    ix.resize(MAX_SOBOL_DEGREE, 0);

    // Initialize pointers to allow both 1D and 2D addressing.
    std::vector<unsigned int *> iu;
    iu.resize(MAX_SOBOL_BITS, nullptr);
    int k = 0;
    for (int j = 0; j < MAX_SOBOL_BITS; j++, k += MAX_SOBOL_DEGREE)
    {
        iu[j] = &int_vec[k];
    }
    // 
    for (int k = 0; k < MAX_SOBOL_DEGREE; k++)
    {
        for (int j = 0; j < degree[k]; j++)
        {
                iu[j][k] <<= (MAX_SOBOL_BITS-1-j);
        }
        
        // Stored values only require normalization.
        for (int j = degree[k]; j < MAX_SOBOL_BITS; j++)
        {
            // Use recurrence to get other values.
            unsigned int ipp  = polynomial[k];
            unsigned int i    = iu[j-degree[k]][k];
            i                ^= (i >> degree[k]);
            for (int l = degree[k]-1; l >= 1; l--)
            {
                if (ipp & 1)
                {
                    i ^= iu[j-l][k];
                }
                ipp >>= 1;
            }
            iu[j][k] = i;
        }
    }
}

void SobolSequence::seq(int n, std::vector<double> *x)
{
    double       fac = 1.0/(1 << MAX_SOBOL_BITS);
    int          j;
    unsigned int im = index++;
    for(j = 0; j < MAX_SOBOL_BITS; j++)
    {
        if (!(im & 1))
        {
            break;
        }
        im = im / 2;
    }
    im = j * MAX_SOBOL_DEGREE;
    for(int k = 0; k < std::min(n, MAX_SOBOL_DEGREE); k++)
    {
        // Bitwise exclusive or, yuck!
        ix[k] = ix[k] ^ int_vec[im + k];
        (*x)[k]  = ix[k]*fac;
    }
}

