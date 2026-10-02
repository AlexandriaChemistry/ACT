#include "sobol.h"

#define MAX_SOBOL_BITS 30
#define MAX_SOBOL_DIM 6

SobolSequence::SobolSequence()
{
    const int          degree[MAX_SOBOL_DIM]     = { 1, 2, 3, 3, 4, 4 };
    const unsigned int polynomial[MAX_SOBOL_DIM] = { 0, 1, 1, 2, 1, 4 };
    // Initiate the int_vec array
    int_vec.resize(MAX_SOBOL_BITS*MAX_SOBOL_DIM, 0);
    std::vector<unsigned int> int_vec_help = {
        1, 1, 1, 1, 1, 1, 3, 1, 3, 3, 1, 1, 5, 7, 7, 3, 3, 5, 15, 11, 5, 15, 13, 9
    };
    for(size_t i = 0; i < int_vec_help.size(); i++)
    {
        int_vec[i] = int_vec_help[i];
    }
    
    sobol_dim.resize(MAX_SOBOL_DIM, 0);

    // Initialize pointers to allow both 1D and 2D addressing.
    std::vector<unsigned int *> int_vec_ptr;
    int_vec_ptr.resize(MAX_SOBOL_BITS, nullptr);
    int k = 0;
    for (int j = 0; j < MAX_SOBOL_BITS; j++, k += MAX_SOBOL_DIM)
    {
        int_vec_ptr[j] = &(int_vec[k]);
    }
    for (int kdim = 0; kdim < MAX_SOBOL_DIM; kdim++)
    {
        // Update values in int_vec by shifting to the left
        for (int jdeg = 0; jdeg < degree[kdim]; jdeg++)
        {
            int_vec_ptr[jdeg][kdim] = int_vec_ptr[jdeg][kdim] << (MAX_SOBOL_BITS-1-jdeg);
        }
        
        // Stored values require normalization.
        for (int jbit = degree[kdim]; jbit < MAX_SOBOL_BITS; jbit++)
        {
            // Use recurrence to get other values.
            unsigned int this_poly = polynomial[kdim];
            unsigned int this_iv   = int_vec_ptr[jbit-degree[kdim]][kdim];
            this_iv                = this_iv ^ (this_iv >> degree[kdim]);
            for (int ldeg = degree[kdim]-1; ldeg >= 1; ldeg--)
            {
                if (this_poly & 1)
                {
                    this_iv = this_iv ^ int_vec_ptr[jbit-ldeg][kdim];
                }
                this_poly = this_poly * 2;
            }
            int_vec_ptr[jbit][kdim] = this_iv;
        }
    }
}

void SobolSequence::seq(int ndim, std::vector<double> *q)
{
    double       factor = 1.0/(1 << MAX_SOBOL_BITS);
    int          jbit;
    unsigned int this_index = index++;
    for(jbit = 0; jbit < MAX_SOBOL_BITS; jbit++)
    {
        if (!(this_index & 1))
        {
            break;
        }
        this_index = this_index / 2;
    }
    this_index = jbit * MAX_SOBOL_DIM;
    for(int kdim = 0; kdim < std::min(ndim, MAX_SOBOL_DIM); kdim++)
    {
        // Bitwise exclusive or.
        sobol_dim[kdim] = sobol_dim[kdim] ^ int_vec[this_index + kdim];
        // Since our integers are never more than 30 bits (MAX_SOBOL_BITS)
        // dividing by that (multiplying by factor) will give a float between 0 and 1.
        (*q)[kdim]  = sobol_dim[kdim]*factor;
    }
}

