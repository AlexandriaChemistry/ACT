#ifndef ACT_INITIALIZER_H
#define ACT_INITIALIZER_H

#include "aliases.h"

#include <time.h>
#import <random>


namespace ga {

    /*!
     * Abstract class for initializing individuals
     */
    class Initializer {

    public:
        /*!
         * Initialize an individual
         * @param individual    the individual to initialize
         * @param length        length of the chromosome
         */
        virtual void initialize(      vector&   individual,
                                const int       length) {};
    };


    /*!
     * Toy initializer. Initializes values randomly in range [min, max]
     */
    class SimpleInitializer : public Initializer {

        std::random_device                      rd;
        std::mt19937                            gen;
        std::uniform_real_distribution<double>  dis;

    public:
        /*!
         * Create a new SimpleInitializer object
         * @param min   minimum value to give to a gene
         * @param max   maximum value to give to a gene
         */
        SimpleInitializer(const double min,
                          const double max)
        : gen(rd()), dis(std::uniform_real_distribution<>(min, max)) {
            gen.seed(::time(NULL));
        }

        void initialize(      vector&   individual,
                        const int       length);

    };

}

#endif //ACT_INITIALIZER_H
