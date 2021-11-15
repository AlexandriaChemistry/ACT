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


    /*!
     * Initializes parameters randomly between their lower and upper bound
     */
    class ACTRandomInitializer : public Initializer {

        std::random_device                      rd;
        std::mt19937                            gen;
        std::uniform_real_distribution<double>  dis;

        vector lb, ub;

    public:
        /*!
         * Create a new ACMRandomInitializer instance
         * @param lb    lower bound for each parameter
         * @param ub    upper bound for each parameter
         */
        ACMRandomInitializer(const vector& lb, const vector& ub)
        : gen(rd()), dis(std::uniform_real_distribution<>(0.0, 1.0)) {
            gen.seed(::time(NULL));
            this->lb = lb;
            this->ub = ub;
        }

        void initialize(      vector&   individual,
                        const int       length);

    };

}

#endif //ACT_INITIALIZER_H
