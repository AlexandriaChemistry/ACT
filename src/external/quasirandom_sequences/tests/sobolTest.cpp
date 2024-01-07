#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "../sobol.h"

#include <vector>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class SobolTest : public gmx::test::CommandLineTestBase
{
protected:
    void test08()
    {
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests I8_SOBOL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
# define DIM_MAX 4

            int dim_num = 3;
            long long int seed = 0;
            //long long int seed_in;
            //long long int seed_out;
            
            gmx::test::TestReferenceChecker checker_(this->rootChecker());
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
            
            checker_.checkInt64(dim_num, "dim_num");
            
            std::vector<double> values;
            for (int i = 0; i <= 11; i++ )
            {
                double r[DIM_MAX];
                //seed_in = seed;
                i8_sobol ( dim_num, &seed, r );
                //seed_out = seed;
                for(int j = 0; j < dim_num; j++)
                {
                    values.push_back(r[j]);
                }
            }
            checker_.checkSequence(values.begin(), values.end(), "Value");
        }
        return;
    }
# undef DIM_MAX

    //****************************************************************************80

    void test09()
    //****************************************************************************80
    //
    //  Purpose:
    //
    //    TEST09 tests I8_SOBOL.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license. 
    //
    //  Modified:
    //
    //    12 May 2007
    //
    //  Author:
    //
    //    John Burkardt
    //
    {
# define DIM_NUM 3

        double r[DIM_NUM];
        long long int seed = 0;
        //long long int seed_in;
        //long long int seed_out;
        
        gmx::test::TestReferenceChecker checker_(this->rootChecker());
        auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
        checker_.setDefaultTolerance(tolerance);
        
        std::vector<double> values;
        for (int i = 1; i <= 11; i++ )
        {
            //seed_in = seed;
            i8_sobol ( DIM_NUM, &seed, r );
            //seed_out = seed;
            
            for (int j = 0; j < DIM_NUM; j++ )
            {
                values.push_back(r[j]);
            }
        }
        
        seed = 100;
        
        for (int i = 1; i <= 5; i++ )
        {
            //seed_in = seed;
            i8_sobol ( DIM_NUM, &seed, r );
            //seed_out = seed;
            for (int j = 0; j < DIM_NUM; j++ )
            {
                values.push_back(r[j]);
            }
        }
        
        seed = 3;
        
        for (int i = 1; i <= 11; i++ )
        {
            //seed_in = seed;
            i8_sobol ( DIM_NUM, &seed, r );
            //seed_out = seed;
            for (int j = 0; j < DIM_NUM; j++ )
            {
                values.push_back(r[j]);
            }
        }
        checker_.checkSequence(values.begin(), values.end(), "Value");
        return;
    }
# undef DIM_NUM

};

TEST_F (SobolTest, Test08)
{
    test08();
}

TEST_F (SobolTest, Test09)
{
    test09();
}

} // namespace

} // namespace alexandria
