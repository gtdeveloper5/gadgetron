/** \file       fatwater_test.cpp
    \brief      Test case for fat water seperation and its solvers

    \author     Hui Xue
*/

#include "hoNDArray.h"
#include "hoNDArray_math.h"
#include "simplexLagariaSolver.h"
#include "twoParaExpDecayOperator.h"
#include "curveFittingCostFunction.h"
#include "fatWaterMultiPeakR2StarOperator.h"
#include <gtest/gtest.h>
#include <boost/random.hpp>

using namespace Gadgetron;
using testing::Types;

template<typename T> class fatwater_test : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
    }
};

typedef Types< std::complex<double> > cxImplementations;
TYPED_TEST_CASE(fatwater_test, cxImplementations);

TYPED_TEST(fatwater_test, SignalModel)
{
    Gadgetron::fatWaterMultiPeakR2StarOperator< std::vector<TypeParam> > fw;

    std::vector<TypeParam> b(4);
    b[0] = 1.0;
    b[1] = 0.1;
    b[2] = 100;
    b[3] = 1000.0/300;

    std::vector<TypeParam> te(4);
    te[0] = 1.56;
    te[1] = 2.72;
    te[2] = 3.88;
    te[3] = 5.04;

    std::vector<TypeParam> y(4);
    fw.magnitude(te, b, y);

    EXPECT_NEAR(std::abs(y[0]), 0.00552042845193872, 0.00001);
    EXPECT_NEAR(std::abs(y[1]), 0.000116839405178112 , 0.00001);
}
