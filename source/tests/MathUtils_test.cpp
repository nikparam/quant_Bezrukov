#include "../MathUtils.hpp"

#include "gtest/gtest.h"

// Test DoubleFactorial function of negative numbers
TEST(DoubleFactorialTest, Negative)
{
	EXPECT_EQ(1, MathUtils::doubleFactorial(-5));
	EXPECT_EQ(1, MathUtils::doubleFactorial(-3));
	EXPECT_EQ(1, MathUtils::doubleFactorial(-1));
}

// Test DoubleFactorial function of Zero 
TEST(DoubleFactorialTest, Zero)
{
	EXPECT_EQ(1, MathUtils::doubleFactorial(0));
}

// Test DoubleFactorial function of postive even numbers 
TEST(DoubleFactorialTest, PositiveEven)
{
	EXPECT_EQ(2, MathUtils::doubleFactorial(2));
	EXPECT_EQ(8, MathUtils::doubleFactorial(4));
	EXPECT_EQ(48, MathUtils::doubleFactorial(6));
}

// Test DoubleFactorial function of postive odd numbers 
TEST(DoubleFactorialTest, PositiveOdd)
{
	EXPECT_EQ(1, MathUtils::doubleFactorial(1));
	EXPECT_EQ(3, MathUtils::doubleFactorial(3));
	EXPECT_EQ(15, MathUtils::doubleFactorial(5));
}
