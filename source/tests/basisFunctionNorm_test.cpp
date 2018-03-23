#include "basisFunctionNorm.hpp"

#include "../Basis.hpp"

#include "gtest/gtest.h"

TEST(checkNorm, UnitNorm)
{
	Basis basis;
	basis.read("../basis/cc-pvdz.gamess-us.dat");

	EXPECT_EQ(1, checkNorm(basis.getElement(0)->getBasisFunction(0)) ); 
}

