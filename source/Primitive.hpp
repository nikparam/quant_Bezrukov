#pragma once

#include <cmath>
#include <iostream>
#include "QuantumNumbers.hpp"
#include "MathUtils.hpp"

class Primitive
{
public:
	Primitive( double exponent, double coefficient ) :
		exponent(exponent), coefficient(coefficient) 
	{
	}
	
	~Primitive() { }

	void normalize( QuantumNumbers const & qNumbers );

	double get_exponent() { return exponent; }
	double get_coefficient() { return coefficient; }

	void show();

private:
	double exponent;
	double coefficient;
};


