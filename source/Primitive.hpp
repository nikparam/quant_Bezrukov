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

    void multiply_exponent( double mult ) { exponent *= mult; }
	void normalize( QuantumNumbers const & qNumbers );

    double get_exponent() const { return exponent; }
    double get_coefficient() const { return coefficient; }
    double get_norm( QuantumNumbers const& qNumbers ) const;

	void show();

private:
	double exponent;
	double coefficient;
};


