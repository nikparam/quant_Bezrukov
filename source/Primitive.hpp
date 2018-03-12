#pragma once

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

class Primitive
{
public:
	Primitive( int number, double exponent, double coefficient ) :
		number(number), exponent(exponent), coefficient(coefficient)
		{
		}

	void renormalize( )
	{
		double N = std::pow( 2 * exponent / M_PI, 0.75 );
		coefficient *= N;
		//cout << number << " " << exponent << " " << coefficient << endl;
	}

	double get_coefficient()
	{
		return coefficient;
	}

	~Primitive() { }

private:
	int number;
	double exponent;
	double coefficient;
};


