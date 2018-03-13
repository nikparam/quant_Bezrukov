#pragma once

#include <cmath>
#include <iostream>
#include "quantumNumbers.hpp"
#include "MathUtils.hpp"

using std::cout;
using std::endl;

class Primitive
{
public:
	Primitive( double exponent, double coefficient, QuantumNumbers qNumbers ) :
		exponent(exponent), coefficient(coefficient), qNumbers(qNumbers)
	{
		renormalize();
	}

	void renormalize( )
	{
		int L = qNumbers.getAngularMomentum();
		double norm = std::pow( 2.0, 2 * L + 1.5 ) * std::pow( exponent, L + 1.5 ) / MathUtils::doubleFactorial(2 * qNumbers.i - 1) / MathUtils::doubleFactorial(2 * qNumbers.j - 1) / MathUtils::doubleFactorial(2 * qNumbers.k - 1) / pow(M_PI, 1.5);
	   	coefficient *= std::pow(norm, 0.5);	
	}

	double get_exponent() { return exponent; }
	double get_coefficient() { return coefficient; }
	QuantumNumbers get_qNumbers() { return qNumbers; }

	void show()
	{
		cout << "Primitive: (exponent) " << exponent << "; (coefficient): " << coefficient << "; (quantumNumbers) " << qNumbers.i << " " << qNumbers.j << " " << qNumbers.k << endl;
	}

	~Primitive() { }

private:
	double exponent;
	double coefficient;

	QuantumNumbers qNumbers;
};


