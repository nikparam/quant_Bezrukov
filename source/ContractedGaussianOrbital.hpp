#pragma once

#include <vector>
#include <iostream>

#include "Primitive.hpp"
#include "QuantumNumbers.hpp"

class ContractedGaussianOrbital
{
public:
    ContractedGaussianOrbital( QuantumNumbers qNumbers ) : qNumbers(qNumbers)
	{
	}

    ~ContractedGaussianOrbital()
	{
	}

	void add_primitive( double alpha, double exponent );

	void show();
	
	size_t getPrimitivesCount() const { return primitives.size(); }
	QuantumNumbers & getQuantumNumbers() { return qNumbers; }

	Primitive & getPrimitive( int n ) { return primitives.at(n); }

private:
	vector <Primitive> primitives;
	QuantumNumbers qNumbers;
};

