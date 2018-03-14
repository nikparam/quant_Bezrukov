#pragma once

#include "QuantumNumbers.hpp"
#include "BasisFunction.hpp"
#include "Primitive.hpp"

#include <stdexcept> // std::invalid_argument

class BasisFunctionBuilder 
{
public:
	BasisFunctionBuilder( const char angularPart ); 
	~BasisFunctionBuilder() { }
	
	std::vector<QuantumNumbers> generateQuantumNumbers( const char angularPart );

	void add_primitive( double alpha, double coeff );

	std::vector<BasisFunction*> & getBasisFunctions() { return basisFunctions; }

private:
	char angularPart;
	std::vector<BasisFunction*> basisFunctions;
};
