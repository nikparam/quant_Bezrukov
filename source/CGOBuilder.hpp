#pragma once

#include "QuantumNumbers.hpp"
#include "ContractedGaussianOrbital.hpp"
#include "Primitive.hpp"

#include <stdexcept> // std::invalid_argument

class CGOBuilder
{
public:
    CGOBuilder( const char angularPart );
    ~CGOBuilder() { }
	
	std::vector<QuantumNumbers> generateQuantumNumbers( const char angularPart );

    void add_primitive( double alpha, double coeff, double coeff2 = 0 );

    std::vector<ContractedGaussianOrbital*> & getCGOs() { return CGOs; }

private:
	char angularPart;
    std::vector<ContractedGaussianOrbital*> CGOs;
};
